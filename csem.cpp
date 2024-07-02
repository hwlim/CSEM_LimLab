#include<ctime>
#include<cmath>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<pthread.h>

#include "utils.h"
#include "my_assert.h"

#include "ChrMap.h"
#include "BamAlignment.h"
#include "SamParser.h"
#include "BamWriter.h"
#include "Alignment.h"
#include "ChromTable.h"
#include "PThreadWrapper.h"

using namespace std;

struct Params {
  int no;
  vector<READ_INT_TYPE> reads;

  Params(int no) { this->no = no; reads.clear(); }
};

bool extendReads;

int ROUND = 0;
int UPPERBOUND = 200; // 200 by default

int nThreads = 1; // 1 by default

int fragment_length, halfws;

CHR_ID_TYPE m; // m chromosomes
READ_INT_TYPE n; // n reads
HIT_INT_TYPE nAmts; // nAmts tot # of alignments

ChrMap *chrMap;

SamParser *samParser;
BamWriter *bamWriter;

char inpType;
char inpF[STRLEN], outName[STRLEN];
char priorF[STRLEN];

vector<HIT_INT_TYPE> s;
vector<Alignment> alignments;

ChromTable *chromTable;

// for multi-threading
vector<Params> paramsArray;
PThreadWrapper pthreadWrapper;

READ_INT_TYPE nUniqe, nMulti;

void loadData() {
  string currentReadName, readName;
  BamAlignment b;

  HIT_INT_TYPE cnt = 0;

  samParser = new SamParser(inpType, inpF);

  s.clear();
  alignments.clear();
  currentReadName = "";
  while (samParser->next(b)) {
    ++cnt;
    if (cnt % 1000000 == 0) fprintf(stderr, "%u FIN\n", cnt);

    if (!b.isAligned()) continue;
    readName = b.getName();
    if (currentReadName != readName) {
      // push the start position into s
      s.push_back(alignments.size());
      currentReadName = readName;
    }
    alignments.push_back(Alignment(b.getCid(), (extendReads? b.getMidPos(fragment_length) : b.getPos()), b.getDir()));  // extend reads or not
  }

  chrMap = new ChrMap(samParser->getHeader());
  
  m = chrMap->size();
  n = s.size(); 
  nAmts = alignments.size();
  s.push_back(nAmts);

  delete samParser;

  fprintf(stderr, "Loading data is finished!\n");
}

void splitJobs_and_Init() {
  int cur_thread;
  double tot;
  CHR_LEN_TYPE lb, ub;

  cur_thread = 0;
  paramsArray.clear();
  for (READ_INT_TYPE i = 0; i < n; i++) {
    if (s[i + 1] - s[i] == 1) { alignments[s[i]].frac = 1.0; continue; }

    // initialization, for each multi-read, distribute it uniformly
    tot = 0.0;
    for (HIT_INT_TYPE j = s[i]; j < s[i + 1]; j++) {
      lb = max(-1, alignments[j].pos - halfws - 1);
      ub = min(chrMap->getLen(alignments[j].cid) - 1, alignments[j].pos + halfws);
      alignments[j].frac = ub - lb;
      tot += alignments[j].frac;
    }

    for (HIT_INT_TYPE j = s[i]; j < s[i + 1]; j++) {
      alignments[j].frac /= tot;
      alignments[j].isMulti = true;
    }

    // assigning reads to threads
    if (pthreadWrapper.num_threads_used < nThreads) { 
      paramsArray.push_back(Params(pthreadWrapper.num_threads_used));
      ++pthreadWrapper.num_threads_used;
    }
    paramsArray[cur_thread].reads.push_back(i);
    ++cur_thread; if (cur_thread == nThreads) cur_thread = 0;
  }

  pthreadWrapper.threads.assign(pthreadWrapper.num_threads_used, pthread_t());

  chromTable = new ChromTable(chrMap, alignments, halfws, nThreads, priorF);
  fprintf(stderr, "Splitting jobs and initialization are finished!\n");
}

void* allocateMultiReads_per_thread(void* arg) {
  Params *params = (Params*)arg;
  double tot;

  for (size_t i = 0; i < params->reads.size(); i++) {
    READ_INT_TYPE rid = params->reads[i];

    tot = 0.0;
    for (HIT_INT_TYPE j = s[rid]; j < s[rid + 1]; j++) {
      assert(alignments[j].frac >= 0.0);
      tot += alignments[j].frac;
    }

    if (tot <= 0.0) tot = s[rid + 1] - s[rid]; // if adding prior leads to all fracs be 0, allocate the read uniformly

    for (HIT_INT_TYPE j = s[rid]; j < s[rid + 1]; j++) 
      alignments[j].frac /= tot;
  
  }

  return NULL;
}

void allocateMultiReads() {
  // update chromTable
  chromTable->update(UPPERBOUND > 0);

  for (ROUND = 1; ROUND <= UPPERBOUND; ROUND++) {
    // allocate muti-reads    
    // create threads
    for (int i = 0; i < pthreadWrapper.num_threads_used; i++) {
      pthreadWrapper.rc = pthread_create(&pthreadWrapper.threads[i], &pthreadWrapper.attr, allocateMultiReads_per_thread, (void*)(&paramsArray[i]));
      pthread_assert(pthreadWrapper.rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) at ROUND " + itos(ROUND) + "!");
    }
    // join threads
    for (int i = 0; i < pthreadWrapper.num_threads_used; i++) {
      pthreadWrapper.rc = pthread_join(pthreadWrapper.threads[i], NULL);
      pthread_assert(pthreadWrapper.rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) at ROUND " + itos(ROUND) + "!");
    }

    // update chromTable
    chromTable->update(ROUND < UPPERBOUND);

    fprintf(stderr, "ROUND = %d, MAX_DELTA = %.6g\n", ROUND, chromTable->getMaxDelta());
  }
}

void output() {
  HIT_INT_TYPE p;
  BamAlignment b;

  char outF[STRLEN];

  sprintf(outF, "%s.bam", outName);

  samParser = new SamParser(inpType, inpF);
  bamWriter = new BamWriter(outF, samParser->getHeader());

  HIT_INT_TYPE cnt = 0;

  p = 0;
  while (samParser->next(b)) {
    if (b.isAligned()) {
      b.setFrac(alignments[p++].frac);
    }
    bamWriter->write(b);

    ++cnt;
    if (cnt % 1000000 == 0) fprintf(stderr, "%u FIN\n", cnt);
  }

  delete samParser;
  delete bamWriter;

  fprintf(stderr, "Writing output is finished!\n");
}

int main(int argc, char* argv[]) {
  if (argc < 7 || argc > 10) {
    fprintf(stderr, "Usage : csem input_type input_file fragment_length UPPERBOUND output_name number_of_threads [--extend-reads] [--prior prior_file]\n");
    exit(-1);
  }

  assert(strlen(argv[1]) == 1);

  inpType = argv[1][0];
  strcpy(inpF, argv[2]);
  fragment_length = atoi(argv[3]); 
  UPPERBOUND = atoi(argv[4]);
  strcpy(outName, argv[5]);
  nThreads = atoi(argv[6]);

  extendReads = false;
  priorF[0] = 0;

  for (int i = 7; i < argc; i++) {
    if (!strcmp(argv[i], "--extend-reads")) { extendReads = true; }
    if (!strcmp(argv[i], "--prior")) { assert(strlen(argv[i + 1]) > 0); strcpy(priorF, argv[i + 1]); }
  }
 
  halfws = fragment_length / 2;

  chrMap = NULL;
  chromTable = NULL;

  loadData();
  splitJobs_and_Init();
  allocateMultiReads();

  delete chrMap;
  delete chromTable;
  
  output();

  return 0;
}
