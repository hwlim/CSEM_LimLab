#ifndef CHROMTABLE_H_
#define CHROMTABLE_H_

#include<cstdio>
#include<cassert>
#include<string>
#include<vector>
#include<fstream>
#include<algorithm>
#include<pthread.h>

#include "utils.h"
#include "my_assert.h"

#include "ChrMap.h"
#include "Alignment.h"
#include "Chromosome.h"
#include "PThreadWrapper.h"

class ChromTable {
 public:
  ChromTable(ChrMap*, std::vector<Alignment>&, int, int, const char*);
  ~ChromTable();

  void update(bool);

  double getMaxDelta() { return max_delta; }

 private:
  CHR_ID_TYPE m;
  HIT_INT_TYPE nAmts;
  int halfws, nThreads;
  double max_delta;

  bool updateFracs; // if update the frac fields of the alignments vector

  ChrMap* chrMap;

  std::vector<Alignment> &alignments;
  std::vector<Chromosome*> chroms_multi; 

  struct Params {
    int no;
    ChromTable *pointer;
    std::vector<CHR_ID_TYPE> chroms;
    
    Params(int no, ChromTable *pointer) { this->no = no; this->pointer = pointer; chroms.clear(); }
  };

  std::vector<Params> paramsArray;
  PThreadWrapper pthreadWrapper;

  void loadPrior(const char*);
  void assign_chromosomes_to_threads();
  void update_per_thread(Params*);

  static void* update_per_thread_wrapper(void* args) {
    Params *params = (Params*)args;
    params->pointer->update_per_thread(params);
    return NULL;
  }
};

ChromTable::ChromTable(ChrMap* chrMap, std::vector<Alignment>& alignments, int halfws, int nThreads, const char* priorF) : halfws(halfws), nThreads(nThreads), chrMap(chrMap), alignments(alignments) {

  m = chrMap->size();
  nAmts = alignments.size();

  max_delta = 0.0;
  updateFracs = true;

  // initialize chroms_multi
  for (CHR_ID_TYPE i = 0; i < m; i++) chroms_multi.push_back(new Chromosome(halfws, chrMap->getLen(i), alignments));
  for (HIT_INT_TYPE i = 0; i < nAmts; i++) chroms_multi[alignments[i].cid]->addPos(i, alignments[i].isMulti);
  for (CHR_ID_TYPE i = 0; i < m; i++) chroms_multi[i]->init();

  printf("Discretization is performed!\n");

  if (priorF[0] != 0) loadPrior(priorF);

  assign_chromosomes_to_threads();

  printf("ChromTable is constructed!\n");
}

ChromTable::~ChromTable() {
  for (CHR_ID_TYPE i = 0; i < m; i++) delete chroms_multi[i];
}

void ChromTable::loadPrior(const char* priorF) {
  int ngroups; // number of groups
  std::vector<double> groupvalues;
  std::string cname, line;

  std::ifstream fin(priorF);

  fin>>ngroups;
  groupvalues.assign(ngroups, 0.0);
  for (int i = 0; i < ngroups; i++) {
    fin>>groupvalues[i];
    --groupvalues[i]; // deduct one
  }

  while (fin>>cname) {
    getline(fin, line);
    CHR_ID_TYPE cid = chrMap->getCid(cname);
    chroms_multi[cid]->processPriorInfo(line, ngroups, groupvalues);
  }

  fin.close();
  printf("Prior information are loaded and processed!\n");
}

void ChromTable::assign_chromosomes_to_threads() {
  HIT_INT_TYPE avghits, curhits;
  int curthread;

  avghits = 0;
  for (CHR_ID_TYPE i = 0; i < m; i++) avghits += chroms_multi[i]->getSize();
  avghits /= nThreads;

  paramsArray.clear();
  curhits = 0; curthread = 0;
  paramsArray.push_back(Params(curthread, this));
  for (CHR_ID_TYPE i = 0; i < m; i++) {
    paramsArray[curthread].chroms.push_back(i);
    curhits += chroms_multi[i]->getSize();

    if (curthread + 1 < nThreads && i + 1 < m && (curhits >= avghits || nThreads - curthread - 1 >= m - i - 1)) {
      curthread++; curhits = 0;
      paramsArray.push_back(Params(curthread, this));
    }
  }

  pthreadWrapper.num_threads_used = paramsArray.size();
  pthreadWrapper.threads.assign(pthreadWrapper.num_threads_used, pthread_t());

  printf("Jobs are assigned!\n");
}

void ChromTable::update_per_thread(Params* params) {
  for (size_t i = 0; i < params->chroms.size(); i++) {
    CHR_ID_TYPE chrom_id = params->chroms[i];
    chroms_multi[chrom_id]->update(updateFracs);
  }
}

//multi-threading
void ChromTable::update(bool updateFracs = true) {
  this->updateFracs = updateFracs; 

  // create threads
  for (int i = 0; i < pthreadWrapper.num_threads_used; i++) {
    pthreadWrapper.rc = pthread_create(&pthreadWrapper.threads[i], &pthreadWrapper.attr, update_per_thread_wrapper, (void*)(&paramsArray[i]));
    pthread_assert(pthreadWrapper.rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) when ChromTable is being updated!");
  }
  // join threads
  for (int i = 0; i < pthreadWrapper.num_threads_used; i++) {
    pthreadWrapper.rc = pthread_join(pthreadWrapper.threads[i], NULL);
    pthread_assert(pthreadWrapper.rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) when ChromTable is being updated!");
  }

  // update max_delta
  max_delta = 0.0;
  for (CHR_ID_TYPE i = 0; i < m; i++) max_delta = std::max(max_delta, chroms_multi[i]->getMaxDelta());
}

#endif
