#include<ctime>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cassert>
#include<string>
#include<vector>

#include "utils.h"
#include "my_assert.h"

#include "sampling.h"
#include "BamAlignment.h"
#include "SamParser.h"
#include "BamWriter.h"

using namespace std;

int choice, outputFormat;

char outF[STRLEN];

SamParser *samParser;
BamWriter *bamWriter;
FILE *fo;

vector<BamAlignment> alignments;

uniform01 *rg;

// If output format is bam and choice is not to keep original, ZW tag will be removed but the QUAL field is not changed
void printOut(BamAlignment& b) {
  switch(outputFormat) {
  case 0: if (choice != 0) b.removeZWTag(); bamWriter->write(b); break;
  case 1: b.writeToBED(fo, samParser->getHeader()); break;
  case 2: b.writeToTagAlign(fo, samParser->getHeader()); break;
  default: assert(false);
  }
}

void process_a_read() {
  vector<double> arr;
  size_t id;

  switch(choice) {
  case 1 : 
    if (alignments.size() == 1) printOut(alignments[0]);
    break;
  case 2 : 
    arr.assign(alignments.size(), 0.0);
    for (size_t i = 0; i < alignments.size(); i++) {
      arr[i] = alignments[i].getFrac();
      if (i > 0) arr[i] += arr[i - 1];
    }
    id = sample(*rg, arr, (int)alignments.size());
    printOut(alignments[id]);
    break;
  default : assert(false);
  }
}

int main(int argc, char* argv[]) {
  if (argc != 5) {
    printf("Usage: csem-bam-processor input.bam output_name <keep orignal bam 0; unique only 1; sampling 2> <bam 0; bed 1; tagAlign 2>\n");
    exit(-1);
  }

  choice = atoi(argv[3]); 
  assert(choice >= 0 && choice <= 2);
  outputFormat = atoi(argv[4]);
  assert(outputFormat >= 0 && outputFormat <= 2);

  string currentReadName, readName;
  BamAlignment b;

  HIT_INT_TYPE cnt = 0;

  samParser = new SamParser('b', argv[1]);
  bamWriter = NULL; fo = NULL;
  rg = NULL;

  if (choice == 2) rg = new uniform01(engine_type(time(NULL)));

  // decide output format
  switch(outputFormat) {
  case 0 : 
    sprintf(outF, "%s.bam", argv[2]);
    bamWriter = new BamWriter(outF, samParser->getHeader());
    break;
  case 1 :
    sprintf(outF, "%s.bed", argv[2]);
    fo = fopen(outF, "w");
    break;
  case 2 :
    sprintf(outF, "%s.tagAlign", argv[2]);
    fo = fopen(outF, "w");
    break;
  default : assert(false);
  }

  alignments.clear();
  currentReadName = "";

  while (samParser->next(b)) {
    ++cnt;
    if (cnt % 1000000 == 0) printf("%u FIN\n", cnt);

    if (!b.isAligned()) continue;

    if (choice == 0) { printOut(b); continue; }

    readName = b.getName();
    if (currentReadName != readName) {
      if (currentReadName != "") process_a_read();
      currentReadName = readName;
      alignments.clear();
    }
    alignments.push_back(b);
  }
  if (currentReadName != "") process_a_read();

  delete samParser;
  if (choice == 2) delete rg;
  if (outputFormat == 0) delete bamWriter;
  else fclose(fo);
  
  return 0;
}
