// SamParser for csem

#ifndef SAMPARSER_H_
#define SAMPARSER_H_

#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<vector>

#include "sam/bam.h"
#include "sam/sam.h"

#include "utils.h"
#include "my_assert.h"

#include "BamAlignment.h"

class SamParser {
 public:
  SamParser(char, const char*, const char* = 0);
  ~SamParser();

  const bam_header_t* getHeader() const { 
    return header;
  }

  bool next(BamAlignment& b) { return b.read(sam_in); }

 private:
  samfile_t *sam_in;
  bam_header_t *header;
};

// aux, if not 0, points to the file name of fn_list
SamParser::SamParser(char inpType, const char* inpF, const char* aux) {
  switch(inpType) {
  case 'b': sam_in = samopen(inpF, "rb", aux); break;
  case 's': sam_in = samopen(inpF, "r", aux); break;
  default: assert(false);
  }

  general_assert(sam_in != 0, "Cannot open " + cstrtos(inpF) + "! It may not exist.");
  header = sam_in->header;
  general_assert(header != 0, "Fail to parse the header!");
}

SamParser::~SamParser() {
  samclose(sam_in);
}

#endif

