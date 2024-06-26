// BamWriter for csem

#ifndef BAMWRITER_H_
#define BAmWRITER_H_

#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<vector>
#include<sstream>

#include "sam/bam.h"
#include "sam/sam.h"

#include "utils.h"
#include "my_assert.h"

#include "sam_csem_aux.h"

#include "BamAlignment.h"

class BamWriter {
 public:
  BamWriter(const char*, const bam_header_t*) ;
  ~BamWriter();
  
  void write(BamAlignment& b) { b.write(bam_out); }
  
 private:
  samfile_t *bam_out;
};

BamWriter::BamWriter(const char* outF, const bam_header_t* header) {
  bam_header_t *out_header = bam_header_dwt(header);
  
  std::ostringstream strout;
  strout<<"@HD\tVN:1.4\tSO:unknown\n@PG\tID:CSEM\n";
  std::string content = strout.str();
  append_header_text(out_header, content.c_str(), content.length());

  bam_out = samopen(outF, "wb", out_header);
  general_assert(bam_out != 0, "Cannot write to " + cstrtos(outF), "!");

  bam_header_destroy(out_header);
}

BamWriter::~BamWriter() {
  samclose(bam_out);
}

#endif
