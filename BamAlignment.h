#ifndef BAMALIGNMENT_H_
#define BAMALIGNMENT_H_

#include<cmath>
#include<cassert>
#include<string>
#include<algorithm>

#include "stdint.h"
#include "sam/bam.h"
#include "sam/sam.h"

#include "utils.h"
#include "my_assert.h"

class BamAlignment {
 public:
  BamAlignment();
  BamAlignment(const BamAlignment&);
  ~BamAlignment();

  bool read(samfile_t *);
  bool write(samfile_t *);
  
  void writeToBED(FILE*, const bam_header_t*);
  void writeToTagAlign(FILE*, const bam_header_t*);

  bool isPaired() const { return is_paired; }

  // isAligned should be called first in order to call all get and set functions
  bool isAligned() const {
    if (b->core.flag & 0x0004) return false;
    if (is_paired && (b2->core.flag & 0x0004)) return false;
    return true;
  }

  std::string getName() const { return std::string((char*)bam1_qname(b)); }
  
  CHR_ID_TYPE getCid() const { return b->core.tid; }

  CHR_LEN_TYPE getPos() const { return b->core.pos; }

  // length of the query sequence
  int getSeqLength() const { return b->core.l_qseq; } 

  int getISize() const { return b->core.isize; }

  // if length = 2k, midpoint is k - 1
  CHR_LEN_TYPE getMidPos(int fragment_length) const {
    if (is_paired && getISize() > 0) return getPos() + (getISize() - 1) / 2;
    return (getDir() == '+' ? getPos() + (fragment_length - 1) / 2 : getPos() + getSeqLength() - fragment_length / 2 - 1);
  }

  char getDir() const { 
    if (is_paired && !(b->core.flag & 0x0004) && !(b2->core.flag & 0x0004)) 
      return (b->core.flag & 0x0040) ? '+' : '-';
    return (b->core.flag & 0x0010) == 0 ? '+' : '-'; 
  }

  float getFrac() const { 
    uint8_t *p_tag = bam_aux_get(b, "ZW");
    return p_tag != NULL ? bam_aux2f(p_tag) : 1.0;
  }
  
  void setFrac(float frac) {
    uint8_t *p_tag;

    b->core.qual = getMAPQ(frac);

    p_tag = bam_aux_get(b, "ZW");
    if (p_tag != NULL) {
      memcpy(p_tag + 1, (uint8_t*)&(frac), bam_aux_type2size('f'));
    } else {
      bam_aux_append(b, "ZW", 'f', bam_aux_type2size('f'), (uint8_t*)&frac);
    }

    if (is_paired && !(b2->core.flag & 0x0004)) {
      b2->core.qual = b->core.qual;

      p_tag = bam_aux_get(b2, "ZW");
      if (p_tag != NULL) {
	memcpy(p_tag + 1, (uint8_t*)&(frac), bam_aux_type2size('f'));
      } else {
	bam_aux_append(b2, "ZW", 'f', bam_aux_type2size('f'), (uint8_t*)&frac);
      }
    }

  }

  void removeZWTag() {
    uint8_t *p_tag = bam_aux_get(b, "ZW");
    if (p_tag != NULL) bam_aux_del(b, p_tag);
  }

  bool isFiltered() {
    uint8_t *p = bam_aux_get(b, "XM");
    if (p != NULL && bam_aux2i(p) > 0) return true;

    if (is_paired) {
      p = bam_aux_get(b2, "XM");
      if (p != NULL && bam_aux2i(p) > 0) return true; 
    }

    return false;
  }

  std::string getFQScore() { return getQScore(b); }
  
  std::string getRQScore() { assert(is_paired); return getQScore(b2); }

  std::string getFSeq() { return getSeq(b); }
  
  std::string getRSeq() { assert(is_paired); return getSeq(b2); }

 private:

  bool is_paired;
  bam1_t *b, *b2;

  uint8_t getMAPQ(float val) {
    float err = 1.0 - val;
    if (err <= 1e-10) return 100;
    return (uint8_t)(-10 * log10(err) + .5); // round it
  }

  std::string getQScore(const bam1_t* b) {
    uint8_t *p = bam1_qual(b);
    std::string qscore = "";

    if (!(b->core.flag & 0x0010)) {
      for (int i = 0; i < b->core.l_qseq; i++) {
	qscore.append(1, (char)(*p + 33));
	++p;
      }
    }
    else {
      p = p + b->core.l_qseq - 1;
      for (int i = 0; i < b->core.l_qseq; i++) {
	qscore.append(1, (char)(*p + 33));
	--p;
      }
    }

    return qscore;
  }

  std::string getSeq(const bam1_t* b) {
    uint8_t *p = bam1_seq(b);
    std::string readseq = "";
    char base = 0;

    // b must represent an alignment of an aligned read, otherwise, 0x0010 is not reliable.
    if (b->core.flag & 0x0010) {
      for (int i = b->core.l_qseq - 1; i >= 0; i--) {
	switch(bam1_seqi(p, i)) {
	  //case 0 : base = '='; break;
	case 1 : base = 'T'; break;
	case 2 : base = 'G'; break;
	case 4 : base = 'C'; break;
	case 8 : base = 'A'; break;
	case 15 : base = 'N'; break;
	default : assert(false);
	}
	readseq.append(1, base);
      }
    }
    else {
      for (int i = 0; i < b->core.l_qseq; i++) {
	switch(bam1_seqi(p, i)) {
	  //case 0 : base = '='; break;
	case 1 : base = 'A'; break;
	case 2 : base = 'C'; break;
	case 4 : base = 'G'; break;
	case 8 : base = 'T'; break;
	case 15 : base = 'N'; break;
	default : assert(false);
	}
	readseq.append(1, base);
      }
    }

    return readseq;
  }

};

BamAlignment::BamAlignment() : is_paired(false), b(bam_init1()), b2(bam_init1()) {
}

BamAlignment::BamAlignment(const BamAlignment& o) : b(NULL), b2(NULL) {
  is_paired = o.is_paired;
  b = bam_dup1(o.b);
  b2 = bam_dup1(o.b2);
}
  
BamAlignment::~BamAlignment() {
  bam_destroy1(b);
  bam_destroy1(b2);
}

bool BamAlignment::read(samfile_t *in) {
  if (samread(in, b) < 0) return false;
  is_paired = (b->core.flag & 0x0001) > 0;
  if (is_paired) { 
    general_assert(samread(in, b2) >= 0 && (b2->core.flag & 0x0001), "Fail to read the other mate for a paired-end alignment!");
    general_assert(((b->core.flag & 0x00C0) == 0x0040 && (b2->core.flag & 0x00C0) == 0x0080) || 
		   ((b->core.flag & 0x00C0) == 0x0080 && (b2->core.flag & 0x00C0) == 0x0040), 
		   "Cannot detect both mates of a paired-end alignment!");
    
    if ((b->core.flag & 0x0004) && !(b2->core.flag & 0x0004)) { bam1_t* tmp = b; b = b2; b2 = tmp; }  // if one of the mate can be aligned but not the other, switch them

    if (!(b->core.flag & 0x0004) && !(b2->core.flag & 0x0004) && b->core.pos > b2->core.pos) {
      bam1_t* tmp = b; b = b2; b2 = tmp; // switch to make the leftmost segment as b
      assert(b->core.tid == b2->core.tid && b->core.isize > 0);
    }
  }
  
  return true;
}

bool BamAlignment::write(samfile_t *out) {
  general_assert(samwrite(out, b) >= 0, "Fail to write alignments to BAM file!");
  if (is_paired) general_assert(samwrite(out, b2) >= 0, "Fail to write alignments to BAM file!");
  return true;
}

void BamAlignment::writeToBED(FILE *fo, const bam_header_t *header) {
  if (b->core.flag & 0x0004) return;
  fprintf(fo, "%s\t%d\t%d\t%s\t%.2f\t%c\n", header->target_name[b->core.tid], b->core.pos, b->core.pos + b->core.l_qseq, (char*)bam1_qname(b), getFrac() * 1000.0, ((b->core.flag & 0x0010) == 0 ? '+' : '-'));
  if (!is_paired || (b2->core.flag & 0x0004)) return;
  fprintf(fo, "%s\t%d\t%d\t%s\t%.2f\t%c\n", header->target_name[b2->core.tid], b2->core.pos, b2->core.pos + b2->core.l_qseq, (char*)bam1_qname(b2), getFrac() * 1000.0, ((b2->core.flag & 0x0010) == 0 ? '+' : '-'));
}

void BamAlignment::writeToTagAlign(FILE *fo, const bam_header_t *header) {
  if (b->core.flag & 0x0004) return;
  general_assert(!is_paired, "Paired-end alignments are detected for tagAlign format!");

  fprintf(fo, "%s\t%d\t%d\t%s\t%.2f\t%c\n", header->target_name[b->core.tid], b->core.pos, b->core.pos + b->core.l_qseq, getSeq(b).c_str(), getFrac() * 1000.0, ((b->core.flag & 0x0010) == 0 ? '+' : '-'));
}

#endif
