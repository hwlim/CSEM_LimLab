#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include "utils.h"

struct Alignment {
  CHR_ID_TYPE cid;
  CHR_LEN_TYPE pos;
  double frac;

  char dir; // + or -
  bool isMulti; 

  Alignment() { cid = -1; pos = -1; frac = 0.0; dir = 0; isMulti = false; }

  Alignment(CHR_ID_TYPE cid, CHR_LEN_TYPE pos, char dir, double frac = 0.0, bool isMulti = false) {
    this->cid = cid;
    this->pos = pos;
    this->dir = dir;
    this->frac = frac;
    this->isMulti = isMulti;
  }

  bool operator< (const Alignment& o) const {
    return cid < o.cid || (cid == o.cid && pos < o.pos);
  }

  bool operator== (const Alignment& o) const {
    return cid == o.cid && pos == o.pos;
  }
};

#endif
