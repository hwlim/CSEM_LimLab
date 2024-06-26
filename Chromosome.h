#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

#include<cmath>
#include<cassert>
#include<sstream>
#include<vector>
#include<algorithm>

#include "utils.h"

#include "Alignment.h"
#include "ArrayScan.h"

class Chromosome {
 public:
  Chromosome(int, CHR_LEN_TYPE, std::vector<Alignment>&);
  Chromosome(const Chromosome&); // this is only used for sorting!

  HIT_INT_TYPE getSize() const { return size; }

  void addPos(HIT_INT_TYPE, bool);
  void init();
  void processPriorInfo(const std::string&, int, const std::vector<double>&);
  void update(bool);

  double getMaxDelta() const { return max_delta; }

  bool operator()(HIT_INT_TYPE a, HIT_INT_TYPE b) const {
    return alignments[a].pos < alignments[b].pos;
  }

  
 private:
  int halfws;
  CHR_LEN_TYPE clen; // chromosome length
  std::vector<Alignment>& alignments; 


  HIT_INT_TYPE size; // size, total number of alignments
  std::vector<HIT_INT_TYPE> alignPos; // positions in "alignments" vector for multi-reads

  CHR_LEN_TYPE s, offset; // s, total number of unique multi-read alignment positions; offset, where genomic coordinate >= 0
  std::vector<CHR_LEN_TYPE> coords; // discretized coordinates for multi-read alignments

  std::vector<CHR_LEN_TYPE> lengths; // this vector is used by ArrayScan
  std::vector<double> values; // multiread fractions

  std::vector<HIT_INT_TYPE> uniqPos; // positions in "alignments" vector for unique reads

  std::vector<double> baseWindowSums; // constant part of sum in a window, including unique reads and prior counts 
  std::vector<double> basePointValues; // point values at multi-read positions, including unique reads and prior info

  double max_delta;
};

Chromosome::Chromosome(int halfws, CHR_LEN_TYPE clen, std::vector<Alignment>& alignments) : halfws(halfws), clen(clen), alignments(alignments) { 
  size = 0; 
  alignPos.clear();
  uniqPos.clear();
  max_delta = 0.0;
}

Chromosome::Chromosome(const Chromosome& o) : alignments(o.alignments) { }

inline void Chromosome::addPos(HIT_INT_TYPE pos, bool isMulti) {
  if (isMulti) { alignPos.push_back(pos); ++size; }
  else { uniqPos.push_back(pos); }
}

void Chromosome::init() {
  CHR_LEN_TYPE pos; 
  CHR_LEN_TYPE prevpos, curpos, curidx; // these two are for genomic coordinates >= 0 && < clen only
  HIT_INT_TYPE usize; // usize, size of uniqPos

  std::vector<CHR_LEN_TYPE> lens;
  std::vector<double> vals;

  // for multi-read alignments
  assert(size == (HIT_INT_TYPE)alignPos.size());
  std::sort(alignPos.begin(), alignPos.end(), *this);

  coords.clear();
  s = 0; offset = -1; 
  prevpos = -1;
  lengths.clear();
  for (HIT_INT_TYPE i = 0; i < size; i++) {
    pos = alignments[alignPos[i]].pos;
    if (i + 1 == size || pos != alignments[alignPos[i + 1]].pos) {
      if (pos >= 0 && pos < clen) {
	if (offset < 0) offset = coords.size();
	assert(prevpos < pos);
	lengths.push_back(pos - prevpos);
	prevpos = pos;
      }
      coords.push_back(pos);
    }
  }
  s = coords.size();
  values.assign(lengths.size(), 0.0);

  assert(offset < 0 || ((offset < 1 || coords[offset - 1] < 0) && (coords[offset] >= 0 && coords[offset] < clen)));

  // for unique-read alignments
  std::sort(uniqPos.begin(), uniqPos.end(), *this);

  usize = uniqPos.size(); 
  curpos = -1; curidx = -1;
  lens.clear(); vals.clear();
  for (HIT_INT_TYPE i = 0; i < usize; i++) {
    pos = alignments[uniqPos[i]].pos;
    if (pos < 0) continue;
    if (pos >= clen) break;

    if (curpos < pos) {
      lens.push_back(pos - curpos);
      vals.push_back(0.0);
      curpos = pos; 
      ++curidx;
    }
    vals[curidx] += alignments[uniqPos[i]].frac;
  }


  ArrayScan arrScanL(0, lens, vals);
  ArrayScan arrScanU(0, lens, vals);
  ArrayScan arrScanM(0, lens, vals);

  baseWindowSums.assign(s, 0.0);
  basePointValues.assign(s, 0.0);

  for (CHR_LEN_TYPE i = 0; i < s; i++) {
    baseWindowSums[i] = arrScanU.getSumBy(coords[i] + halfws) - arrScanL.getSumBy(coords[i] - halfws - 1);
    if (coords[i] >= 0 && coords[i] < clen) basePointValues[i] = arrScanM.getValueAt(coords[i]);
  }
}

void Chromosome::processPriorInfo(const std::string& line, int ngroups, const std::vector<double>& groupvalues) {
  int gid;
  CHR_LEN_TYPE len;
  std::istringstream strin(line);
  std::vector<CHR_LEN_TYPE> lens;
  std::vector<double> vals;

  CHR_LEN_TYPE sum = 0;

  lens.clear(); vals.clear();
  while (strin>> len>> gid) {
    lens.push_back(len);
    assert(gid >= 0 && gid < ngroups);
    vals.push_back(groupvalues[gid]);
    sum += len;
  }

  assert(sum == clen);

  ArrayScan arrScanL(1, lens, vals);
  ArrayScan arrScanU(1, lens, vals);
  ArrayScan arrScanM(1, lens, vals);

  for (CHR_LEN_TYPE i = 0; i < s; i++) {
    baseWindowSums[i] += arrScanU.getSumBy(coords[i] + halfws) - arrScanL.getSumBy(coords[i] - halfws - 1);
    if (coords[i] >= 0 && coords[i] < clen) basePointValues[i] += arrScanM.getValueAt(coords[i]);
  }
}

void Chromosome::update(bool updateFrac = true) {
  CHR_LEN_TYPE pos, curidx;
  double value;

  // update values
  max_delta = 0.0;

  curidx = offset;
  value = 0.0;
  for (HIT_INT_TYPE i = 0; i < size; i++) {
    pos = alignments[alignPos[i]].pos;
    if (pos < 0) continue;
    if (pos >= clen) break;

    if (pos > coords[curidx]) {
      if (value + basePointValues[curidx] < 0.0) value = -basePointValues[curidx];
      max_delta = std::max(max_delta, fabs(values[curidx - offset] - value));
      values[curidx - offset] = value;      

      ++curidx;
      value = 0.0;
    }

    value += alignments[alignPos[i]].frac;
  }
  if (curidx >= 0) {
    assert(curidx < s && coords[curidx] >= 0 && coords[curidx] < clen);
    if (value + basePointValues[curidx] < 0.0) value = -basePointValues[curidx];
    max_delta = std::max(max_delta, fabs(values[curidx - offset] - value));
    values[curidx - offset] = value;
  }
 
  if (!updateFrac) return;

  // update alignments.frac for the multi-read alignments in this chromosome
  ArrayScan arrScanL(0, lengths, values);
  ArrayScan arrScanU(0, lengths, values);

  curidx = -1; value = 0.0;
  for (HIT_INT_TYPE i = 0; i < size; i++) {
    pos = alignments[alignPos[i]].pos;
    if (curidx < 0 || pos > coords[curidx]) {
      ++curidx;
      assert(curidx < s && pos == coords[curidx]);
      value = baseWindowSums[curidx] + (arrScanU.getSumBy(coords[curidx] + halfws) - arrScanL.getSumBy(coords[curidx] - halfws - 1));     
      assert(value >= 0.0);
    }
    alignments[alignPos[i]].frac = value;
  }
}

#endif
