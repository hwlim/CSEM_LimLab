#ifndef ARRAYSCAN_H_
#define ARRAySCAN_H_

#include<cassert>
#include<vector>

#include "utils.h"

class ArrayScan {
 public:  
  ArrayScan(int, const std::vector<CHR_LEN_TYPE>&, const std::vector<double>&);
  
  // positions should be larger than the previous ones
  double getSumBy(CHR_LEN_TYPE);
  double getValueAt(CHR_LEN_TYPE);

 private:
  CHR_LEN_TYPE s; // number of elements 
  int type; // 0, end-value; 1, continue-value
  const std::vector<CHR_LEN_TYPE>& lengths;
  const std::vector<double>& values;

  CHR_LEN_TYPE nextp, curpos; // nextp, next position in lengths/values; curpos, current postion in the chromosome
  double sum;
};

ArrayScan::ArrayScan(int type, const std::vector<CHR_LEN_TYPE>& lengths, const std::vector<double>& values) : type(type), lengths(lengths), values(values) {
  s = (CHR_LEN_TYPE)lengths.size();
  nextp = 0; 
  curpos = -1;
  sum = 0.0;
}

double ArrayScan::getSumBy(CHR_LEN_TYPE pos) {
  if (curpos < 0 && pos <= curpos) return 0.0;
  assert(pos > curpos);
  while (nextp < s && curpos + lengths[nextp] < pos) {
    curpos += lengths[nextp];
    sum += (type == 0 ? values[nextp] : values[nextp] * lengths[nextp]);
    ++nextp;
  }
  if (nextp == s) return sum;

  double result = sum;
  if (type == 0) result += (curpos + lengths[nextp] == pos ? values[nextp] : 0.0);
  else result += values[nextp] * (pos - curpos);

  return result;
}

double ArrayScan::getValueAt(CHR_LEN_TYPE pos) {
  assert(pos > curpos);
  if (nextp == s) return 0.0;

  do {
    curpos += lengths[nextp++];
  } while (curpos < pos && nextp < s);

  if (curpos < pos) return 0.0;

  --nextp;
  double result = (type == 1 || curpos == pos ? values[nextp] : 0.0); 
  curpos -= lengths[nextp];

  return result;
}

#endif
