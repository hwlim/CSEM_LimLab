#ifndef BCTREE_H_
#define BCTREE_H_

#include<cstdio>
#include<cstring>
#include<cassert>
#include<vector>

#include "utils.h"

class bcTree {
 public:
  bcTree() { n = 0; axis.clear(); sum.clear(); sum_uniq.clear(); }
  ~bcTree() {}
  void addPos(CHR_LEN_TYPE);
  void clear();
  double count(CHR_LEN_TYPE);
  void update(CHR_LEN_TYPE, double);
  void showTree();

  CHR_LEN_TYPE getN() const { return n; }

  // back up the bcTree counts with only unique reads
  void backup() { sum_uniq = sum; }
  // new round of update begins, copy unique reads only counts to sum
  void init() { sum = sum_uniq; }

 private:
  CHR_LEN_TYPE n; // n elements;
  std::vector<CHR_LEN_TYPE> axis; // the axis of elements
  std::vector<double> sum; // the sum of corresponding substrees;
  std::vector<double> sum_uniq; // a copy of bcTree counts with unique reads only
};

// must be added in ascending order
void bcTree::addPos(CHR_LEN_TYPE pos) {
  axis.push_back(pos);
  sum.push_back(0.0);
  ++n;
}

void bcTree::clear() {
  if (n > 0) sum.assign(n, 0.0);
}

double bcTree::count(CHR_LEN_TYPE pos) {
  CHR_LEN_TYPE l, r, mid;
  double res = 0.0;

  l = 0; r = n - 1;
  while (l <= r) {
    mid = (l + r) >> 1;
    if (axis[mid] == pos) {
      res += sum[mid];
      break;
    }
    else if (axis[mid] < pos) {
      res += sum[mid];
      l = mid + 1;
    }
    else {
      r = mid - 1;
    }
  }

  return res;
}

// val is the difference to update
void bcTree::update(CHR_LEN_TYPE pos, double val) {
  CHR_LEN_TYPE l, r, mid;
  
  l = 0; r = n - 1;
  while (l <= r) {
    mid = (l + r) >> 1;
    
    if (axis[mid] == pos) {
      sum[mid] += val;
      return;
    }
    else if (axis[mid] < pos) {
      l = mid + 1;
    }
    else {
      sum[mid] += val;
      r = mid - 1;
    }
  }
  assert(false); // if cannot find pos, assert(false) to report error
}

void bcTree::showTree() {
  if (n == 0) return;
  for (CHR_LEN_TYPE i = 0; i < n; i++) 
    printf("%.6lf\t", sum[i]);
  printf("\n");
}

#endif
