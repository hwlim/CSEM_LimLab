// sampling.h for CSEM

#ifndef SAMPLING
#define SAMPLING

#include<cstdio>
#include<cassert>
#include<vector>

#include "boost/random.hpp"

typedef boost::mt19937 engine_type;
typedef boost::uniform_01<engine_type> uniform01;

// arr should be cumulative!
// interval : [,)
// random number should be in [0, arr[len - 1])
// If by chance arr[len - 1] == 0.0, one possibility is to sample uniformly from 0...len-1
int sample(uniform01& rg, std::vector<double>& arr, int len) {
  int l, r, mid;
  double prb = rg() * arr[len - 1];

  l = 0; r = len - 1;
  while (l <= r) {
    mid = (l + r) / 2;
    if (arr[mid] <= prb) l = mid + 1;
    else r = mid - 1;
  }

  if (l >= len) { printf("%d %lf %lf\n", len, arr[len - 1], prb); }
  assert(l < len);

  return l;
}

#endif
