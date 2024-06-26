#ifndef UTILS_H_
#define UTILS_H_

#include<stdint.h>

// change to uint64_t if too many data are there
typedef uint32_t HIT_INT_TYPE;
typedef uint32_t READ_INT_TYPE;
typedef int32_t CHR_ID_TYPE;
typedef int32_t CHR_LEN_TYPE; // must be signed type , the coordinates can be negative for some cases

const int STRLEN = 10005 ;

#endif
