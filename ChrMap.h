#ifndef CHRMAP_H_
#define CHRMAP_H_

#include<cassert>
#include<map>
#include<string>
#include<vector>

#include "sam/bam.h"

#include "utils.h"

class ChrMap {
 public:
  
  ChrMap(const bam_header_t* header) { 
    s = header->n_targets;
    chrLens.assign(header->target_len, header->target_len + s);
    for (CHR_ID_TYPE cid = 0; cid < s; cid++) {
      std::string cname(header->target_name[cid]);
      assert(cname2cid.find(cname) == cname2cid.end());
      cname2cid[cname] = cid;
    }
  }
 
  CHR_ID_TYPE size() { return s; }

  CHR_LEN_TYPE getLen(CHR_ID_TYPE cid) { 
    assert(cid >= 0 && cid < s);
    return chrLens[cid]; 
  }

  CHR_ID_TYPE getCid(std::string cname) {
    iter = cname2cid.find(cname);
    assert(iter != cname2cid.end());
    return iter->second;
  }

 private:
  CHR_ID_TYPE s;
  std::vector<CHR_LEN_TYPE> chrLens;
  std::map<std::string, CHR_ID_TYPE> cname2cid;
  std::map<std::string, CHR_ID_TYPE>::iterator iter;

};

#endif
