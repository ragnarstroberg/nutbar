#ifndef NuProj_h

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdint>


typedef int16_t part_type;
typedef int64_t mvec_type;

struct NuProj
{
  int32_t nptt,ngood;
  std::vector<int32_t> no_spart,pindx,j,t,dim,max_cM2,max_cTz2;
  std::vector<float> x;
  std::vector<std::vector<float>> coef_st;
  static const int delimiter_length;
  
  NuProj();
  void ReadFile(std::string fname);
  void PrintProj();
  void Clear();

};

#define NuProj_h
#endif
