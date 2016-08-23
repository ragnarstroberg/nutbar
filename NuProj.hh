#ifndef NuProj_h

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdint>

#define DELIMITER_LENGTH 4

using namespace std;

typedef int16_t part_type;
typedef int64_t mvec_type;

class NuProj
{
 public:

  int32_t nptt,ngood;
  vector<int32_t> no_spart,pindx,j,t,dim,max_cM2,max_cTz2;
  vector<float> x;
  vector<vector<float>> coef_st;
  
  NuProj();
  void ReadFile(string fname);
  void PrintProj();
  void Clear();

};

#define NuProj_h
#endif
