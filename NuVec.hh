#ifndef NuVec_h

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

class NuVec
{
 public:

  int32_t no_state,no_level;
  vector<float> alpha; //eigenvalues
  vector<vector<float>> coefT;
  
  NuVec();
  void ReadFile(string fname);
  void PrintVectors();
  void PrintDetailedVectors();
  bool CheckOrthoNormal();

};

#define NuVec_h
#endif
