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

//typedef int16_t part_type;
//typedef int64_t mvec_type;

class NuVec
{
 public:

  int32_t no_state,no_level; // states are the basis vectors, levels are the eigenvectors.
  int J2;
  vector<float> alpha; //eigenvalues
  vector<vector<float>> coefT; // eigenvectors.
  
  NuVec();
  NuVec(int J2);
  void ReadFile(string fname);
  void PrintVectors();
  void PrintDetailedVectors();
  bool CheckOrthoNormal();

};

#define NuVec_h
#endif
