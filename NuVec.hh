#ifndef NuVec_h

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdint>


struct NuVec
{

  int32_t no_state,no_level; // states are the basis vectors, levels are the eigenvectors.
  int J2;
  std::vector<float> alpha; //eigenvalues
  std::vector<std::vector<float>> coefT; // eigenvectors.

  static const int delimiter_length;
  
  NuVec();
  NuVec(int J2);
//  void ReadFile(std::string fname);
  void ReadFile(std::string fname, int32_t max_levels=1e9);
  void PrintVectors();
  void PrintDetailedVectors();
  bool CheckOrthoNormal();

};

#define NuVec_h
#endif
