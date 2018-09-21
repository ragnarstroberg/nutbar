#ifndef NuBasis_h

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdint>
#include <bitset>

#define KEY_BITS 64  // number of bits used in the definition of key_type below

typedef int16_t part_type;
typedef int64_t mvec_type;
typedef std::bitset<KEY_BITS> key_type;


//struct MschemeOrbit;
struct MschemeOrbit
{
  int n,l2,j2,mj2,tz2;

  MschemeOrbit(int n, int l2, int j2, int mj2, int tz2): n(n),l2(l2),j2(j2),mj2(mj2),tz2(tz2){};
};



struct NuBasis
{

  int32_t no_spart,sum_nJT;
  std::vector<int32_t> pindx,nJT,j,t,ibf,max_cM2,max_cTz2;
  std::vector<std::vector<part_type>> partition;
  std::vector<std::vector<key_type>> vec;
  static const int gwords;
  std::vector<MschemeOrbit> m_orbits;

  NuBasis();
  void ReadFile(std::string fname);
  void ReadSPS(std::string fname);
  void PrintBasis();
  void PrintSPS();
  void Clear();


};





#define NuBasis_h
#endif

