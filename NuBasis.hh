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

using namespace std;

typedef int16_t part_type;
typedef int64_t mvec_type;


struct MschemeOrbit;

class NuBasis
{
 public:

  int32_t no_spart,sum_nJT;
  vector<int32_t> pindx,nJT,j,t,ibf,max_cM2,max_cTz2;
  vector<vector<part_type>> partition;
  vector<vector<vector<mvec_type>>> vec;
  static const int gwords;
  vector<MschemeOrbit> m_orbits;

  NuBasis();
  void ReadFile(string fname);
  void ReadSPS(string fname);
  void PrintBasis();
  void PrintSPS();
  void Clear();

  //TODO: These should be applied to a J state
  void LoweringOperator( vector<mvec_type>& mvec_in, vector<vector<mvec_type>>& mvec_out, vector<float>& coefs);
  void RaisingOperator( vector<mvec_type>& mvec_in, vector<vector<mvec_type>>& mvec_out, vector<float>& coefs);
  vector<mvec_type> TimeReverse( vector<mvec_type>& mvec_in);

};


struct MschemeOrbit
{
  int n,l2,j2,mj2,tz2;

  MschemeOrbit(int n, int l2, int j2, int mj2, int tz2): n(n),l2(l2),j2(j2),mj2(mj2),tz2(tz2){};
};



#define NuBasis_h
#endif

