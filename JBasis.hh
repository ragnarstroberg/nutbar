#ifndef JBasis_h

#include <vector>

#include "NuBasis.hh"
#include "NuProj.hh"
#include "JMState.hh"

using namespace std;

struct ab_pair{  size_t ia;  size_t ib; ab_pair(int a,int b): ia(a), ib(b){}; };

class JBasis
{
 public:
  int J2,M2;
  int N_p, N_n; // these are not used.
//  vector<JMState> basis_states;
//  vector<array<int,4>> basis_states;
//  vector<array<int,3>> basis_states;
  map<int,vector<ab_pair>> basis_states;
  vector<JMState> jmstates_a, jmstates_b;
  vector<MschemeOrbit> m_orbits;

  JBasis();
//  JBasis(int j2, int m2);
//  JBasis( string sps_file,  vector<string> proton_files, vector<string> neutron_files, int j2, int m2 );
  JBasis( string sps_file,  vector<string> proton_files, vector<string> neutron_files, vector<int>& J2list );

  void SetupBasis( string sps_file,  vector<string> A_files, vector<string> B_files, vector<int>& J2list);
  void AddBasisStates_J( size_t nA, size_t nB, vector<int>& offsets_a, vector<int>& offsets_b,  int J2  );
  void SetUpJMState_ab(  NuBasis& nubasis, NuProj& nuproj,  vector<string> filenames, vector<JMState>& jmstates,  vector<int>& offsets);
//  JMState GetBasisState(size_t index) const;
  JMState GetBasisState(size_t index, int J2, int M2) const;

  int GetNaiveMschemeDimension(int J2) const;

};




#define JBasis_h
#endif
