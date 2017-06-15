#ifndef JBasis_h

#include <vector>

#include "NuBasis.hh"
#include "NuProj.hh"
#include "JMState.hh"

using namespace std;

class JBasis
{
 public:
  int J2,M2;
  int N_p, N_n;
//  vector<JMState> basis_states;
  vector<array<int,4>> basis_states;
  vector<JMState> jmstates_a, jmstates_b;
  vector<MschemeOrbit> m_orbits;

  JBasis();
  JBasis(int j2, int m2);
  JBasis( string sps_file,  vector<string> proton_files, vector<string> neutron_files, int j2, int m2 );

  void SetupBasis( string sps_file,  vector<string> A_files, vector<string> B_files);
  void SetUpJMState_ab(  string sps_file,  vector<string> filenames, vector<JMState>& jmstates,  vector<int>& offsets);
  JMState GetBasisState(size_t index) const;

};




#define JBasis_h
#endif
