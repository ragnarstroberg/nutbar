#ifndef JBasis_h

#include <vector>
#include <map>
#include <string>

#include "NuBasis.hh"
#include "NuProj.hh"
#include "JMState.hh"


struct ab_pair{  size_t ia;  size_t ib; ab_pair(int a,int b): ia(a), ib(b){}; };

struct JBasis
{
  int J2,M2;
  int N_p, N_n; // these are not used.
  std::map<int,std::vector<ab_pair>> basis_states;
  std::vector<JMState> jmstates_a, jmstates_b;
  std::vector<MschemeOrbit> m_orbits;

  JBasis();
  JBasis( std::string sps_file,  std::vector<std::string> proton_files, std::vector<std::string> neutron_files, std::vector<int>& J2list );

  void SetupBasis( std::string sps_file,  std::vector<std::string> A_files, std::vector<std::string> B_files, std::vector<int>& J2list);
  void AddBasisStates_J( size_t nA, size_t nB, std::vector<int>& offsets_a, std::vector<int>& offsets_b,  int J2  );
  void SetUpJMState_ab(  NuBasis& nubasis, NuProj& nuproj,  std::vector<std::string> filenames, std::vector<JMState>& jmstates,  std::vector<int>& offsets);
  JMState GetBasisState(size_t index, int J2, int M2) const;

  int GetNaiveMschemeDimension(int J2) const;

};




#define JBasis_h
#endif
