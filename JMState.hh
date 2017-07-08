// A JM state is a proton or neutron state with good J and M, expressed in an M-scheme basis
// It should be formed by combining NuProj and NuBasis in a straightforward way
// These then form the basis in which the eigenvectors in NuVec are expressed.
#ifndef JMState_h
#include "NuBasis.hh"
#include "NuProj.hh"
#include "Profiler.hh"
#include <vector>
#include <string>
#include <unordered_map>

using namespace std;



class JMState
{
 public:
  int J2, T2, M2;//, pindx;
  vector<int> pindx;
  vector<part_type> partition;
  unordered_map<key_type,float> m_coefs;
  vector<MschemeOrbit> m_orbits;
  Profiler profiler;

  JMState();
  JMState(const JMState&)=default;
  JMState(JMState&&)=default;
  JMState(const NuBasis& nubasis, const NuProj& nuproj, int istate);

  void Print() const;
  string PrintMstate(key_type vec_in) const;

  void Normalize();
  void EliminateZeros();
  JMState Jplus();
  JMState Jminus();
  JMState TimeReverse();
  void RotateToM(int M);
  void SetJ(int j){J2 = j;};
  void SetM(int m){M2 = m;};

  float Norm() const;

  JMState OuterProduct( const JMState&) const;

  JMState& operator=(const JMState&);
  JMState& operator=( JMState&&);
  JMState& operator+=(const JMState&);
  JMState operator+(const JMState&);
  JMState& operator*=(const double);
  JMState operator*(const double);

};


JMState TensorProduct( JMState jm1, JMState jm2, int J2, int M2);
JMState operator*(const double, const JMState&);

key_type operator+( const key_type& lhs, const key_type& rhs);

//double CG(int j2a, int m2a, int j2b, int m2b, int J2, int M2);

#define JMState_h
#endif
