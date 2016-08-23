// A JM state is a proton or neutron state with good J and M, expressed in an M-scheme basis
// It should be formed by combining NuProj and NuBasis in a straightforward way
// These then form the basis in which the eigenvectors in NuVec are expressed.
#ifndef JMState_h
#include "NuBasis.hh"
#include "NuProj.hh"
#include <vector>
#include <string>
#include <unordered_map>

using namespace std;


struct KeyHash
{
  size_t operator() (const vector<mvec_type>& key) const
  {
    size_t hash_out = 0;
    for (size_t i=0;i<key.size();++i)  hash_out |= (key[i]<<i);
    return hash_out;
  }
};


class JMState
{
 public:
  int J2, T2, M2, pindx;
  unordered_map<vector<mvec_type>,float,KeyHash> m_coefs;
  vector<MschemeOrbit> m_orbits;

  JMState();
  JMState(const JMState&)=default;
  JMState(JMState&&)=default;
  JMState(const NuBasis& nubasis, const NuProj& nuproj, int istate);

  void Print() const;
  string PrintMstate(vector<mvec_type> vec_in) const;

  void Normalize();
  void EliminateZeros();
  JMState Jplus();
  JMState Jminus();
  JMState TimeReverse();
  void RotateToM(int M);
  void SetJ(int j){J2 = j;};
  void SetM(int m){M2 = m;};

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

vector<mvec_type> operator+( const vector<mvec_type>& lhs, const vector<mvec_type>& rhs);

#define JMState_h
#endif
