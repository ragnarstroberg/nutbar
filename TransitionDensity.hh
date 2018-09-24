
#ifndef TransitionDensity_h

#include "NuBasis.hh"
#include "NuProj.hh"
#include "NuVec.hh"
#include "JMState.hh"
#include "JBasis.hh"
#include "Profiler.hh"
#include "Settings.hh"

#include <vector>
#include <armadillo>
#include <unordered_map>


class TransitionDensity
{
 public:
  
  std::vector<std::vector<float>> blank_vector_i;
  std::vector<std::vector<float>> blank_vector_f;
  std::vector<MschemeOrbit> m_orbits;
  std::vector<int> Jlist_i;
  std::vector<int> Jlist_f;
  int MJtot_i;
  int MJtot_f;
  std::vector<NuVec> nuvec_list_i;
  std::vector<NuVec> nuvec_list_f;
  std::unordered_map< key_type, std::vector<std::vector<float>> > amplitudes_i;  // this is the moneymaker.
  std::unordered_map< key_type, std::vector<std::vector<float>> > amplitudes_f;


  TransitionDensity();
  TransitionDensity(std::vector<int> jlist);
  TransitionDensity(std::vector<int> jlist_i, std::vector<int> jlist_f);

  void CalculateMschemeAmplitudes(Settings& settings);
  void CalculateMschemeAmplitudes_fi(std::vector<NuVec>& , JBasis& ,  std::vector<std::vector<float>>& blank_vector, std::unordered_map< key_type, std::vector<std::vector<float>> >& amplitudes);


  double OBTD(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a, int m_index_b, int Lambda2 );
  double TBTD(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a, int m_index_b, int m_index_c, int m_index_d, int J2ab, int J2cd, int Lambda2 );


  double TD_ax(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a );
  double TD_axaxa(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a, int m_index_b, int m_index_c, int J2ab, int Lambda2);

  arma::mat CalcOBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, Settings& settings);
  arma::mat CalcTBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, Settings& settings);

  arma::vec CalcTransitionDensity_ax( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, Settings& settings);
  arma::mat CalcTransitionDensity_axaxa( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, Settings& settings);


  void SetupKets();
  void Jplus(std::vector<key_type>& mvecs_in, std::vector<double>& amp_in, int J2, int M2);

};



// Some useful operator overrides
std::vector<float> operator*(const float lhs, const std::vector<float>& rhs);
std::vector<float> operator*(const std::vector<float>& lhs, const std::vector<float>& rhs);
std::vector<float>& operator+=(std::vector<float>& lhs, const std::vector<float>& rhs);




#define TransitionDensity_h
#endif
