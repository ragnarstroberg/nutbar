
#ifndef TransitionDensity_h

#include "NuBasis.hh"
#include "NuProj.hh"
#include "NuVec.hh"
#include "JMState.hh"
#include "JBasis.hh"
#include "Profiler.hh"

#include <vector>
#include <armadillo>
#include <unordered_map>


class TransitionDensity
{
 public:
  
  std::string basename_i; // used in ReadFiles, GetAZFromFileName, CalculateMschemeAmplitudes, CalcOBTD, CalcTBTD, CalcTransitionDensity_ax, ReadOBTD, ReadTBTD.  All but the 1st 2 check initial==final.
  std::string basename_f;
  std::string sps_file_name;
  size_t total_number_levels_i;
  size_t total_number_levels_f;
  std::vector<std::vector<float>> blank_vector_i;
  std::vector<std::vector<float>> blank_vector_f;
  std::vector<MschemeOrbit> m_orbits;
  std::vector<int> jorbits;
  std::vector<int> Jlist_i;
  std::vector<int> Jlist_f;
  int MJtot_i;
  int MJtot_f;
  int Nshell;
  int A_i,Z_i,A_f,Z_f,Acore,Zcore;
  JBasis jbasis_i;
  JBasis jbasis_f;
  std::vector<NuVec> nuvec_list_i;
  std::vector<NuVec> nuvec_list_f;
  std::unordered_map< key_type, std::vector<std::vector<float>> > amplitudes_i;  // this is the moneymaker.
  std::unordered_map< key_type, std::vector<std::vector<float>> > amplitudes_f;
  std::unordered_map<int,int> max_states_per_J_i;
  std::unordered_map<int,int> max_states_per_J_f;
  std::vector<int> ket_a;
  std::vector<int> ket_b;
  std::vector<int> ket_J;
  static std::vector<char> an_code;   // only needed for file reading
  static std::vector<std::string> periodic_table;   // only needed for file reading
  Profiler profiler;
  std::string densfile_name;
  bool same_basename_f_i;


  TransitionDensity();
  TransitionDensity(std::vector<int> jlist);
  TransitionDensity(std::vector<int> jlist_i, std::vector<int> jlist_f);

  void ReadFiles( );
  void CalculateMschemeAmplitudes();
  void CalculateMschemeAmplitudes_fi(std::vector<NuVec>& , JBasis& , std::unordered_map<int,int>&, std::vector<std::vector<float>>&, std::unordered_map< key_type, std::vector<std::vector<float>> >& );
  void SetAZ(int a, int z){A_i=a;Z_i=z;A_f=a;Z_f=a;};
  void SetAZ_i(int a, int z){A_i=a;Z_i=z;};
  void SetAZ_f(int a, int z){A_f=a;Z_f=z;};
  void SetAZcore(int a, int z){Acore=a;Zcore=z;};
  double OBTD(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a, int m_index_b, int Lambda2 );
  double TBTD(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a, int m_index_b, int m_index_c, int m_index_d, int J2ab, int J2cd, int Lambda2 );
//  arma::mat GetOneBodyTransitionOperator( std::string filename, int& Lambda, int& RankT, int& parity );
//  arma::mat GetTwoBodyTransitionOperator( std::string filename, int& Lambda, int& RankT, int& parity );
//  void GetScalarTransitionOperator( std::string filename, double& Op0b, arma::mat& Op1b, arma::mat& Op2b );
  arma::mat CalcOBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2);
  arma::mat CalcTBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2);
  arma::mat ReadOBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, std::string fname);
  arma::mat ReadTBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, std::string fname);

  arma::vec CalcTransitionDensity_ax( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f);
  arma::mat CalcTransitionDensity_axaxa( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f);
  double TD_ax(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a );
  double TD_axaxa(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a, int m_index_b, int m_index_c, int J2ab );

  arma::vec GetDaggerOperator_ax( std::string filename );
  arma::mat GetDaggerOperator_axaxa( std::string filename );

  void SetupKets();
  void Jplus(std::vector<key_type>& mvecs_in, std::vector<double>& amp_in, int J2, int M2);
//  void WriteEGV( std::string fname);
//  void WriteTRDENS_input(std::string fname);
  void ReadDensities( );
  void GetAZFromFileName( );
  void ReadSPfile();
  void SetDensFile( std::string name ) ;


};



// Some useful operator overrides
std::vector<float> operator*(const float lhs, const std::vector<float>& rhs);
std::vector<float> operator*(const std::vector<float>& lhs, const std::vector<float>& rhs);
std::vector<float>& operator+=(std::vector<float>& lhs, const std::vector<float>& rhs);




#define TransitionDensity_h
#endif
