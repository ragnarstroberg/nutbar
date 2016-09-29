
#ifndef TransitionDensity_h

#include "NuBasis.hh"
#include "NuProj.hh"
#include "NuVec.hh"
#include "JMState.hh"
#include "JBasis.hh"

#include <vector>
#include <armadillo>
#include <unordered_map>

using namespace std;



class TransitionDensity
{
 public:
  
  string basename;
  string sps_file_name;
  size_t total_number_levels;
  vector<vector<float>> blank_vector;
  vector<MschemeOrbit> m_orbits;
  vector<int> jorbits;
  vector<int> Jlist;
  int MJtot;
  int Nshell;
  int A,Z,Acore,Zcore;
  vector<JBasis> jbasis_list;
  vector<NuVec> nuvec_list;
  unordered_map< vector<mvec_type>, vector<vector<float>>, KeyHash > amplitudes;
  static vector<char> an_code;
  static vector<string> periodic_table;
  unordered_map<int,int> max_states_per_J;
  vector<int> ket_a;
  vector<int> ket_b;
  vector<int> ket_J;


  TransitionDensity();
  TransitionDensity(vector<int> jlist);
  void ReadInputFromFile(string filename);
  void ReadInputInteractive();
  void ReadFiles( );
  void CalculateMschemeAmplitudes();
  void SetAZ(int a, int z){A=a;Z=z;};
  void SetAZcore(int a, int z){Acore=a;Zcore=z;};
  double OBTD(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a, int m_index_b, int Lambda2 );
  double TBTD(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a, int m_index_b, int m_index_c, int m_index_d, int J2ab, int J2cd, int Lambda2 );
  arma::mat GetOneBodyTransitionOperator( string filename, int& Lambda, int& RankT, int& parity );
  arma::mat GetTwoBodyTransitionOperator( string filename, int& Lambda, int& RankT, int& parity );
  void GetScalarTransitionOperator( string filename, double& Op0b, arma::mat& Op1b, arma::mat& Op2b );
  arma::mat CalcOBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2);
  arma::mat CalcTBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2);
  void SetupKets();
  void Jplus(vector<vector<mvec_type>>& mvecs_in, vector<double>& amp_in, int J2, int M2);
  void WriteEGV( string fname);
  void WriteTRDENS_input(string fname);
  void SetMaxStatesPerJ( int J2, int imax){ max_states_per_J[J2] = imax;};
  void GetAZFromFileName( );
  void ReadSPfile();


};



// Some useful operator overrides
vector<float> operator*(const float lhs, const vector<float>& rhs);
vector<float> operator*(const vector<float>& lhs, const vector<float>& rhs);
vector<float>& operator+=(vector<float>& lhs, const vector<float>& rhs);




#define TransitionDensity_h
#endif
