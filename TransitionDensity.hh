
#ifndef TransitionDensity_h

#include "NuBasis.hh"
#include "NuProj.hh"
#include "NuVec.hh"
#include "JMState.hh"
#include "JBasis.hh"

#include <vector>
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


  TransitionDensity();
  TransitionDensity(vector<int> jlist);
  void ReadInputFromFile(string filename);
  void ReadInputInteractive();
  void ReadFiles( );
//  void ReadFiles( string spsfile, string abfile_base );
  void CalculateMschemeAmplitudes();
  void SetAZ(int a, int z){A=a;Z=z;};
  void SetAZcore(int a, int z){Acore=a;Zcore=z;};
  double OBTD(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a, int m_index_b, int Lambda2 );
  double TBTD(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a, int m_index_b, int m_index_c, int m_index_d, int J2ab, int J2cd, int Lambda2 );
  void WriteEGV( string fname);
  void WriteTRDENS_input(string fname);
  void SetMaxStatesPerJ( int J2, int imax){ max_states_per_J[J2] = imax;};
  void GetAZFromFileName( );


};



// Some useful operator overrides
vector<float> operator*(const float lhs, const vector<float>& rhs);
vector<float> operator*(const vector<float>& lhs, const vector<float>& rhs);
vector<float>& operator+=(vector<float>& lhs, const vector<float>& rhs);




#define TransitionDensity_h
#endif
