#ifndef Settings_hh
#define Settings_hh


#include <vector>
#include <string>

//#include "NuBasis.hh"
#include "JBasis.hh"





// Define a struct to hold all the input information
// and some basic info for the calculation
struct Settings
{
  std::string basename_vectors_i; // base name for the vectors, e.g. ne200
  std::string basename_vectors_f; // base name for the vectors, e.g. ne200
  std::string basename_sps;     // base name for the sp and sps files, e.g. sdpn
  std::vector<std::string> scalar_op_files; // list of scalar operator file names
  std::vector<std::string> tensor_op_files; // list of tensor operator file names
  std::vector<std::string> dagger_op_files; // list of dagger operator file names
  std::vector<int> J2_i;  // 2*J for initial states
  std::vector<int> NJ_i;  // number of states for each initial J
  std::vector<int> J2_f;  // 2*J for final states
  std::vector<int> NJ_f;  // number of states for each final J
  std::vector<std::string> options;
  bool densities_from_file;
  bool write_egv;
  bool write_log;
  bool diagonal_only;
  bool same_basename_fi;
  bool same_states_fi;


  int Acore,Zcore;
  int A_i, Z_i, N_i;
  int A_f, Z_f, N_f;
  int total_number_levels_i;
  int total_number_levels_f;

  std::vector<int> ket_a;
  std::vector<int> ket_b;
  std::vector<int> ket_J;


  JBasis jbasis_i;
  JBasis jbasis_f;
  std::vector<MschemeOrbit> m_orbits;
  std::vector<int> jorbits;

  Settings() : densities_from_file(false), write_egv(false), write_log(false),
               same_basename_fi(false), same_states_fi(false),
               Acore(0),Zcore(0), A_i(0), Z_i(0), A_f(0), Z_f(0),
               total_number_levels_i(0), total_number_levels_f(0) 
    {};

};




#endif

