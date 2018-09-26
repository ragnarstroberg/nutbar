#ifndef ReadWrite_hh
#define ReadWrite_hh 1

#include <tuple>
#include <vector>
#include <string>
#include <fstream>

#include "TransitionDensity.hh"
#include "Operators.hh"
#include "Settings.hh"




class ReadWrite
{
 public:

  Settings settings;
  static std::vector<char> an_code;
  static std::vector<std::string> periodic_table;

  int Acore,Zcore,A_i,Z_i,A_f,Z_f;
  std::ofstream logfile;
  std::fstream densityfile;


  ReadWrite() : Acore(0),Zcore(0), A_i(0), Z_i(0), A_f(0), Z_f(0) {};

  void GetCoreFromSPfile();
  void GetAZFromFileName();
  std::vector<std::string> FindNBAFiles( std::string basename, int Zval, int Nval, std::string a_or_b );
  std::string FindXVCFile( std::string basename, int Jtot, int Zval, int Nval );
  void ReadInput(std::istream& input, std::string mode);
  void ReadNuShellFiles( TransitionDensity& trans );

//  ScalarOperator ReadScalarOperator( std::string filename, TransitionDensity& trans);
//  TensorOperator ReadTensorOperator( std::vector<std::string> filenames, TransitionDensity& trans);
//  DaggerOperator ReadDaggerOperator( std::string filename, TransitionDensity& trans);
  ScalarOperator ReadScalarOperator( std::string filename);
  TensorOperator ReadTensorOperator( std::vector<std::string> filenames);
  DaggerOperator ReadDaggerOperator( std::string filename);

//  void ReadScalarTransitionOperator( std::string filename, TransitionDensity& trans, double& Op0b, arma::mat& Op1b, arma::mat& Op2b);
//  arma::mat ReadOneBodyTransitionOperator( std::string filename, TransitionDensity& trans, int& Rank_J, int& Rank_T, int& parity );
//  arma::mat ReadTwoBodyTransitionOperator( std::string filename, TransitionDensity& trans, int& Rank_J, int& Rank_T, int& parity );



  arma::mat ReadOBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, TransitionDensity& trans);
  arma::mat ReadTBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, TransitionDensity& trans);

  void WriteDensityHeader( );
  void WriteOBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, arma::mat& obtd);
  void WriteTBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, arma::mat& tbtd);


  void WriteLogHeader();
  void WriteLog_Tensor1b( int indexJi, int ivec, int indexJf, int fvec, TensorOperator& TensorOp,  arma::mat& obtd );
  void WriteLog_Tensor2b( TensorOperator& TensorOp , arma::mat& tbtd );
  void WriteLog_Dagger_ax( int indexJi, int ivec, int indexJf, int fvec, DaggerOperator& DaggerOp, arma::vec& td_ax );
  void WriteLog_Dagger_axaxa( DaggerOperator& DaggerOp, arma::mat& td_axaxa );
//  void WriteLog_Tensor1b( int indexJi, int ivec, int indexJf, int fvec, int Lambda, arma::mat& TensorOp1b, arma::mat& obtd );
//  void WriteLog_Tensor2b( int Lambda, arma::mat& TensorOp2b, arma::mat& tbtd );

//  void WriteScalarResults( TransitionDensity& trans, std::vector<ScalarOperator>& ScalarOps, std::vector<arma::field<arma::mat>>& ScalarME1, std::vector<arma::field<arma::mat>>& ScalarME2 );
//  void WriteTensorResults( TransitionDensity& trans, TensorOperator& TensorOp, arma::field<arma::mat>& TensorME1, arma::field<arma::mat>& TensorME2 );
  void WriteScalarResults( TransitionDensity& trans, std::vector<ScalarOperator>& ScalarOps, std::vector<ScalarNME>& scalarnme );
  void WriteTensorResults( TransitionDensity& trans, TensorOperator& TensorOp, TensorNME& tensornme );
  void WriteDaggerResults( TransitionDensity& trans, std::vector<DaggerOperator>& DaggerOp, std::vector<DaggerNME>& daggernme );

  void WriteEGV(TransitionDensity& trans, std::string fname);
  void WriteTRDENS_input( TransitionDensity& trans, std::string fname);

  void PrintOptions();
  void PrintOperatorNames();


};



#endif
