
#ifndef Operators_hh
#define Operators_hh 1


#include <armadillo>
#include "Settings.hh"


struct ScalarOperator
{
  double ZeroBody;
  arma::mat OneBody;
  arma::mat TwoBody;
};


struct TensorOperator
{
  int Lambda;
  int Rank_T;
  int Parity;
  arma::mat OneBody;
  arma::mat TwoBody;
  TensorOperator() : Lambda(0),Rank_T(0),Parity(0) {};
};

struct DaggerOperator
{
  int Lambda2;
  int Rank_T2;
  int Parity;
  arma::mat Op_ax;
  arma::mat Op_axaxa;
};



struct ScalarNME
{
  double ZeroBody;
  arma::field<arma::mat> OneBody;
  arma::field<arma::mat> TwoBody;

  ScalarNME(Settings& settings);
};

struct TensorNME
{
  double ZeroBody;
  arma::field<arma::mat> OneBody;
  arma::field<arma::mat> TwoBody;
  int Lambda;

  TensorNME(Settings& settings, TensorOperator& TensorOp);
};


struct DaggerNME
{
  arma::field<arma::mat> ax;
  arma::field<arma::mat> axaxa;
  int Lambda2;

  DaggerNME(Settings& settings, DaggerOperator& DaggerOp);
};

#endif
