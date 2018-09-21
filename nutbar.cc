//
//  Program nutbar  (NuShellX Transtitions from Binary Arrays)
//
//

#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <string>
#ifndef NOBOOST
#include <boost/filesystem.hpp>
#endif

#include "NuBasis.hh"
#include "NuProj.hh"
#include "NuVec.hh"
#include "JMState.hh"
#include "JBasis.hh"
#include "TransitionDensity.hh"
#include "Profiler.hh"
#include "ReadWrite.hh"
#include "Operators.hh"
#include "Settings.hh"



// forward declaration of overloaded operator "+" so we can add fields of matrices
arma::field<arma::mat> operator+(arma::field<arma::mat>&,arma::field<arma::mat>&);




int main(int argc, char** argv)
{

  // TransitionDensity class does all the heavy lifting
  Profiler profiler;
  ReadWrite readwrite;
  TransitionDensity trans;

  
  
  if (argc>1)
  {
    std::ifstream infile(argv[1]);
    readwrite.ReadInput(infile, "batch");
  }
  else
  {
    readwrite.ReadInput(std::cin, "interactive");
  }

  Settings& settings = readwrite.settings;
  

  // Write out the options we selected, as a sanity check.
  readwrite.PrintOptions();
  readwrite.PrintOperatorNames();
 

  // Read a bunch of NuShellX files, and put things into
  // the elaborate data structures contained in TransitionDensity
  readwrite.ReadNuShellFiles( trans );


  // if we're not reading the transition densities from file, we have to go ahead and calculate them.
  if ( not settings.densities_from_file  )
  {
    trans.CalculateMschemeAmplitudes(settings);
  }

  if (readwrite.settings.write_egv)
  {
    readwrite.WriteEGV(trans, "mbpt.egv");
    readwrite.WriteTRDENS_input(trans, "trdens.in");
  }
  



  std::vector<ScalarOperator> ScalarOps;
  std::vector<ScalarNME> scalarnme;
  for (auto opfile : settings.scalar_op_files )
  {
    ScalarOps.push_back(  readwrite.ReadScalarOperator(opfile) );
    scalarnme.push_back( ScalarNME(settings) );
  }


  TensorOperator TensorOp = readwrite.ReadTensorOperator( settings.tensor_op_files);
  TensorNME tensornme(settings, TensorOp);


  std::vector<DaggerOperator> DaggerOps;
  std::vector<DaggerNME> daggernme;
  for ( auto opfile : settings.dagger_op_files )
  {
    DaggerOps.push_back( readwrite.ReadDaggerOperator( opfile) );
    daggernme.push_back( DaggerNME(settings, DaggerOps.back()) );
  }


  // construct a (nJi x nJf) matrix of many-body matrix elements for each Ji, Jf


  if ( not readwrite.settings.densities_from_file )
  {
    readwrite.WriteDensityHeader( trans );
  }


  for (size_t indexJi=0;indexJi<settings.J2_i.size();++indexJi )
  {
   for (size_t indexJf=0;indexJf<settings.J2_f.size();++indexJf )
   {
     for (int ivec = 0; ivec<settings.NJ_i[indexJi]; ++ivec)
     {
      for (int fvec = 0; fvec<settings.NJ_f[indexJf]; ++fvec)
      {
       if (settings.basename_vectors_i == settings.basename_vectors_f)
       {
        if ( settings.J2_i[indexJi]==settings.J2_f[indexJf] and fvec>ivec) continue;
        if ( (settings.diagonal_only) and ( (settings.J2_i[indexJi]!=settings.J2_f[indexJf]) or (ivec!=fvec)) )  continue;
       }
       arma::mat obtd,tbtd;
       arma::vec td_ax; // a+ operator
       arma::mat td_axaxa; // a+a+a operator
       if ( settings.densities_from_file )
       {
         std::cout << "Reading densities from file..." << std::endl;
         obtd = readwrite.ReadOBTD(indexJi,ivec,indexJf,fvec,TensorOp.Lambda*2, trans);
         tbtd = readwrite.ReadTBTD(indexJi,ivec,indexJf,fvec,TensorOp.Lambda*2, trans);
       }
       else
       {
         obtd = trans.CalcOBTD(indexJi,ivec,indexJf,fvec,TensorOp.Lambda*2,settings);
         readwrite.WriteOBTD(indexJi,ivec,indexJf,fvec,TensorOp.Lambda*2, obtd);
         if (settings.tensor_op_files.size() > 1 or settings.scalar_op_files.size()>0)
         {
           tbtd = trans.CalcTBTD(indexJi,ivec,indexJf,fvec,TensorOp.Lambda*2,settings);
           readwrite.WriteTBTD(indexJi,ivec,indexJf,fvec,TensorOp.Lambda*2, tbtd);
         }
         if (settings.dagger_op_files.size() > 0 )
         {
           td_ax = trans.CalcTransitionDensity_ax( indexJi, ivec, indexJf, fvec, settings);
           td_axaxa = trans.CalcTransitionDensity_axaxa( indexJi, ivec, indexJf, fvec, settings);
         }
       }
 
       if (settings.tensor_op_files.size()>0)
       {
         tensornme.OneBody(indexJf,indexJi)(fvec,ivec) = arma::accu( TensorOp.OneBody % obtd );

         if (settings.write_log)
         {
           readwrite.WriteLog_Tensor1b( indexJi, ivec, indexJf, fvec, TensorOp, obtd );
         }

       }
       
       if (settings.tensor_op_files.size()>1)
       {
         tensornme.TwoBody(indexJf,indexJi)(fvec,ivec) = arma::accu( TensorOp.TwoBody % tbtd );

         if (settings.write_log)
         {
           readwrite.WriteLog_Tensor2b( TensorOp, tbtd );
         }

       }
  
       if (settings.J2_i[indexJi]==settings.J2_f[indexJf])
       {
         for (size_t iop=0; iop<ScalarOps.size(); ++iop)
         {
           scalarnme[iop].OneBody(indexJf,indexJi)(fvec,ivec) = arma::accu( ScalarOps[iop].OneBody % obtd ) ;
           scalarnme[iop].TwoBody(indexJf,indexJi)(fvec,ivec) = arma::accu( ScalarOps[iop].TwoBody % tbtd ) ;
         }
       }
      }
    }
   }
  }

  if (settings.tensor_op_files.size()>0)
  {
    std::cout << "Writing Tensor files..." << std::endl;
    readwrite.WriteTensorResults( trans, TensorOp, tensornme );
  }




  if (settings.scalar_op_files.size()>0)
  {
    std::cout << "WriteScalarResults..." << std::endl;
    readwrite.WriteScalarResults( trans, ScalarOps, scalarnme );
  }




  profiler.PrintAll();

  return 0;
}







// straightforward element-by-element addition
arma::field<arma::mat> operator+(arma::field<arma::mat>& lhs ,arma::field<arma::mat>& rhs)
{
  arma::field<arma::mat> out = lhs;
  for (size_t irow=0; irow<lhs.n_rows; ++irow)
  {
    for (size_t icol=0; icol<lhs.n_cols; ++icol)
    {
      out(irow,icol) += rhs(irow,icol);
    }
  }
  return out;
}
