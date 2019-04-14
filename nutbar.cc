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
 

  std::cout << "Reading Nushell Files" << std::endl;
  // Read a bunch of NuShellX files, and put things into
  // the elaborate data structures contained in TransitionDensity
  readwrite.ReadNuShellFiles( trans );


  std::cout << "Calculating Mscheme Amplitudes" << std::endl;
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
  

  std::cout << "Read in Operators... "<< std::endl;

  std::vector<ScalarOperator> ScalarOps;
  std::vector<ScalarNME> scalarnme;
  for (auto opfile : settings.scalar_op_files )
  {
    ScalarOps.push_back(  readwrite.ReadScalarOperator(opfile) );
    scalarnme.push_back( ScalarNME(settings) );
  }


  TensorOperator TensorOp = readwrite.ReadTensorOperator( settings.tensor_op_files);
  TensorNME tensornme(settings, TensorOp);

  if ( ScalarOps.size()>0 and settings.tensor_op_files.size()>0)
  {
    std::cout << "TROUBLE!! I'm not smart enough to handle both a scalar and a tensor at the same time. Pathetic, I know...  Bailing out." << std::endl;
    return 1;
  }

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
    readwrite.WriteDensityHeader( );
  }
  if (settings.write_log)
  {
    readwrite.WriteLogHeader();
  }

  std::cout << "Calculate densities" << std::endl;

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
       std::vector<arma::vec> td_ax(DaggerOps.size()); // a+ operator
       std::vector<arma::mat> td_axaxa(DaggerOps.size()); // a+a+a operator
       if ( settings.densities_from_file )
       {
         std::cout << "Reading densities from file..." << std::endl;
         obtd = readwrite.ReadOBTD(indexJi,ivec,indexJf,fvec,TensorOp.Lambda*2, trans);
         tbtd = readwrite.ReadTBTD(indexJi,ivec,indexJf,fvec,TensorOp.Lambda*2, trans);
       }
       else  // if we're not reading from file, we need to calculate the densities
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
//           for ( auto& dagnme : daggernme )
           for ( size_t idag=0; idag<DaggerOps.size();idag++ )
           {
             td_ax[idag] = trans.CalcTransitionDensity_ax( indexJi, ivec, indexJf, fvec, DaggerOps[idag].Lambda2, settings);
             td_axaxa[idag] = trans.CalcTransitionDensity_axaxa( indexJi, ivec, indexJf, fvec, DaggerOps[idag].Lambda2, settings);

           }
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
           // We divide by sqrt(2J+1) because for a scalar operator we typically expect to be given < f | Op | i >, rather than <f || Op || i >.
           // I'm not 100% sure this is the best way to go, but it's what we're doing for now.
           scalarnme[iop].OneBody(indexJf,indexJi)(fvec,ivec) = arma::accu( ScalarOps[iop].OneBody % obtd ) / sqrt( settings.J2_i[indexJi] +1.0 );
           scalarnme[iop].TwoBody(indexJf,indexJi)(fvec,ivec) = arma::accu( ScalarOps[iop].TwoBody % tbtd ) / sqrt( settings.J2_i[indexJi] +1.0 );
           if (settings.write_log)
           {
             readwrite.WriteLog_Scalar1b( indexJi, ivec, indexJf, fvec,  ScalarOps[iop], obtd );
             readwrite.WriteLog_Scalar2b( ScalarOps[iop], tbtd );
           }
         }
       }

       for ( size_t idag=0; idag<DaggerOps.size(); ++idag)
       {
             if (settings.write_log)
             {
               readwrite.WriteLog_Dagger_ax( indexJi, ivec, indexJf, fvec,  DaggerOps[idag], td_ax[idag] );
               readwrite.WriteLog_Dagger_axaxa( DaggerOps[idag], td_axaxa[idag] );
             }

             daggernme[idag].ax(indexJf,indexJi)(fvec,ivec) = arma::accu( td_ax[idag] % DaggerOps[idag].Op_ax );
             daggernme[idag].axaxa(indexJf,indexJi)(fvec,ivec) = arma::accu( td_axaxa[idag] % DaggerOps[idag].Op_axaxa );
       }

      } // for fvec
    } // for ivec
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


  if (DaggerOps.size()>0)
  {
    std::cout << "WriteDaggerResults..." << std::endl;
    readwrite.WriteDaggerResults( trans, DaggerOps, daggernme );
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
