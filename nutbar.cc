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



// forward declaration of overloaded operator "+" so we can add fields of matrices
arma::field<arma::mat> operator+(arma::field<arma::mat>&,arma::field<arma::mat>&);




int main(int argc, char** argv)
{

  // TransitionDensity class does all the heavy lifting
  Profiler profiler;
  ReadWrite readwrite;
  TransitionDensity trans;

  
//  Settings settings;
  Settings& settings = readwrite.settings;
  
  if (argc>1)
  {
    std::ifstream infile(argv[1]);
    readwrite.ReadInput(infile, "batch");
  }
  else
  {
    readwrite.ReadInput(std::cin, "interactive");
  }
  

  // Write out the options we selected, as a sanity check.
  readwrite.PrintOptions();
  readwrite.PrintOperatorNames();
 

  // Read a bunch of NuShellX files, and put things into
  // the elaborate data structures contained in TransitionDensity
  readwrite.ReadNuShellFiles( trans );


  // if we're not reading the transition densities from file, we have to go ahead and calculate them.
  if ( not settings.densities_from_file  )
  {
    trans.CalculateMschemeAmplitudes();
  }

  if (readwrite.settings.write_egv)
  {
    readwrite.WriteEGV(trans, "mbpt.egv");
    readwrite.WriteTRDENS_input(trans, "trdens.in");
  }
  

  arma::mat TensorOp1b;
  arma::mat TensorOp2b;


  // Read in the operator files. This should also be handled by ReadWrite
  int Lambda=0,RankT=0,parity=0;
  if (settings.tensor_op_files.size() > 0)
  {
    TensorOp1b = readwrite.ReadOneBodyTransitionOperator(settings.tensor_op_files[0], trans, Lambda, RankT, parity);
  }
  if (settings.tensor_op_files.size() > 1)
  {
    int Lambda2,RankT2,parity2;
    TensorOp2b = readwrite.ReadTwoBodyTransitionOperator(settings.tensor_op_files[1], trans, Lambda2, RankT2, parity2);
    if (Lambda2 != Lambda or RankT2!=RankT or parity2 != parity)
    {
      std::cout << "ERROR Tensor rank mismatch between files " << settings.tensor_op_files[0] << "  and  " << settings.tensor_op_files[1] << std::endl;
      std::cout << " Lambda: " << Lambda << "," << Lambda2 << "   Tz: " << RankT << "," << RankT2 << "  parity: " << parity << "," << parity2 << std::endl;
      return 1;
    }
  }

  arma::vec DaggerOp_ax;
  arma::mat DaggerOp_axaxa;

  if (settings.dagger_op_files.size() > 0)
  {
    DaggerOp_ax    = trans.GetDaggerOperator_ax( settings.tensor_op_files[0] );
    DaggerOp_axaxa = trans.GetDaggerOperator_axaxa( settings.tensor_op_files[0] );
  }
  
  
//  std::vector<ScalarOperator> OpScalar(settings.scalar_op_file.size());
  std::vector<ScalarOperator> ScalarOps;
//  for (size_t i=0; i<settings.scalar_op_files.size();++i)
  std::cout << "Reading scalar ops" << std::endl;
  for (auto opfile : settings.scalar_op_files )
  {
    ScalarOps.push_back(  readwrite.ReadScalarOperator(opfile, trans) );
  }
  std::cout << "done. size of ScalarOps is " << ScalarOps.size() << std::endl;

//  std::vector< double>    OpScalar0b(settings.scalar_op_files.size());
//  std::vector< arma::mat> OpScalar1b(settings.scalar_op_files.size());
//  std::vector< arma::mat> OpScalar2b(settings.scalar_op_files.size());
  


  // construct a (nJi x nJf) matrix of many-body matrix elements for each Ji, Jf
  arma::field<arma::mat> TensorME1(settings.J2_f.size(),settings.J2_i.size() );
  arma::field<arma::mat> TensorME2(settings.J2_f.size(),settings.J2_i.size() );

  // Below, the vector corresponds to the different scalar operators,
  // the arma::field corresponds to the initial and final J 
  // and the arma::mat corresponds to which states for a given Ji, Jf
//  std::vector< arma::field<arma::mat> > ScalarME1(OpScalar1b.size());  
//  std::vector< arma::field<arma::mat> > ScalarME2(OpScalar1b.size());
  std::vector< arma::field<arma::mat> > ScalarME1(ScalarOps.size());  
  std::vector< arma::field<arma::mat> > ScalarME2(ScalarOps.size());

  for (size_t i=0; i<settings.scalar_op_files.size();++i)
  {
//      trans.GetScalarTransitionOperator(settings.scalar_op_files[i],OpScalar0b[i],OpScalar1b[i],OpScalar2b[i]);
//      readwrite.ReadScalarTransitionOperator(settings.scalar_op_files[i], trans, OpScalar0b[i],OpScalar1b[i],OpScalar2b[i]);
      ScalarME1[i] = arma::field<arma::mat>(settings.J2_f.size(),settings.J2_i.size() );
      ScalarME2[i] = arma::field<arma::mat>(settings.J2_f.size(),settings.J2_i.size() );
  }

  arma::field<arma::vec> DaggerME_ax( settings.J2_f.size(), settings.J2_i.size() );
  arma::field<arma::mat> DaggerME_axaxa( settings.J2_f.size(), settings.J2_i.size() );
 



//  trans.SetupKets();
//  std::string densityfilename = "nutbar_densities_" + settings.basename_vectors_f.substr( settings.basename_vectors_f.find_last_of("/")+1 ) + ".dat";
//  std::string logfilename = "nutbar_" + settings.basename_vectors_f.substr( settings.basename_vectors_f.find_last_of("/")+1 ) + ".log";
//  std::ofstream logfile;
//  if (settings.write_log) logfile.open(logfilename);

//  if ( settings.densities_from_file )
//    trans.SetDensFile(densityfilename);

  if ( not readwrite.settings.densities_from_file )
  {
    readwrite.WriteDensityHeader( trans );
  }

 
  for (size_t indexJi=0;indexJi<settings.J2_i.size();++indexJi )
  {
   for (size_t indexJf=0;indexJf<settings.J2_f.size();++indexJf )
   {

     if (settings.dagger_op_files.size() > 0)
     {
       DaggerME_ax(indexJf,indexJi).zeros(settings.NJ_f[indexJf],settings.NJ_i[indexJi]);
       DaggerME_axaxa(indexJf,indexJi).zeros(settings.NJ_f[indexJf],settings.NJ_i[indexJi]);
     }
     if ( std::abs(settings.J2_i[indexJi] - settings.J2_f[indexJf])> 2*Lambda) continue;
     if (settings.tensor_op_files.size() > 0)
     {
       TensorME1(indexJf,indexJi).zeros(settings.NJ_f[indexJf],settings.NJ_i[indexJi]);
       TensorME2(indexJf,indexJi).zeros(settings.NJ_f[indexJf],settings.NJ_i[indexJi]);
     }

     if (settings.J2_f[indexJf]==settings.J2_i[indexJi])
     {
//       for (size_t isc=0; isc<OpScalar1b.size();++isc)
       for (size_t isc=0; isc<ScalarOps.size();++isc)
       {
          ScalarME1[isc](indexJf,indexJi).zeros(settings.NJ_f[indexJf],settings.NJ_i[indexJi]);
          ScalarME2[isc](indexJf,indexJi).zeros(settings.NJ_f[indexJf],settings.NJ_i[indexJi]);
       }
     }
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
         obtd = readwrite.ReadOBTD(indexJi,ivec,indexJf,fvec,Lambda*2, trans);
         tbtd = readwrite.ReadTBTD(indexJi,ivec,indexJf,fvec,Lambda*2, trans);
//         obtd = trans.ReadOBTD(indexJi,ivec,indexJf,fvec,Lambda*2,densityfilename);
//         tbtd = trans.ReadTBTD(indexJi,ivec,indexJf,fvec,Lambda*2,densityfilename);
       }
       else
       {
         obtd = trans.CalcOBTD(indexJi,ivec,indexJf,fvec,Lambda*2);
         readwrite.WriteOBTD(indexJi,ivec,indexJf,fvec,Lambda*2, obtd);
         if (settings.tensor_op_files.size() > 1 or settings.scalar_op_files.size()>0)
         {
           tbtd = trans.CalcTBTD(indexJi,ivec,indexJf,fvec,Lambda*2);
           readwrite.WriteTBTD(indexJi,ivec,indexJf,fvec,Lambda*2, tbtd);
         }
         if (settings.dagger_op_files.size() > 0 )
         {
           td_ax = trans.CalcTransitionDensity_ax( indexJi, ivec, indexJf, fvec);
           td_axaxa = trans.CalcTransitionDensity_axaxa( indexJi, ivec, indexJf, fvec);
         }
       }
 
       if (settings.tensor_op_files.size()>0)
       {
         double obme = arma::accu( TensorOp1b % obtd );
         TensorME1(indexJf,indexJi)(fvec,ivec) = obme;

         if (settings.write_log)
         {
           readwrite.WriteLog_Tensor1b( indexJi, ivec, indexJf, fvec, Lambda, TensorOp1b, obtd );
         }

       }
       
       if (settings.tensor_op_files.size()>1)
       {
         double tbme = arma::accu( TensorOp2b % tbtd);
         TensorME2(indexJf,indexJi)(fvec,ivec) = tbme;

         if (settings.write_log)
         {
           readwrite.WriteLog_Tensor2b( Lambda, TensorOp2b, tbtd );
         }

       }
  
//       if (settings.J2_i[indexJi]==settings.J2_f[indexJf] and OpScalar1b.size()>0)
       if (settings.J2_i[indexJi]==settings.J2_f[indexJf])
       {
//         for (size_t iop=0; iop<OpScalar1b.size(); ++iop)
         for (size_t iop=0; iop<ScalarOps.size(); ++iop)
         {
//           double scop1 = arma::accu( OpScalar1b[iop] % obtd ) ;
//           double scop2 = arma::accu( OpScalar2b[iop] % tbtd ) ;
           double scop1 = arma::accu( ScalarOps[iop].OneBody % obtd ) ;
           double scop2 = arma::accu( ScalarOps[iop].TwoBody % tbtd ) ;
           ScalarME1[iop](indexJf,indexJi)(fvec,ivec) = scop1;
           ScalarME2[iop](indexJf,indexJi)(fvec,ivec) = scop2;
         }
       }
      }
    }
   }
  }

  std::ostringstream oss;
  if ( settings.tensor_op_files.size()>0)
  {
    oss << Lambda;
    std::string filename = "nutbar_tensor" + oss.str() + "_" + settings.basename_vectors_f.substr( settings.basename_vectors_f.find_last_of("/")+1 );
    if (settings.basename_vectors_i != settings.basename_vectors_f)
    {
       filename +  "_" + settings.basename_vectors_i.substr( settings.basename_vectors_i.find_last_of("/")+1 );
    }
    filename += ".dat";
    std::ofstream tensor_out(filename);
  
    tensor_out << "##########################################################################################" << std::endl;
    tensor_out << "# initial vector basename: " << settings.basename_vectors_i << std::endl;
    tensor_out << "# final   vector basename: " << settings.basename_vectors_f << std::endl;
    tensor_out << "# One body file: " << settings.tensor_op_files[0] << std::endl;
    if (settings.tensor_op_files.size()>1)
      tensor_out << "# Two body file: " << settings.tensor_op_files[1] << std::endl;
    tensor_out << "# Jf  nJf     Ji  nJi      Ef          Ei         <Op1b>          <Op2b>         <Op1b+2b> " << std::endl;
    tensor_out << "###########################################################################################" << std::endl;
  
    for (size_t indexJf=0; indexJf<settings.J2_f.size();++indexJf)
    {
      double Jf = settings.J2_f[indexJf]*0.5;
      for (size_t njf=0; njf<settings.NJ_f[indexJf];++njf)
      {
        for (size_t indexJi=0; indexJi<settings.J2_i.size();++indexJi)
        {
         double Ji = settings.J2_i[indexJi]*0.5;
         for (size_t nji=0; nji<settings.NJ_i[indexJi];++nji)
         {
           tensor_out << std::fixed << std::setw(4) << std::setprecision(1) << Jf << " " << std::setw(3) << njf+1 << "    "
                      << std::fixed << std::setw(4) << std::setprecision(1) << Ji << " " << std::setw(3) << nji+1 << " "
                      << std::fixed << std::setw(10) << std::setprecision(3) << trans.nuvec_list_f[indexJf].alpha[njf] << "  "
                      << std::fixed << std::setw(10) << std::setprecision(3) << trans.nuvec_list_i[indexJi].alpha[nji] << "  "
                      << std::scientific << std::setw(14) << std::setprecision(6) << TensorME1(indexJf,indexJi)(njf,nji) << "  "
                      << std::scientific << std::setw(14) << std::setprecision(6) << TensorME2(indexJf,indexJi)(njf,nji) << "  "
                      << std::scientific << std::setw(14) << std::setprecision(6) << TensorME1(indexJf,indexJi)(njf,nji) + TensorME2(indexJf,indexJi)(njf,nji) << "  "
                      << std::endl;
         }
        }
      }
    }
  }


  std::cout << "WriteScalarResults..." << std::endl;
  readwrite.WriteScalarResults( trans, ScalarOps, ScalarME1, ScalarME2 );
//  for (size_t isc=0; isc<settings.scalar_op_files.size(); ++isc)
//  for (size_t isc=0; isc<ScalarOps.size(); ++isc)
//  {
//    oss.str("");
//    oss << isc;
//    std::string filename = "nutbar_scalar" + oss.str() + "_" + settings.basename_vectors_f.substr( settings.basename_vectors_f.find_last_of("/")+1 );
//    if (settings.basename_vectors_i != settings.basename_vectors_f)
//    {
//       filename +  "_" + settings.basename_vectors_i.substr( settings.basename_vectors_i.find_last_of("/")+1 );
//    }
//    filename += ".dat";
//    std::ofstream scalar_out(filename); 
//
//    scalar_out << "############################################################################################################" << std::endl;
//    scalar_out << "# initial vector basename: " << settings.basename_vectors_i << std::endl;
//    scalar_out << "# final   vector basename: " << settings.basename_vectors_f << std::endl;
//    scalar_out << "# Operator file: " << settings.scalar_op_files[isc] << std::endl;
////    scalar_out << "# Zero body term: " << OpScalar0b[isc] << std::endl;
//    scalar_out << "# Zero body term: " << ScalarOps[isc].ZeroBody << std::endl;
//    scalar_out << "# Jf  nJf     Ji  nJi      Ei          Ef         <Op1b>          <Op2b>         <Op1b+2b>     <Op0b+1b+2b>" << std::endl;
//    scalar_out << "###########################################################################################################" << std::endl;
//  
//    for (size_t indexJf=0; indexJf<settings.J2_f.size();++indexJf)
//    {
//      double Jf = settings.J2_f[indexJf]*0.5;
//      for (size_t njf=0; njf<settings.NJ_f[indexJf];++njf)
//      {
//        for (size_t indexJi=0; indexJi<settings.J2_i.size();++indexJi)
//        {
//         double Ji = settings.J2_i[indexJi]*0.5;
//         if (std::abs(Ji-Jf)>0.1) continue;
//         for (size_t nji=0; nji<settings.NJ_i[indexJi];++nji)
//         {
//           if ( ( std::find( begin(settings.options), end(settings.options), "diag") != end(settings.options)) and not( Ji==Jf and nji==njf))  continue;
//           double zerobody = (nji==njf) ? ScalarOps[isc].ZeroBody : 0;
//           scalar_out << std::fixed << std::setw(4) << std::setprecision(1) << Jf << " " << std::setw(3) << njf+1 << "    "
//                      << std::fixed << std::setw(4) << std::setprecision(1) << Ji << " " << std::setw(3) << nji+1 << " "
//                      << std::fixed << std::setw(10) << std::setprecision(3) << trans.nuvec_list_f[indexJf].alpha[njf] << "  "
//                      << std::fixed << std::setw(10) << std::setprecision(3) << trans.nuvec_list_i[indexJi].alpha[nji] << "  "
//                      << std::scientific << std::setw(14) << std::setprecision(6) << ScalarME1[isc](indexJf,indexJi)(njf,nji) << "  "
//                      << std::scientific << std::setw(14) << std::setprecision(6) << ScalarME2[isc](indexJf,indexJi)(njf,nji) << "  "
//                      << std::scientific << std::setw(14) << std::setprecision(6) << ScalarME1[isc](indexJf,indexJi)(njf,nji) + ScalarME2[isc](indexJf,indexJi)(njf,nji) << "  "
//                      << std::scientific << std::setw(14) << std::setprecision(6) << ScalarME1[isc](indexJf,indexJi)(njf,nji) + ScalarME2[isc](indexJf,indexJi)(njf,nji) + zerobody 
//                      << std::endl;
//         }
//        }
//      }
//    }
//
//
//
//  }



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
