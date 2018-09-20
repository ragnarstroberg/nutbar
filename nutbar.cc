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



// forward declaration so we can add fields of matrices
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
//    settings = ReadInput(infile, "batch");
    readwrite.ReadInput(infile, "batch");
  }
  else
  {
//    settings = ReadInput(std::cin, "interactive");
    readwrite.ReadInput(std::cin, "interactive");
  }
  
  std::cout << "options: [ ";
  for (std::string opt : settings.options) std::cout << opt << " ";
  std::cout << " ]" << std::endl; 
 

  // Read a bunch of NuShellX files, and put things into
  // the elaborate data structures contained in TransitionDensity
  readwrite.ReadNuShellFiles( trans );




  // Some shuffling around to put things in the TransitionDensity class format
  // this should really be done in a much neater way
  // if all goes to plan, this can be eliminated.
//  trans.basename_i = settings.basename_vectors_i;
//  trans.basename_f = settings.basename_vectors_f;
//  trans.sps_file_name = settings.basename_sps + ".sps";
//  trans.Jlist_i = settings.J2_i;
//  trans.Jlist_f = settings.J2_f;
//  for (size_t i=0;i<settings.J2_i.size();++i)
//  {
//    trans.max_states_per_J_i[settings.J2_i[i]] = settings.NJ_i[i];
//  }
//  for (size_t i=0;i<settings.J2_f.size();++i)
//  {
//    trans.max_states_per_J_f[settings.J2_f[i]] = settings.NJ_f[i];
//  }


//  trans.ReadFiles();
  
  // if we're not reading the transition densities from file, we have to go ahead and calculate them.
//  if ( std::find( begin(settings.options), end(settings.options), "read_dens") == end(settings.options)  )
  if ( not settings.densities_from_file  )
  {
    trans.CalculateMschemeAmplitudes();
  }

//#ifdef VERBOSE
//    std::cout << "done calculating Mscheme amplitudes." << std::endl;
//    std::cout << "number of initial m-scheme states: " << trans.amplitudes_i.size() << " x " << trans.blank_vector_i.size()<< std::endl;
//    std::cout << "number of final m-scheme states: " << trans.amplitudes_f.size() << " x " << trans.blank_vector_f.size()<< std::endl;
//#endif
  

  if (settings.write_egv)
  {
    readwrite.WriteEGV(trans, "mbpt.egv");
    readwrite.WriteTRDENS_input(trans, "trdens.in");
//    trans.WriteEGV("mbpt.egv");
//    trans.WriteTRDENS_input("trdens.in");
  }
  

  std::cout << "operators: [ " ;
  for (auto ops : settings.scalar_op_files) std::cout << ops << " ";
  for (auto ops : settings.tensor_op_files) std::cout << ops << " ";
  for (auto ops : settings.dagger_op_files) std::cout << ops << " ";
  std::cout << "] " << std::endl;

  arma::mat TensorOp1b;
  arma::mat TensorOp2b;


  // Read in the operator files. This should also be handled by ReadWrite
  int Lambda=0,RankT=0,parity=0;
  if (settings.tensor_op_files.size() > 0)
  {
//    TensorOp1b = trans.GetOneBodyTransitionOperator(settings.tensor_op_files[0], Lambda, RankT, parity);
    TensorOp1b = readwrite.ReadOneBodyTransitionOperator(settings.tensor_op_files[0], trans, Lambda, RankT, parity);
  }
  if (settings.tensor_op_files.size() > 1)
  {
    int Lambda2,RankT2,parity2;
//    TensorOp2b = trans.GetTwoBodyTransitionOperator(settings.tensor_op_files[1], Lambda2, RankT2, parity2);
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
  
  
  std::vector< double>    OpScalar0b(settings.scalar_op_files.size());
  std::vector< arma::mat> OpScalar1b(settings.scalar_op_files.size());
  std::vector< arma::mat> OpScalar2b(settings.scalar_op_files.size());
  


  // construct a (nJi x nJf) matrix of many-body matrix elements for each Ji, Jf
  arma::field<arma::mat> TensorME1(settings.J2_f.size(),settings.J2_i.size() );
  arma::field<arma::mat> TensorME2(settings.J2_f.size(),settings.J2_i.size() );

  // Below, the vector corresponds to the different scalar operators,
  // the arma::field corresponds to the initial and final J 
  // and the arma::mat corresponds to which states for a given Ji, Jf
  std::vector< arma::field<arma::mat> > ScalarME1(OpScalar1b.size());  
  std::vector< arma::field<arma::mat> > ScalarME2(OpScalar1b.size());

  for (size_t i=0; i<settings.scalar_op_files.size();++i)
  {
//      trans.GetScalarTransitionOperator(settings.scalar_op_files[i],OpScalar0b[i],OpScalar1b[i],OpScalar2b[i]);
      readwrite.ReadScalarTransitionOperator(settings.scalar_op_files[i], trans, OpScalar0b[i],OpScalar1b[i],OpScalar2b[i]);
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

 
// Terrible nomenclature. Here Ji is an index for the array of J values of the initial state.
// This really should be changed.
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
       for (size_t isc=0; isc<OpScalar1b.size();++isc)
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
//             logfile << "2Jf nJf  2Ji nJi  2Lambda = " << std::setw(3) << std::setprecision(1) << settings.J2_f[Jf] << " " << fvec+1
//             << "    " << std::setw(3) << std::setprecision(1) << settings.J2_i[Ji] << " " << ivec+1
//             << "    " << std::setw(3) << std::setprecision(1) << Lambda*2  << std::endl;
//             logfile << "======= One Body Terms ======" << std::endl;
//             logfile << std::setw(4) << "a" << " " << std::setw(4) << "b"
//                     << std::setw(14) << "obtd(a,b)" << " "
//                     << std::setw(14) << "<a||Op||b>" << " "
//                     << std::setw(14) << "obtd * Op" << " " 
//                     << std::setw(14) << "Sum obtd * Op" << std::endl;
//              int nc = TensorOp1b.n_cols;
//              int nr = TensorOp1b.n_rows;
//              float runningsum = 0;
//              for (int c=0;c<nc;c++)
//              {
//                for (int r=0;r<nr;r++)
//                {
//                  runningsum += obtd(r,c) * TensorOp1b(r,c);
//                  logfile << std::setw(4) << std::fixed << c << " " << std::setw(4) << std::fixed << r
//                          << std::setw(14) << std::fixed << std::setprecision(6) << obtd(r,c) << " "
//                          << std::setw(14) << std::fixed << std::setprecision(6) << TensorOp1b(r,c) << " "
//                          << std::setw(14) << std::fixed << std::setprecision(6) << obtd(r,c) * TensorOp1b(r,c) << " " 
//                          << std::setw(14) << std::fixed << std::setprecision(6) << runningsum << std::endl;
//                }
//              }
           }

       }
       
       if (settings.tensor_op_files.size()>1)
       {
         double tbme = arma::accu( TensorOp2b % tbtd);
         TensorME2(indexJf,indexJi)(fvec,ivec) = tbme;

           if (settings.write_log)
           {
             readwrite.WriteLog_Tensor2b( Lambda, TensorOp2b, tbtd );
//             logfile << "======= Two Body Terms ======" << std::endl;
//             logfile << std::setw(4) << "a" << " " << std::setw(4) << "b"
//                     << std::setw(14) << "tbtd(a,b)" << " "
//                     << std::setw(14) << "<a||Op||b>" << " "
//                     << std::setw(14) << "tbtd * Op" << " " 
//                     << std::setw(14) << "Sum tbtd * Op" << std::endl;
//              int nc = TensorOp2b.n_cols;
//              int nr = TensorOp2b.n_rows;
//              float runningsum = 0;
//              for (int c=0;c<nc;c++)
//              {
//                for (int r=0;r<nr;r++)
//                {
//                  runningsum += tbtd(r,c) * TensorOp2b(r,c);
//                  logfile << std::setw(4) << std::fixed << c << " " << std::setw(4) << std::fixed << r
//                          << std::setw(14) << std::fixed << std::setprecision(6) << tbtd(r,c) << " "
//                          << std::setw(14) << std::fixed << std::setprecision(6) << TensorOp2b(r,c) << " "
//                          << std::setw(14) << std::fixed << std::setprecision(6) << tbtd(r,c) * TensorOp2b(r,c) << " " 
//                          << std::setw(14) << std::fixed << std::setprecision(6) << runningsum << std::endl;
//                }
//              }
//              logfile << "#" << std::endl;
//              logfile << "#" << std::endl;
           }

       }
  
       if (settings.J2_i[indexJi]==settings.J2_f[indexJf] and OpScalar1b.size()>0)
       {

         for (size_t iop=0; iop<OpScalar1b.size(); ++iop)
         {
           double scop1 = arma::accu( OpScalar1b[iop] % obtd ) ;
           double scop2 = arma::accu( OpScalar2b[iop] % tbtd ) ;
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


  for (size_t isc=0; isc<settings.scalar_op_files.size(); ++isc)
  {
    oss.str("");
    oss << isc;
    std::string filename = "nutbar_scalar" + oss.str() + "_" + settings.basename_vectors_f.substr( settings.basename_vectors_f.find_last_of("/")+1 );
    if (settings.basename_vectors_i != settings.basename_vectors_f)
    {
       filename +  "_" + settings.basename_vectors_i.substr( settings.basename_vectors_i.find_last_of("/")+1 );
    }
    filename += ".dat";
    std::ofstream scalar_out(filename); 

    scalar_out << "############################################################################################################" << std::endl;
    scalar_out << "# initial vector basename: " << settings.basename_vectors_i << std::endl;
    scalar_out << "# final   vector basename: " << settings.basename_vectors_f << std::endl;
    scalar_out << "# Operator file: " << settings.scalar_op_files[isc] << std::endl;
    scalar_out << "# Zero body term: " << OpScalar0b[isc] << std::endl;
    scalar_out << "# Jf  nJf     Ji  nJi      Ei          Ef         <Op1b>          <Op2b>         <Op1b+2b>     <Op0b+1b+2b>" << std::endl;
    scalar_out << "###########################################################################################################" << std::endl;
  
    for (size_t indexJf=0; indexJf<settings.J2_f.size();++indexJf)
    {
      double Jf = settings.J2_f[indexJf]*0.5;
      for (size_t njf=0; njf<settings.NJ_f[indexJf];++njf)
      {
        for (size_t indexJi=0; indexJi<settings.J2_i.size();++indexJi)
        {
         double Ji = settings.J2_i[indexJi]*0.5;
         if (std::abs(Ji-Jf)>0.1) continue;
         for (size_t nji=0; nji<settings.NJ_i[indexJi];++nji)
         {
           if ( ( std::find( begin(settings.options), end(settings.options), "diag") != end(settings.options)) and not( Ji==Jf and nji==njf))  continue;
           double zerobody = (nji==njf) ? OpScalar0b[isc] : 0;
           scalar_out << std::fixed << std::setw(4) << std::setprecision(1) << Jf << " " << std::setw(3) << njf+1 << "    "
                      << std::fixed << std::setw(4) << std::setprecision(1) << Ji << " " << std::setw(3) << nji+1 << " "
                      << std::fixed << std::setw(10) << std::setprecision(3) << trans.nuvec_list_f[indexJf].alpha[njf] << "  "
                      << std::fixed << std::setw(10) << std::setprecision(3) << trans.nuvec_list_i[indexJi].alpha[nji] << "  "
                      << std::scientific << std::setw(14) << std::setprecision(6) << ScalarME1[isc](indexJf,indexJi)(njf,nji) << "  "
                      << std::scientific << std::setw(14) << std::setprecision(6) << ScalarME2[isc](indexJf,indexJi)(njf,nji) << "  "
                      << std::scientific << std::setw(14) << std::setprecision(6) << ScalarME1[isc](indexJf,indexJi)(njf,nji) + ScalarME2[isc](indexJf,indexJi)(njf,nji) << "  "
                      << std::scientific << std::setw(14) << std::setprecision(6) << ScalarME1[isc](indexJf,indexJi)(njf,nji) + ScalarME2[isc](indexJf,indexJi)(njf,nji) + zerobody 
                      << std::endl;
         }
        }
      }
    }



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
