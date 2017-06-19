//
//  Program nutbar  (NuShellX Transtitions from Binary Arrays)
//
//

#include <iostream>
#include <iomanip>
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


using namespace std;


// Define a struct to hold all the input information
struct Settings
{
  string basename_vectors_i; // base name for the vectors, e.g. ne200
  string basename_vectors_f; // base name for the vectors, e.g. ne200
  string basename_sps;     // base name for the sp and sps files, e.g. sdpn
  vector<string> scalar_op_files; // list of scalar operator file names
  vector<string> tensor_op_files; // list of tensor operator file names
  vector<int> J2_i;  // 2*J for initial states
  vector<int> NJ_i;  // number of states for each initial J
  vector<int> J2_f;  // 2*J for final states
  vector<int> NJ_f;  // number of states for each final J
  vector<string> options;
//  bool write_egv; // option to write out eigenvector in format readable by Petr Navratil's TRDENS code
};

// forward declaration. implementation is below main().
Settings ReadInput(istream& input, string mode);
arma::field<arma::mat> operator+(arma::field<arma::mat>&,arma::field<arma::mat>&);




int main(int argc, char** argv)
{

  // TransitionDensity class does all the heavy lifting
  Profiler profiler;
  TransitionDensity trans;
  
  Settings settings;
  
  if (argc>1)
  {
    ifstream infile(argv[1]);
    settings = ReadInput(infile, "batch");
  }
  else
  {
    settings = ReadInput(cin, "interactive");
  }
  
  cout << "options: [ ";
  for (string opt : settings.options) cout << opt << " ";
  cout << " ]" << endl; 
 

  // Some shuffling around to put things in the TransitionDensity class format
  // this should really be done in a much neater way
  trans.basename_i = settings.basename_vectors_i;
  trans.basename_f = settings.basename_vectors_f;
  trans.sps_file_name = settings.basename_sps + ".sps";
  trans.Jlist_i = settings.J2_i;
  trans.Jlist_f = settings.J2_f;
  for (size_t i=0;i<settings.J2_i.size();++i)
  {
    trans.max_states_per_J_i[settings.J2_i[i]] = settings.NJ_i[i];
  }
  for (size_t i=0;i<settings.J2_f.size();++i)
  {
    trans.max_states_per_J_f[settings.J2_f[i]] = settings.NJ_f[i];
  }


#ifdef VERBOSE 
    cout << "Reading Files" << endl;
#endif
  trans.ReadFiles();
#ifdef VERBOSE
    cout << "done reading files." << endl;
    cout << "CalculateMschemeAmplitudes()" << endl;
#endif
  
  if ( find( begin(settings.options), end(settings.options), "read_dens") == end(settings.options)  )
    trans.CalculateMschemeAmplitudes();

#ifdef VERBOSE
    cout << "done calculating Mscheme amplitudes." << endl;
    cout << "number of initial m-scheme states: " << trans.amplitudes_i.size() << " x " << trans.blank_vector_i.size()<< endl;
    cout << "number of final m-scheme states: " << trans.amplitudes_f.size() << " x " << trans.blank_vector_f.size()<< endl;
#endif
  
//  if (settings.write_egv)
  if ( find( begin(settings.options), end(settings.options), "egv")  != end(settings.options))
  {
    trans.WriteEGV("mbpt.egv");
    trans.WriteTRDENS_input("trdens.in");
  }
  
  cout << "operators: [ " ;
  for (auto ops : settings.scalar_op_files) cout << ops << " ";
  for (auto ops : settings.tensor_op_files) cout << ops << " ";
  cout << "] " << endl;

  arma::mat TensorOp1b;
  arma::mat TensorOp2b;


  int Lambda=0,RankT=0,parity=0;
  if (settings.tensor_op_files.size() > 0)
  {
    TensorOp1b = trans.GetOneBodyTransitionOperator(settings.tensor_op_files[0], Lambda, RankT, parity);
//    cout << "Op1b " << endl << TensorOp1b << endl;
  }
  if (settings.tensor_op_files.size() > 1)
  {
    int Lambda2,RankT2,parity2;
    TensorOp2b = trans.GetTwoBodyTransitionOperator(settings.tensor_op_files[1], Lambda2, RankT2, parity2);
    if (Lambda2 != Lambda or RankT2!=RankT or parity2 != parity)
    {
      cout << "ERROR Tensor rank mismatch between files " << settings.tensor_op_files[0] << "  and  " << settings.tensor_op_files[1] << endl;
      return 1;
    }
  }
  
  
  vector<double> OpScalar0b(settings.scalar_op_files.size());
  vector<arma::mat> OpScalar1b(settings.scalar_op_files.size());
  vector<arma::mat> OpScalar2b(settings.scalar_op_files.size());
  


  // construct a (nJi x nJf) matrix of many-body matrix elements for each Ji, Jf
  arma::field<arma::mat> TensorME1(settings.J2_f.size(),settings.J2_i.size() );
  arma::field<arma::mat> TensorME2(settings.J2_f.size(),settings.J2_i.size() );
  vector<arma::field<arma::mat> > ScalarME1(OpScalar1b.size());
  vector<arma::field<arma::mat> > ScalarME2(OpScalar1b.size());



  for (size_t i=0; i<settings.scalar_op_files.size();++i)
  {
      trans.GetScalarTransitionOperator(settings.scalar_op_files[i],OpScalar0b[i],OpScalar1b[i],OpScalar2b[i]);
      ScalarME1[i] = arma::field<arma::mat>(settings.J2_f.size(),settings.J2_i.size() );
      ScalarME2[i] = arma::field<arma::mat>(settings.J2_f.size(),settings.J2_i.size() );
  }
 


  trans.SetupKets();

  if ( find( begin(settings.options), end(settings.options), "read_dens") == end(settings.options)  )
    trans.SetDensFile("nutbar_densities.dat");
 
  for (size_t Ji=0;Ji<settings.J2_i.size();++Ji )
  {
   for (size_t Jf=0;Jf<settings.J2_f.size();++Jf )
   {
     TensorME1(Jf,Ji).zeros(settings.NJ_f[Jf],settings.NJ_i[Ji]);
     TensorME2(Jf,Ji).zeros(settings.NJ_f[Jf],settings.NJ_i[Ji]);
     if ( abs(settings.J2_i[Ji] - settings.J2_f[Jf])> 2*Lambda) continue;
     if (Jf==Ji)
     {
       for (size_t isc=0; isc<OpScalar1b.size();++isc)
       {
          ScalarME1[isc](Jf,Ji).zeros(settings.NJ_f[Jf],settings.NJ_i[Ji]);
          ScalarME2[isc](Jf,Ji).zeros(settings.NJ_f[Jf],settings.NJ_i[Ji]);
       }
     }
     for (int ivec = 0; ivec<settings.NJ_i[Ji]; ++ivec)
     {
      for (int fvec = 0; fvec<settings.NJ_f[Jf]; ++fvec)
      {
       if (settings.basename_vectors_i == settings.basename_vectors_f)
       {
        if ( settings.J2_i[Ji]==settings.J2_f[Jf] and fvec>ivec) continue;
        if ( ( find( begin(settings.options), end(settings.options), "diag") != end(settings.options)) and not( settings.J2_i[Ji]==settings.J2_f[Jf] and ivec==fvec))  continue;
       }
       arma::mat obtd,tbtd;
       if ( find( begin(settings.options), end(settings.options), "read_dens") != end(settings.options) )
       {
         cout << "Reading densities from file..." << endl;
         obtd = trans.ReadOBTD(Ji,ivec,Jf,fvec,Lambda*2,"nutbar_densities.dat");
         tbtd = trans.ReadTBTD(Ji,ivec,Jf,fvec,Lambda*2,"nutbar_densities.dat");
       }
       else
       {
         obtd = trans.CalcOBTD(Ji,ivec,Jf,fvec,Lambda*2);
         if (settings.tensor_op_files.size() > 1 or settings.scalar_op_files.size()>0)
           tbtd = trans.CalcTBTD(Ji,ivec,Jf,fvec,Lambda*2);
       }
 
       if (settings.tensor_op_files.size()>0)
       {
         double obme = arma::accu( TensorOp1b % obtd );
         TensorME1(Jf,Ji)(fvec,ivec) = obme;
       }
       
       if (settings.tensor_op_files.size()>1)
       {
         double tbme = arma::accu( TensorOp2b % tbtd);
         TensorME2(Jf,Ji)(fvec,ivec) = tbme;
       }
  
       if (settings.J2_i[Ji]==settings.J2_f[Jf] and OpScalar1b.size()>0)
       {

         for (size_t iop=0; iop<OpScalar1b.size(); ++iop)
         {
           double scop1 = arma::accu( OpScalar1b[iop] % obtd ) ;
           double scop2 = arma::accu( OpScalar2b[iop] % tbtd ) ;
           ScalarME1[iop](Jf,Ji)(fvec,ivec) = scop1;
           ScalarME2[iop](Jf,Ji)(fvec,ivec) = scop2;
         }
       }
      }
    }
   }
  }

  ostringstream oss;
  if ( settings.tensor_op_files.size()>0)
  {
    oss << Lambda;
    string filename = "nutbar_tensor" + oss.str() + "_" + settings.basename_vectors_f.substr( settings.basename_vectors_f.find_last_of("/")+1 );
    if (settings.basename_vectors_i != settings.basename_vectors_f)
    {
       filename +  "_" + settings.basename_vectors_i.substr( settings.basename_vectors_i.find_last_of("/")+1 );
    }
    filename += ".dat";
    ofstream tensor_out(filename);
  
    tensor_out << "##########################################################################################" << endl;
    tensor_out << "# initial vector basename: " << settings.basename_vectors_i << endl;
    tensor_out << "# final   vector basename: " << settings.basename_vectors_f << endl;
    tensor_out << "# One body file: " << settings.tensor_op_files[0] << endl;
    if (settings.tensor_op_files.size()>1)
      tensor_out << "# Two body file: " << settings.tensor_op_files[1] << endl;
    tensor_out << "# Jf  nJf     Ji  nJi      Ef          Ei         <Op1b>          <Op2b>         <Op1b+2b> " << endl;
    tensor_out << "###########################################################################################" << endl;
  
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
           tensor_out << fixed << setw(4) << setprecision(1) << Jf << " " << setw(3) << njf+1 << "    "
                      << fixed << setw(4) << setprecision(1) << Ji << " " << setw(3) << nji+1 << " "
                      << fixed << setw(10) << setprecision(3) << trans.nuvec_list_f[indexJf].alpha[njf] << "  "
                      << fixed << setw(10) << setprecision(3) << trans.nuvec_list_i[indexJi].alpha[nji] << "  "
                      << scientific << setw(14) << setprecision(6) << TensorME1(indexJf,indexJi)(njf,nji) << "  "
                      << scientific << setw(14) << setprecision(6) << TensorME2(indexJf,indexJi)(njf,nji) << "  "
                      << scientific << setw(14) << setprecision(6) << TensorME1(indexJf,indexJi)(njf,nji) + TensorME2(indexJf,indexJi)(njf,nji) << "  "
                      << endl;
         }
        }
      }
    }
  }


  for (size_t isc=0; isc<settings.scalar_op_files.size(); ++isc)
  {
    oss.str("");
    oss << isc;
    string filename = "nutbar_scalar" + oss.str() + "_" + settings.basename_vectors_f.substr( settings.basename_vectors_f.find_last_of("/")+1 );
    if (settings.basename_vectors_i != settings.basename_vectors_f)
    {
       filename +  "_" + settings.basename_vectors_i.substr( settings.basename_vectors_i.find_last_of("/")+1 );
    }
    filename += ".dat";
    ofstream scalar_out(filename); 

    scalar_out << "############################################################################################################" << endl;
    scalar_out << "# initial vector basename: " << settings.basename_vectors_i << endl;
    scalar_out << "# final   vector basename: " << settings.basename_vectors_f << endl;
    scalar_out << "# Operator file: " << settings.scalar_op_files[isc] << endl;
    scalar_out << "# Zero body term: " << OpScalar0b[isc] << endl;
    scalar_out << "# Jf  nJf     Ji  nJi      Ei          Ef         <Op1b>          <Op2b>         <Op1b+2b>     <Op0b+1b+2b>" << endl;
    scalar_out << "###########################################################################################################" << endl;
  
    for (size_t indexJf=0; indexJf<settings.J2_f.size();++indexJf)
    {
      double Jf = settings.J2_f[indexJf]*0.5;
      for (size_t njf=0; njf<settings.NJ_f[indexJf];++njf)
      {
        for (size_t indexJi=0; indexJi<settings.J2_i.size();++indexJi)
        {
         double Ji = settings.J2_i[indexJi]*0.5;
         if (abs(Ji-Jf)>0.1) continue;
         for (size_t nji=0; nji<settings.NJ_i[indexJi];++nji)
         {
           if ( ( find( begin(settings.options), end(settings.options), "diag") != end(settings.options)) and not( Ji==Jf and nji==njf))  continue;
           double zerobody = (nji==njf) ? OpScalar0b[isc] : 0;
           scalar_out << fixed << setw(4) << setprecision(1) << Jf << " " << setw(3) << njf+1 << "    "
                      << fixed << setw(4) << setprecision(1) << Ji << " " << setw(3) << nji+1 << " "
                      << fixed << setw(10) << setprecision(3) << trans.nuvec_list_f[indexJf].alpha[njf] << "  "
                      << fixed << setw(10) << setprecision(3) << trans.nuvec_list_i[indexJi].alpha[nji] << "  "
                      << scientific << setw(14) << setprecision(6) << ScalarME1[isc](indexJf,indexJi)(njf,nji) << "  "
                      << scientific << setw(14) << setprecision(6) << ScalarME2[isc](indexJf,indexJi)(njf,nji) << "  "
                      << scientific << setw(14) << setprecision(6) << ScalarME1[isc](indexJf,indexJi)(njf,nji) + ScalarME2[isc](indexJf,indexJi)(njf,nji) << "  "
                      << scientific << setw(14) << setprecision(6) << ScalarME1[isc](indexJf,indexJi)(njf,nji) + ScalarME2[isc](indexJf,indexJi)(njf,nji) + zerobody 
                      << endl;
         }
        }
      }
    }



  }



  profiler.PrintAll();

  return 0;
}




Settings ReadInput(istream& input, string mode)
{

  Settings settings;
  string line;
  istringstream iss;

  if (mode == "interactive")
  {
    settings.basename_sps = "(none found)"; // set default
#ifndef NOBOOST
    for (auto& path : boost::filesystem::directory_iterator("."))
    {
      string pathstring = path.path().string();
      if (pathstring.substr( pathstring.find_last_of(".") ) == ".sps")
      {
        if (pathstring.find("bab") == string::npos)
        {
           settings.basename_sps = pathstring.substr( 0,pathstring.find_last_of(".") );
           break;
        }
      }
    }
#endif
    cout << "basename for *.sps files, default [" << settings.basename_sps << "]: " << flush;
  }

  getline(input,line);
  if (mode !="interactive" or line !="")
  {
    iss.clear();
    iss.str(line);
    iss >> settings.basename_sps;
  }

  if (mode == "interactive")
  {
    settings.basename_vectors_i = "(none found)"; // set default
#ifndef NOBOOST
    for (auto& path : boost::filesystem::directory_iterator("."))
    {
      string pathstring = path.path().string();
      if (pathstring.substr( pathstring.find_last_of(".") ) == ".xvc")
      {
        settings.basename_vectors_i = pathstring.substr( 0, pathstring.find_last_of("/")+6 );
        break;
      }
    }
#endif
    cout << "basename for initial state *.xvc files, default [" << settings.basename_vectors_i << "]: " << flush;
  }

  getline(input,line);
  if (mode!="interactive" or line !="")
  {
     iss.clear();
     iss.str(line);
     iss >> settings.basename_vectors_i;
  }

  if (mode == "interactive")
  {
    settings.basename_vectors_f = settings.basename_vectors_i; // set default
    cout << "basename for final state *.xvc files, default [" << settings.basename_vectors_f << "]: " << flush;
  }

  getline(input,line);
  if (mode!="interactive" or line !="")
  {
    iss.clear();
    iss.str(line);
    iss >> settings.basename_vectors_f;
  }




  if (mode == "interactive")
  {
    cout << "list operator files (scalars end in .int, tensors end in .op, separate with spaces): " << flush;
  }

  getline(input,line);
  iss.clear();
  iss.str(line);
  string opfilename;
  while (iss >> opfilename)
  {
    if (opfilename.find_first_of("#!") != string::npos) break;
    if (opfilename.substr(opfilename.find_last_of(".")) == ".int")
    {
      settings.scalar_op_files.push_back(opfilename);
    }
    else if (opfilename.substr(opfilename.find_last_of(".")) == ".op")
    {
      settings.tensor_op_files.push_back(opfilename);
    }
    else
    {
      cout << "ERROR: Unrecognized file extension for file " << opfilename << endl;
      return settings;
    }
  }

  if (mode == "interactive")
  {
    cout << "list J values for initial states: " << flush;
  }
  getline(input,line);
  iss.clear();
  iss.str(line);
  double J;
  while ( iss>>J )
  {
    settings.J2_i.push_back(int(J*2));
  }

  if (mode == "interactive")
  {
    cout << "list number of states for each initial J: " << flush;
  }
  getline(input,line);
  iss.clear();
  iss.str(line);
  double nJ;
  while ( iss>>nJ )
  {
    settings.NJ_i.push_back(int(nJ));
  }
  // if we didn't read enough values, just repeat the last entered value until we have enough
  if (settings.NJ_i.size()<1) settings.NJ_i.push_back(1);
  while (settings.J2_i.size()>settings.NJ_i.size()) settings.NJ_i.push_back(settings.NJ_i.back());

  
  if (mode == "interactive")
  {
    cout << "list J values for final states: " << flush;
  }
  getline(input,line);
  iss.clear();
  iss.str(line);
  while ( iss>>J )
  {
    settings.J2_f.push_back(int(J*2));
  }

  if (mode == "interactive")
  {
    cout << "list number of states for each final J: " << flush;
  }
  getline(input,line);
  iss.clear();
  iss.str(line);
  while ( iss>>nJ )
  {
    settings.NJ_f.push_back(int(nJ));
  }
  if (settings.NJ_f.size()<1) settings.NJ_f.push_back(1);
  while (settings.J2_f.size()>settings.NJ_f.size()) settings.NJ_f.push_back(settings.NJ_f.back());

  if (mode == "interactive")
  {
//    cout << "output eigenvectors as mbpt.egv? (Y/[n]): " << flush;
    cout << "list any other option keywords (egv diag): " << flush;
  }
  getline(input,line);
  iss.clear();
  iss.str(line);
  string s;
  while (iss >> s) settings.options.push_back(s);

  string outfilename = "nutbar_" + settings.basename_vectors_f.substr( settings.basename_vectors_f.find_last_of("/")+1 );
  if (settings.basename_vectors_i != settings.basename_vectors_f)
  {
     outfilename +  "_" + settings.basename_vectors_i.substr( settings.basename_vectors_i.find_last_of("/")+1 );
  }
  outfilename += ".input";
  ofstream outfile(outfilename);
  outfile << settings.basename_sps << endl;
  outfile << settings.basename_vectors_i << endl;
  outfile << settings.basename_vectors_f << endl;
  for (auto op : settings.scalar_op_files ) outfile << op << " ";
  for (auto op : settings.tensor_op_files ) outfile << op << " ";
  outfile << endl;
  for (auto J : settings.J2_i)  outfile << fixed << setprecision(1) << float(J)/2 << " ";
  outfile << endl;
  for (auto nJ : settings.NJ_i)  outfile << nJ << " ";
  outfile << endl;
  for (auto J : settings.J2_f)  outfile << fixed << setprecision(1) << float(J)/2 << " ";
  outfile << endl;
  for (auto nJ : settings.NJ_f)  outfile << nJ << " ";
  outfile << endl;
  for (auto opt : settings.options) outfile << opt << " ";
  outfile << endl;
//  if (settings.write_egv) outfile << "y" << endl;
//  else outfile << "n" << endl;

  return settings;

}


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
