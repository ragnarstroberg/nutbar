

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#ifndef NOBOOST
#include <boost/filesystem.hpp>
#endif


#include "ReadWrite.hh"
#include "TransitionDensity.hh"
#include "NuVec.hh"
#include "Operators.hh"
#include "Settings.hh"

#define SQRT2 1.4142135623730950488


// a and b are presumably skipped here because a and b are used to label proton/neutron
std::vector<char> ReadWrite::an_code =    {  '0','1','2','3','4','5','6','7','8','9',
                                             '_','-','c','d','e','f','g','h','i','j',  
                                             'k','l','m','n','o','p','q','r','s','t',  
                                             'u','v','w','x','y','z','A','B','C','D',  
                                             'E','F','G','H','I','J','K','L','M','N',  
                                             'O','P','Q','R','S','T','U','V','W','X',  
                                             '0','1','2','3','4','5','6','7','8','9',  
                                             '_','-','c'};

std::vector<std::string> ReadWrite::periodic_table = {
  "xx","h_","he",
       "li","be","b_","c_","n_","o_","f_","ne",
       "na","mg","al","si","p_","s_","cl","ar",
       "k_","ca","sc","ti","v_","cr","mn","fe","co","ni","cu","zn","ga","ge","as","se","br","kr",
       "rb","sr","y_","zr","nb","mo","tc","ru","rh","pd","ag","cd","in","sn","sb","te","i_","xe",
       "cs","ba",
         "la","ce","pr","nd","pm","sm","eu","gd","tb","dy","ho","er","tm","yb",
                 "lu","hf","ta","w_","re","os","ir","pr","au","hg","tl","pb","bi","po","at","rn",
       "fr","ra",
         "ac","th","u_","np","pu","am","cm","bk","cf","es","fm","md","no","lr",
                 "rf","db","sg","bh","hs","mt","ds","rg","cn"};










void ReadWrite::GetAZFromFileName(  )
{


  GetCoreFromSPfile(); // make sure we know what the core is.

  std::string basename_i = settings.basename_vectors_i;
  std::string basename_f = settings.basename_vectors_f;
  std::string trimmed_basename_i = (basename_i.find("/")==std::string::npos) ? basename_i : basename_i.substr( basename_i.find_last_of("/")+1 );
  std::string trimmed_basename_f = (basename_f.find("/")==std::string::npos) ? basename_f : basename_f.substr( basename_f.find_last_of("/")+1 );
  std::string element_i = trimmed_basename_i.substr( 0,2);
  std::string element_f = trimmed_basename_f.substr( 0,2);

  auto el_position_i = std::find( periodic_table.begin(),periodic_table.end(), element_i);
  auto el_position_f = std::find( periodic_table.begin(),periodic_table.end(), element_f);
  if (el_position_i == periodic_table.end())
  {
   std::cout << "ERROR! : could not find " << element_i << " in periodic table" << std::endl;
   return;
  }
  if (el_position_f == periodic_table.end())
  {
   std::cout << "ERROR! : could not find " << element_i << " in periodic table" << std::endl;
   return;
  }
  settings.Z_i = el_position_i - periodic_table.begin();
  settings.Z_f = el_position_f - periodic_table.begin();
  std::istringstream( trimmed_basename_i.substr(2,2) ) >> settings.A_i;
  std::istringstream( trimmed_basename_f.substr(2,2) ) >> settings.A_f;


  if ( settings.A_i < settings.Acore ) settings.A_i +=100;
  if ( settings.A_f < settings.Acore ) settings.A_f +=100;
  
}





// Read the A and Z of the core from the *.sp file
// Is this really all we need from that?
void ReadWrite::GetCoreFromSPfile()
{
  std::string sp_file_name = settings.basename_sps + ".sp"; // make sure it's .sp and not .sps
  std::ifstream infile(sp_file_name);
  if (not infile.good())
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "TransitionDensity::ReadSPfile -- error reading " << sp_file_name << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    exit(EXIT_FAILURE);
  }
  infile.ignore(256,'\n').ignore(256,'\n'); // skip two lines
  infile >> settings.Acore >> settings.Zcore;
}








void ReadWrite::ReadInput(std::istream& input, std::string mode)
{

  std::string line;
  std::istringstream iss;

  if (mode == "interactive")
  {
    settings.basename_sps = "(none found)"; // set default
#ifndef NOBOOST
    for (auto& path : boost::filesystem::directory_iterator("."))
    {
      std::string pathstring = path.path().string();
      if (pathstring.substr( pathstring.find_last_of(".") ) == ".sps")
      {
        if (pathstring.find("bab") == std::string::npos)
        {
           settings.basename_sps = pathstring.substr( 0,pathstring.find_last_of(".") );
           break;
        }
      }
    }
#endif
    std::cout << "basename for *.sps files, default [" << settings.basename_sps << "]: " << std::flush;
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
      std::string pathstring = path.path().string();
      if (pathstring.substr( pathstring.find_last_of(".") ) == ".xvc")
      {
        settings.basename_vectors_i = pathstring.substr( 0, pathstring.find_last_of("/")+6 );
        break;
      }
    }
#endif
    std::cout << "basename for initial state *.xvc files, default [" << settings.basename_vectors_i << "]: " << std::flush;
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
    std::cout << "basename for final state *.xvc files, default [" << settings.basename_vectors_f << "]: " << std::flush;
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
    std::cout << "list operator files (scalars end in .int, tensors end in .op, SF's end in .dag, separate with spaces): " << std::flush;
  }

  getline(input,line);
  iss.clear();
  iss.str(line);
  std::string opfilename;
  while (iss >> opfilename)
  {
    if (opfilename.find_first_of("#!") != std::string::npos) break;
    if (opfilename.substr(opfilename.find_last_of(".")) == ".int")
    {
      settings.scalar_op_files.push_back(opfilename);
    }
    else if (opfilename.substr(opfilename.find_last_of(".")) == ".op")
    {
      settings.tensor_op_files.push_back(opfilename);
    }
    else if (opfilename.substr(opfilename.find_last_of(".")) == ".dag")
    {
      std::cout << "pushing back file name " << opfilename << std::endl;
      settings.dagger_op_files.push_back(opfilename);
    }
    else
    {
      std::cout << "ERROR: Unrecognized file extension for file " << opfilename << std::endl;
      return ;
    }
  }

  if (mode == "interactive")
  {
    std::cout << "list J values for initial states: " << std::flush;
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
    std::cout << "list number of states for each initial J: " << std::flush;
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
    std::cout << "list J values for final states: " << std::flush;
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
    std::cout << "list number of states for each final J: " << std::flush;
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
    std::cout << "list any other option keywords (egv, diag, read_dens, log): " << std::flush;
  }
  getline(input,line);
  iss.clear();
  iss.str(line);
  std::string s;
  while (iss >> s) settings.options.push_back(s);

  std::string outfilename = "nutbar_" + settings.basename_vectors_f.substr( settings.basename_vectors_f.find_last_of("/")+1 );
  if (settings.basename_vectors_i != settings.basename_vectors_f)
  {
     outfilename +  "_" + settings.basename_vectors_i.substr( settings.basename_vectors_i.find_last_of("/")+1 );
  }
  outfilename += ".input";
  std::ofstream outfile(outfilename);
  outfile << settings.basename_sps << std::endl;
  outfile << settings.basename_vectors_i << std::endl;
  outfile << settings.basename_vectors_f << std::endl;
  for (auto op : settings.scalar_op_files ) outfile << op << " ";
  for (auto op : settings.tensor_op_files ) outfile << op << " ";
  for (auto op : settings.dagger_op_files ) outfile << op << " ";
  outfile << std::endl;
  for (auto J : settings.J2_i)  outfile << std::fixed << std::setprecision(1) << float(J)/2 << " ";
  outfile << std::endl;
  for (auto nJ : settings.NJ_i)  outfile << nJ << " ";
  outfile << std::endl;
  for (auto J : settings.J2_f)  outfile << std::fixed << std::setprecision(1) << float(J)/2 << " ";
  outfile << std::endl;
  for (auto nJ : settings.NJ_f)  outfile << nJ << " ";
  outfile << std::endl;
  for (auto opt : settings.options) outfile << opt << " ";
  outfile << std::endl;


  settings.densities_from_file = ( std::find( std::begin(settings.options), std::end(settings.options), "read_dens")  != std::end(settings.options) );
  settings.write_egv =           ( std::find( std::begin(settings.options), std::end(settings.options), "egv")        != std::end(settings.options) );
  settings.write_log =           ( std::find( std::begin(settings.options), std::end(settings.options), "log")        != std::end(settings.options) );
  settings.diagonal_only =       ( std::find( std::begin(settings.options), std::end(settings.options), "diag")       != std::end(settings.options) );
  settings.same_basename_fi = (settings.basename_vectors_i == settings.basename_vectors_f);

  settings.same_states_fi = ( settings.same_basename_fi and (settings.J2_i.size() == settings.J2_f.size()) );
  if ( settings.same_states_fi )
  {
    for (size_t i=0;i<settings.J2_i.size();i++)
    {
      if ( (settings.J2_i[i]!=settings.J2_f[i])  or (settings.NJ_i[i]!=settings.NJ_f[i]) )
      {
        settings.same_states_fi = false;
        break;
      }
    }
  }


  std::string densityfilename = "nutbar_densities_" + settings.basename_vectors_f.substr( settings.basename_vectors_f.find_last_of("/")+1 ) + ".dat";

  if (settings.densities_from_file)
  {
    densityfile.open( densityfilename, std::ios_base::in );
  }
  else
  {
    densityfile.open( densityfilename, std::ios_base::out | std::ios_base::trunc ); // need the trunc flag to create a file if it doesn't exist
  }

  if (settings.write_log)
  {
    std::string logfilename = "nutbar_" + settings.basename_vectors_f.substr( settings.basename_vectors_f.find_last_of("/")+1 ) + ".log";
    logfile.open(logfilename, std::ios_base::out);
    std::cout << "Writing more detailed info to " << logfilename << std::endl;
  }


}








arma::mat ReadWrite::ReadOBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, TransitionDensity& trans)
{
  auto& m_orbits = settings.m_orbits;
  auto& jorbits = settings.jorbits;

  bool found_it = false;
  int Ji = settings.J2_i[J_index_i];
  int Jf = settings.J2_f[J_index_f];
  size_t njorb = jorbits.size();
  densityfile.seekg(std::ios_base::beg); // go back to the beginning of the file (we left it open)
  if (not densityfile.good())
  {
   std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
   std::cout << "!! TransitionDensity::ReadOBTD -- trouble reading density file" << std::endl;
   std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
   exit(EXIT_FAILURE);
  }
  std::string line;
  std::ostringstream line_to_find;
  line_to_find <<  "Jf nJf  Ji nJi  Lambda = "   << std::setw(3) << std::setprecision(1) << Jf*0.5 << " " << eigvec_f+1
          << "    " << std::setw(3) << std::setprecision(1) << Ji*0.5 << " " << eigvec_i+1
          << "    " << std::setw(3) << std::setprecision(1) << Lambda2*0.5;


  while ( getline( densityfile, line) )
  {
     if( line == line_to_find.str() )
     {
       found_it = true;
       break;
     }
  }
  while ( getline( densityfile, line) )
  {
     if( line == "-------------- OBTD ---------------------" ) break;
  }

  bool same_i_f =  ( (settings.basename_vectors_i==settings.basename_vectors_f) and (Ji==Jf) and (eigvec_i==eigvec_f) );
  arma::mat obtd(njorb,njorb,arma::fill::zeros);
  int i,j;
  double obd;
  while ( line.size() > 2 )
  {
    getline( densityfile, line );
    std::istringstream(line) >> i >> j >> obd;
    obtd(i,j) = obd;
    if (same_i_f)
    {
      int j2i = m_orbits[jorbits[i]].j2;
      int j2j = m_orbits[jorbits[j]].j2;
      obtd(j,i) = (1-std::abs(j2j-j2i)%4) * obtd(i,j);
    }

  }
//  densityfile.close();
  if (not found_it)
  {
    std::cout << "WARNING!! Didn't find " << line_to_find.str() << "  in nutbar_densities.dat " << std::endl;
  }
  return obtd;
}


arma::mat ReadWrite::ReadTBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, TransitionDensity& trans)
{

  auto& ket_a = settings.ket_a;
  auto& ket_b = settings.ket_b;
  auto& ket_J = settings.ket_J;

  bool found_it = false;
  int Ji = settings.J2_i[J_index_i];
  int Jf = settings.J2_f[J_index_f];
  densityfile.seekg(std::ios_base::beg); // go back to the beginning of the file (we left it open)
  if (not densityfile.good())
  {
   std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
   std::cout << "!! TransitionDensity::ReadTBTD -- trouble reading density file"  << std::endl;
   std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
   exit(EXIT_FAILURE);
  }
  std::string line;
  std::ostringstream line_to_find;
  line_to_find << "Jf nJf  Ji nJi  Lambda = "  << std::setw(3) << std::setprecision(1) << Jf*0.5 << " " << eigvec_f+1
          << "    " << std::setw(3) << std::setprecision(1) << Ji*0.5 << " " << eigvec_i+1
          << "    " << std::setw(3) << std::setprecision(1) << Lambda2*0.5;


  while ( getline( densityfile, line) )
  {
     if( line == line_to_find.str() )
     {
       found_it = true;
       break;
     }
  }
  while ( getline( densityfile, line) )
  {
     if( line == "-------------- TBTD ---------------------" ) break;
  }
  arma::mat tbtd(ket_J.size(), ket_J.size(), arma::fill::zeros);
  bool same_i_f =  ( (settings.basename_vectors_i==settings.basename_vectors_f) and (Ji==Jf) and (eigvec_i==eigvec_f) );
  int ibra,iket;
  double tbd;
  while ( line.size() > 2 )
  {
    getline( densityfile, line );
    std::istringstream(line) >> ibra >> iket >> tbd;

    int J2ab = ket_J[ibra];
    int J2cd = ket_J[iket];
    tbtd(ibra,iket) = tbd;

    if (same_i_f)
    {
      tbtd(iket,ibra) = tbtd(ibra,iket) * (1-std::abs(J2ab-J2cd)%4); 
    }
  }
  if (not found_it)
  {
    std::cout << "WARNING!! Didn't find " << line_to_find.str() << "  in nutbar_densities.dat " << std::endl;
  }
  return tbtd;
}






void ReadWrite::WriteDensityHeader(  )
{
  auto& m_orbits = settings.m_orbits;
  auto& jorbits = settings.jorbits;
  auto& ket_a = settings.ket_a;
  auto& ket_b = settings.ket_b;
  auto& ket_J = settings.ket_J;
  densityfile << "# One and two body transition densities" << std::endl << "#" << std::endl;
  densityfile << "# One-body basis:" << std::endl;
  densityfile << "# i   n   l   2j  2tz (proton=+1)" << std::endl;
  for (size_t i=0;i<jorbits.size();++i)
  {
    auto morb = m_orbits[jorbits[i]];
    densityfile << std::setw(3) << i << " " << std::setw(3) << morb.n << " " << std::setw(3) << morb.l2/2 << " "
            << std::setw(3) << morb.j2 << " " << std::setw(3) << morb.tz2 << std::endl;
  }

  densityfile << "# Two-body basis: " << std::endl;
  densityfile << "# i   a   b   J" << std::endl;
  for (size_t i=0;i<ket_a.size();++i)
  {
    densityfile << std::setw(3) << i << " " << std::setw(3) << ket_a[i] << " " << std::setw(3) << ket_b[i] << " " << std::setw(3) << ket_J[i]/2 << std::endl;
  }

}





void ReadWrite::WriteOBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, arma::mat& obtd)
{

  int Ji = settings.J2_i[J_index_i];
  int Jf = settings.J2_f[J_index_f];
  size_t njorb = obtd.n_rows;

  if ( densityfile.fail() )
  {
    std::cerr << __func__ << "  ERROR:  couldn't write to densityfile." << std::endl;
    return;
  }
  

  densityfile << std::endl;
  densityfile << "Jf nJf  Ji nJi  Lambda = " << std::setw(3) << std::setprecision(1) << Jf*0.5 << " " << eigvec_f+1
          << "    " << std::setw(3) << std::setprecision(1) << Ji*0.5 << " " << eigvec_i+1
          << "    " << std::setw(3) << std::setprecision(1) << Lambda2*0.5  << std::endl;
  densityfile << "-------------- OBTD ---------------------" << std::endl;
 

  bool same_i_f =  ( (settings.same_basename_fi) and (Ji==Jf) and (eigvec_i==eigvec_f) );

  for (size_t i=0; i<njorb; ++i)
  {
    int jmin = 0;
    if ( same_i_f ) jmin = i;
    for (size_t j=jmin; j<njorb; ++j)
    {
      if ( std::abs(obtd(i,j))>1e-7 )
      {
         densityfile << std::setw(3) << i << " " << std::setw(3) << j << " "  << std::setw(12) << std::fixed << std::setprecision(8) << obtd(i,j)  << std::endl;
      }
    }
  }

}





void ReadWrite::WriteTBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, arma::mat& tbtd)
{

  // generate all the two body states that are needed
  int Ji = settings.J2_i[J_index_i];
  int Jf = settings.J2_f[J_index_f];
  if (std::abs(Ji-Jf)>Lambda2 or Ji+Jf<Lambda2) return;

  densityfile << std::endl;
  densityfile << "-------------- TBTD ---------------------" << std::endl;
  for (size_t i=0;i<tbtd.n_rows;++i)
  {
    for (size_t j=0;j<tbtd.n_cols;++j)
    {
       if (std::abs(tbtd(i,j))>1e-7)
       densityfile << std::setw(3) << i << " " << std::setw(3) << j << " "  << std::setw(12) << std::fixed << std::setprecision(8) << tbtd(i,j) << std::endl;
    }
  }
}














// Make a list of all the *.nba files
// note there's not a 1-to-1 correspondence of *.nba files and *.prj files. Do we need both?
std::vector<std::string> ReadWrite::FindNBAFiles( std::string basename, int Zval, int Nval, std::string a_or_b )
{
  std::ostringstream ostr;
  std::ifstream testread; // for checking if files exist
  std::vector<std::string> file_list;
  if ((Zval<0) or (Nval<0))
  {
    std::cout << __func__ << " Warning: Zval = " << Zval << " and Nval = " << Nval << " ... does that seem right?" << std::endl;
  }
  for (int icode : { Zval, -Zval, Nval, -Nval} )
  {
    for (size_t iJ=0; iJ<an_code.size(); iJ++)
    {
      ostr.str("");
      ostr.str().clear();
      ostr << basename << a_or_b << an_code[iJ] << "0" << an_code[36 + icode/2];
      testread.open(ostr.str()+".nba");
      if ( not testread.fail() )
      {
        if ( std::find(file_list.begin(),file_list.end(), ostr.str() ) == file_list.end())
          file_list.push_back( ostr.str() );
        testread.close();
      }
    }
  }
  return file_list;

}


std::string ReadWrite::FindXVCFile( std::string basename, int Jtot, int Zval, int Nval )
{
  
  // Guess the name of the xvc file
    std::ostringstream ostr;
    ostr << basename << an_code[Jtot/2] << an_code[Zval] << an_code[Nval] << ".xvc";
    std::string vecfile = ostr.str();

    if ( std::ifstream(vecfile).fail() ) // that didn't work. let's try the other way...
    {
      ostr.str("");
      ostr.clear();
      ostr << basename << an_code[Jtot/2] << an_code[Nval] << an_code[Zval] << ".xvc";
      vecfile = ostr.str();

      if ( std::ifstream(vecfile).fail() )
      {
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        std::cout << "ERROR! I cant figure out what the *.xvc file should be. Exiting." << std::endl;
        std::cout << "( " << ostr.str() << " ) didnt work. abfile_base = " << basename << std::endl;
        std::cout << "I think there are " << Zval << " valence protons and " << Nval << " valence neutrons" << std::endl;
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        exit(EXIT_FAILURE) ;
      }
    }
    return vecfile;
 
}












void ReadWrite::ReadNuShellFiles( TransitionDensity& trans )
{

  double t_start;
  GetAZFromFileName();
  std::cout << "Acore, Zcore = " << settings.Acore << " " << settings.Zcore << std::endl;

  trans.Jlist_i = settings.J2_i;
  trans.Jlist_f = settings.J2_f;


  // Use the largest projection M_J that we can get away with
  trans.MJtot_i = (*std::min_element(std::begin(settings.J2_i),std::end(settings.J2_i)))%2;
  trans.MJtot_f = (*std::min_element(std::begin(settings.J2_f),std::end(settings.J2_f)))%2;


  // Sanity check: if we have odd mass number, the J should be half-integer.
  for (auto j : settings.J2_i )
  {
    if (j%2 != settings.A_i%2)
    {
      std::cout << __func__ << "  Warning settings.A_i=" << settings.A_i << " and J*2 = " << j << ".  This isn't good" << std::endl;
    }
  }
  for (auto j : settings.J2_f )
  {
    if (j%2 != settings.A_f%2)
    {
      std::cout << __func__ << "  Warning A_f=" << settings.A_f << " and J*2 = " << j << ".  This isn't good" << std::endl;
    }
  }
  

  int nvalence_protons_i  = settings.Z_i - settings.Zcore;
  int nvalence_neutrons_i = settings.A_i-settings.Z_i - (settings.Acore-settings.Zcore);
  int nvalence_protons_f  = settings.Z_f - settings.Zcore;
  int nvalence_neutrons_f = settings.A_f-settings.Z_f - (settings.Acore-settings.Zcore);
  std::vector<std::string> Afiles_i,Bfiles_i,Afiles_f,Bfiles_f;

  std::ostringstream ostr;
  std::ifstream testread; // for checking if files exist


#ifdef VERBOSE
  std::cout << __func__ << " -- finding all prj and nba files" << std::endl;
#endif
  Afiles_i = FindNBAFiles( settings.basename_vectors_i, nvalence_protons_i, nvalence_neutrons_i, "a");
  Bfiles_i = FindNBAFiles( settings.basename_vectors_i, nvalence_protons_i, nvalence_neutrons_i, "b");
  Afiles_f = FindNBAFiles( settings.basename_vectors_f, nvalence_protons_f, nvalence_neutrons_f, "a");
  Bfiles_f = FindNBAFiles( settings.basename_vectors_f, nvalence_protons_f, nvalence_neutrons_f, "b");


//  t_start = omp_get_wtime();
#ifdef VERBOSE
  std::cout << __func__ << " -- Starting loop over Jlist_i" << std::endl;
#endif


  // check to see if the list of final states is exactly the same as the list of initial states
  bool initial_final_same = true;
  for (auto& afile : Afiles_f)
  {
    if ( find( begin(Afiles_i), end(Afiles_i), afile) == end(Afiles_i) ) initial_final_same = false;
  }
  for (auto& bfile : Bfiles_f)
  {
    if ( find( begin(Bfiles_i), end(Bfiles_i), bfile) == end(Bfiles_i) ) initial_final_same = false;
  }


  for ( auto a : Afiles_i ) std::cout << a << "  ";
  for ( auto a : Bfiles_i ) std::cout << a << "  ";
  for ( auto a : Afiles_f ) std::cout << a << "  ";
  for ( auto a : Bfiles_f ) std::cout << a << "  ";
  std::cout << std::endl;


  settings.jbasis_i = JBasis( settings.basename_sps + ".sps", Afiles_i, Bfiles_i, settings.J2_i ); // Newly added...


  for (size_t iJ=0; iJ<settings.J2_i.size(); iJ++)
  {
    int Jtot = settings.J2_i[iJ];
    int NJ = settings.NJ_i[iJ];


    std::string vecfile = FindXVCFile( settings.basename_vectors_i,  Jtot, nvalence_protons_i, nvalence_neutrons_i );

    trans.nuvec_list_i.emplace_back( NuVec(Jtot) );
    trans.nuvec_list_i.back().ReadFile(vecfile, NJ);

   int imax = std::min( NJ, trans.nuvec_list_i.back().no_level) ; // how many eigenvectors with Jtot are there in the file, and how many are we interested in? 

   trans.blank_vector_i.push_back(std::vector<float>(imax, 0.));
   settings.total_number_levels_i += imax;

  }


  std::cout << "Making jbasis_f" << std::endl;
  settings.jbasis_f = JBasis( settings.basename_sps + ".sps", Afiles_f, Bfiles_f, settings.J2_f ); // Newly added...

  for (size_t iJ=0; iJ<settings.J2_f.size(); iJ++ )
  {
    int Jtot = settings.J2_f[iJ];
    int NJ = settings.NJ_f[iJ];


    std::string vecfile = FindXVCFile(settings.basename_vectors_f,  Jtot, nvalence_protons_f, nvalence_neutrons_f );

    // jesus this is still so obtuse...
    auto ji_iter = std::find( std::begin(settings.J2_i), std::end(settings.J2_i), Jtot);
    if (initial_final_same and ji_iter != end(settings.J2_i) )
    {
      auto index = ji_iter - std::begin(settings.J2_i);
      trans.nuvec_list_f.emplace_back( trans.nuvec_list_i[ index ] );
    }
    else
    {
      trans.nuvec_list_f.emplace_back( NuVec(Jtot) );
      trans.nuvec_list_f.back().ReadFile(vecfile,NJ);
    }
   int imax = std::min( NJ, trans.nuvec_list_f.back().no_level) ; // how many eigenvectors with Jtot are there in the file, and how many are we interested in? 

   trans.blank_vector_f.push_back(std::vector<float>(imax, 0.));
   settings.total_number_levels_f += imax;

  }



  // make a list of j orbits.
  settings.m_orbits = settings.jbasis_i.m_orbits;
  trans.m_orbits = settings.m_orbits;
  settings.jorbits.clear();
  for (size_t i=0;i<settings.m_orbits.size();++i )
  {
    if (settings.m_orbits[i].mj2 == -settings.m_orbits[i].j2)    settings.jorbits.push_back(i);
  }
  

  // set up our numbering scheme for 2-body kets |ab;J>
  for (size_t a=0;a<settings.jorbits.size();++a)
  {
    int ja = settings.m_orbits[settings.jorbits[a]].j2;
    for (size_t b=a; b<settings.jorbits.size();++b)
    {      
      int jb = settings.m_orbits[settings.jorbits[b]].j2;
      int Jmin = std::abs(ja-jb);
      int Jmax = ja+jb;
      for (int J2=Jmin;J2<=Jmax;J2+=2)
      {
        if (a==b and (J2%4)>0) continue;
        settings.ket_a.push_back(a);
        settings.ket_b.push_back(b);
        settings.ket_J.push_back(J2);
      }
    }
  }


}





ScalarOperator ReadWrite::ReadScalarOperator( std::string filename )
{

  ScalarOperator Op;

  auto& m_orbits = settings.m_orbits;
  auto& jorbits = settings.jorbits;
  auto& ket_a = settings.ket_a;
  auto& ket_b = settings.ket_b;
  auto& ket_J = settings.ket_J;


  auto& Op0b = Op.ZeroBody;
  auto& Op1b = Op.OneBody;
  auto& Op2b = Op.TwoBody;

  Op1b.zeros( jorbits.size(), jorbits.size() );
  Op2b.zeros( ket_a.size(), ket_a.size() );
  
  std::ifstream opfile(filename);
  if (! opfile.good() )
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "!!! ERROR: Trouble reading " << filename << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    return Op;
  }
  std::string line = "!";


  while (line.find("!") != std::string::npos)
  {
    getline(opfile,line);
    auto zb_ptr = line.find("Zero body term:");
    if ( zb_ptr != std::string::npos)
       std::istringstream( line.substr(zb_ptr + 16) ) >> Op0b;
  }

  std::istringstream iss( line );
  int dummy;
  iss >> dummy;
  for (size_t i=0; i<jorbits.size(); ++i)
  {
     iss >> Op1b(i,i);
     Op1b(i,i) *= sqrt( m_orbits[jorbits[i]].j2 + 1); // convert to a reduced matrix element
  }


  int a,b,c,d,J,T;
  double ME;
  while( opfile >> a >> b >> c >> d >> J >> T >> ME )
  {
    a--;b--;c--;d--;
    auto& oa = m_orbits[jorbits[a]];
    auto& ob = m_orbits[jorbits[b]];
    auto& oc = m_orbits[jorbits[c]];
    auto& od = m_orbits[jorbits[d]];

    if (a>b)
    {
      ME *= -(1-std::abs(oa.j2 + ob.j2 -J*2 )%4);
      std::swap(a,b);
    }
    if (c>d)
    {
      ME *= -(1-std::abs(oc.j2 + od.j2 -J*2 )%4);
      std::swap(c,d);
    }
    ME *= sqrt(2*J+1.); // Scalar ME's are stored as not-reduced. This converts them to reduced.

    if ( oa.tz2 != ob.tz2 )
    {
      if ( not (oa.n==ob.n and oa.l2==ob.l2 and oa.j2==ob.j2) )  ME /= SQRT2; // pn matrix elements are not normalized
      if ( not (oc.n==od.n and oc.l2==od.l2 and oc.j2==od.j2) )  ME /= SQRT2; // pn matrix elements are not normalized
    }

    size_t ibra=0,iket=0;
    while( ibra<ket_a.size() and not( (ket_a[ibra]==a) and (ket_b[ibra]==b) and ket_J[ibra]==2*J) ) ibra++;
    while( iket<ket_a.size() and not( (ket_a[iket]==c) and (ket_b[iket]==d) and ket_J[iket]==2*J) ) iket++;

    Op2b(ibra,iket)  += ME ; 
    Op2b(iket,ibra) = Op2b(ibra,iket);
  }


  return Op;
}










TensorOperator ReadWrite::ReadTensorOperator( std::vector<std::string> filenames)
{
  TensorOperator Op;

  auto& Op1b = Op.OneBody;
  auto& Op2b = Op.TwoBody;

  auto& jorbits = settings.jorbits;
  auto& m_orbits = settings.m_orbits;
  auto& ket_a = settings.ket_a;
  auto& ket_b = settings.ket_b;
  auto& ket_J = settings.ket_J;

  if ( filenames.size() < 1 ) return Op;
  std::string filename1b = filenames[0];

  // Read the 1-body file
  std::ifstream file1b(filename1b);

  if (not file1b.good())
  {
    std::cout << "Trouble reading file " << filename1b << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  int Rank_J, Rank_T, parity;


  while ( line.find("Rank_J") == std::string::npos)   getline(file1b, line);
  std::istringstream( line.substr( line.rfind(":")+1 ) ) >> Rank_J;
  while ( line.find("Rank_T") == std::string::npos)   getline(file1b, line);
  std::istringstream( line.substr( line.rfind(":")+1 ) ) >> Rank_T;
  while ( line.find("Parity") == std::string::npos)   getline(file1b, line);
  std::istringstream( line.substr( line.rfind(":")+1 ) ) >> parity;

  Op.Lambda = Rank_J;
  Op.Rank_T = Rank_T;
  Op.Parity = parity;

  while ( line.find("index") == std::string::npos)   getline(file1b, line);
  std::vector<int> orbits_in;
  getline(file1b, line);
  while ( line.size() > 5 and line.find("a")==std::string::npos )  // no specific reason why 5. Just looking for the empty comment line
  {
    int index,n,l,j2,tz2;
    std::istringstream(line.substr(1)) >> index >> n >> l >> j2 >> tz2;
    tz2 *= -1; // switch isospin conventions...
    // find the corresponding j-orbit
    for (size_t i=0;i<jorbits.size(); ++i)
    { // we use -tz2 because we switch isospin conventions
      if (    m_orbits[jorbits[i]].n==n  and m_orbits[jorbits[i]].l2==2*l and m_orbits[jorbits[i]].j2==j2 and m_orbits[jorbits[i]].tz2==tz2)
      {
         orbits_in.push_back(i);
      }
    }

    getline(file1b, line);
  }


  Op1b.zeros( jorbits.size(),jorbits.size() );

  getline(file1b, line); // skip final header
  int ain,bin,a,b;
  double Op_ab;
  while ( file1b >> ain >> bin >> Op_ab)
  {
    a = orbits_in[ain-1]; // fortran indexing...
    b = orbits_in[bin-1];
    int j2a = m_orbits[jorbits[a]].j2;
    int j2b = m_orbits[jorbits[b]].j2;
    int tz2a = m_orbits[jorbits[a]].tz2;
    int tz2b = m_orbits[jorbits[b]].tz2;
    int la = m_orbits[jorbits[a]].l2/2;
    int lb = m_orbits[jorbits[b]].l2/2;
    if ( ( (std::abs(j2a-j2b) > 2*Rank_J) or (j2a+j2b < 2*Rank_J) ) and std::abs(Op_ab)>1e-6 )
    {
      std::cout << __func__ << "  WARNING: mat. el. violates triangle contition for rank-" << Rank_J << " operator. -- " << ain << " " << bin << " " << Op_ab << std::endl;
      continue;
    }
    if ( (std::abs(tz2a-tz2b) != 2*Rank_T) and std::abs(Op_ab)>1e-6   )
    {
      std::cout << __func__ << "  WARNING: mat. el. has wrong isospin projection. dTz = " << Rank_T << " operator. -- " << ain << " " << bin << " " << Op_ab << std::endl;
      continue;
    }
    if ( (la+lb)%2 != (1-parity)/2 and std::abs(Op_ab)>1e-6  )
    {
      std::cout << __func__ << "  WARNING: mat. el. has wrong parity (" << parity << ")-- " << ain << " " << bin << " " << Op_ab << std::endl;
      continue;
    }
    Op1b(a,b) = Op_ab;
    Op1b(b,a) = (1 - std::abs(j2a-j2b)%4) * Op_ab; // phase factor (-1)^(ja-jb)
  }



  if (filenames.size()<2) return Op;
  std::string filename2b = filenames[1];

  // Now read the 2-body file
  std::ifstream file2b(filename2b);
  if (not file2b.good())
  {
    std::cout << "Trouble reading file " << filename2b << std::endl;
    exit(EXIT_FAILURE);
  }

  while ( line.find("Rank_J") == std::string::npos)   getline(file2b, line);
  std::istringstream( line.substr( line.rfind(":")+1 ) ) >> Rank_J;
  while ( line.find("Rank_T") == std::string::npos)   getline(file2b, line);
  std::istringstream( line.substr( line.rfind(":")+1 ) ) >> Rank_T;
  while ( line.find("Parity") == std::string::npos)   getline(file2b, line);
  std::istringstream( line.substr( line.rfind(":")+1 ) ) >> parity;

  if ( (Rank_J != Op.Lambda) or (Rank_T != Op.Rank_T) or (parity != Op.Parity) )
  {
    std::cout << __func__ << " ERROR: tensor rank mismatch between files " << filename1b << " and " << filename2b <<std::endl;
    exit(EXIT_FAILURE);
  }


  // Read in the single particle basis used in the file
  while ( line.find("index") == std::string::npos)   getline(file2b, line);
  std::vector<int> orbits_in_2b;
  getline(file2b, line);
  while ( line.size() > 5 and line.find("a")==std::string::npos )  // no specific reason why 5. Just looking for the empty comment line
  {
    int index,n,l,j2,tz2;
    std::istringstream(line.substr(1)) >> index >> n >> l >> j2 >> tz2;
    tz2 *= -1; //  we switch isospin conventions...
    // find the corresponding j-orbit
    for (size_t i=0;i<jorbits.size(); ++i)
    {
      if (    m_orbits[jorbits[i]].n==n  and m_orbits[jorbits[i]].l2==2*l and m_orbits[jorbits[i]].j2==j2 and m_orbits[jorbits[i]].tz2==tz2)
      {
         orbits_in_2b.push_back(i);
      }
    }
    getline(file2b, line);
  }
  for ( size_t i=0; i<orbits_in_2b.size(); i++)
  {
    if (orbits_in_2b[i] != orbits_in[i])
    {
      std::cout << __func__ << "  ERROR:  the 1b and 2b have differing singple particle bases. Something must be wrong. Exiting." << std::endl;
      exit(EXIT_FAILURE);
    }
  }


  Op2b.zeros( ket_J.size(), ket_J.size() );

  if ( line.find("Jab")==std::string::npos )
    getline(file2b, line); // skip final header
  int cin,din,c,d,Jab,Jcd; // already defined a,b,ain,bin above
  double Op_abcd;
  while ( file2b >> ain >> bin >> cin >> din >> Jab >> Jcd >> Op_abcd )
  {
    a = orbits_in[ain-1]; // fortran indexing...
    b = orbits_in[bin-1];
    c = orbits_in[cin-1]; // fortran indexing...
    d = orbits_in[din-1];
    int j2a = m_orbits[jorbits[a]].j2;
    int j2b = m_orbits[jorbits[b]].j2;
    int j2c = m_orbits[jorbits[c]].j2;
    int j2d = m_orbits[jorbits[d]].j2;
    int tz2a = m_orbits[jorbits[a]].tz2;
    int tz2b = m_orbits[jorbits[b]].tz2;
    int tz2c = m_orbits[jorbits[c]].tz2;
    int tz2d = m_orbits[jorbits[d]].tz2;
    int la = m_orbits[jorbits[a]].l2/2;
    int lb = m_orbits[jorbits[b]].l2/2;
    int lc = m_orbits[jorbits[c]].l2/2;
    int ld = m_orbits[jorbits[d]].l2/2;
    if ( ((std::abs(Jab-Jcd) > Rank_J) or (Jab+Jcd < Rank_J)) and std::abs(Op_abcd)>2e-6 )
    {
      std::cout << "WARNING: mat. el. violates triangle contition for rank-" << Rank_J << " operator. -- "
               << ain << " " << bin << " " << cin << " " << din << " " << Jab << " " << Jcd << " " << Op_abcd << std::endl;
      continue;
    }
    if ( (std::abs(tz2a+tz2b-tz2c-tz2d) != 2*Rank_T)  and std::abs(Op_abcd)>1e-6)
    {
      std::cout << "WARNING: mat. el. has wrong isospin projection. dTz = " << Rank_T << " operator. -- " 
               << ain << " " << bin << " " << cin << " " << din << " " << Jab << " " << Jcd << " " << Op_abcd << std::endl;
      continue;
    }
    if ( (la+lb+lc+ld)%2 != (1-parity)/2 and std::abs(Op_abcd)>1e-6)
    {
      std::cout << "WARNING: mat. el. has wrong parity -- " 
               << ain << " " << bin << " " << cin << " " << din << " " << Jab << " " << Jcd << " " << Op_abcd << std::endl;
      continue;
    }
    Jab *=2;
    Jcd *=2;
    if (a>b)
    {
      std::swap(a,b);
      Op_abcd *= -(1-std::abs(j2a+j2b-Jab)%4);
    }
    if (c>d)
    {
      std::swap(c,d);
      Op_abcd *= -(1-std::abs(j2c+j2d-Jcd)%4);
    }

    size_t ibra=0,iket=0;
    while( ibra<ket_a.size() and not( (ket_a[ibra]==a) and (ket_b[ibra]==b) and ket_J[ibra]==Jab) ) ibra++;
    while( iket<ket_a.size() and not( (ket_a[iket]==c) and (ket_b[iket]==d) and ket_J[iket]==Jcd) ) iket++;

    if (ibra >= ket_a.size() or iket >= ket_a.size())
    {
      std::cout << "trouble:  " << a << " " << b << " " << c << " " << d << " " <<Jab << " " << Jcd << " " << Op_abcd << std::endl;
      std::cout << "   ibra = " << ibra << "  iket = " << iket << "  size = " << ket_a.size() << std::endl;
    }

    Op2b(ibra,iket) = Op_abcd;
    if (ibra!=iket)
      Op2b(iket,ibra) = (1 - std::abs(Jab+Jcd )%4) * Op_abcd; // phase factor (-1)^(Jab-Jcd)
  }


  return Op;
}




DaggerOperator ReadWrite::ReadDaggerOperator( std::string filename)
{
  DaggerOperator Op;

  auto& jorbits = settings.jorbits;
  auto& m_orbits = settings.m_orbits;
  auto& ket_a = settings.ket_a;
  auto& ket_b = settings.ket_b;
  auto& ket_J = settings.ket_J;

  auto& Op_ax = Op.Op_ax;
  auto& Op_axaxa = Op.Op_axaxa;


  std::ifstream dagfile(filename);

  if (not dagfile.good())
  {
    std::cout << "Trouble reading file " << filename << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  int Rank_J, Rank_T, parity;


  // This looks a bit ugly, but all it does is read after the ":" and trim off the "/2".
  // For example,  "Rank_J : 5/2"     would be parsed to give 5.
  while ( line.find("Rank_J") == std::string::npos)   getline(dagfile, line);
  std::istringstream( line.substr(line.rfind(":")+1,line.rfind("/2")-line.rfind(":")-1 )) >> Rank_J;
  while ( line.find("Rank_T") == std::string::npos)   getline(dagfile, line);
  std::istringstream( line.substr(line.rfind(":")+1,line.rfind("/2")-line.rfind(":")-1 )) >> Rank_T;
  while ( line.find("Parity") == std::string::npos)   getline(dagfile, line);
  std::istringstream( line.substr( line.rfind(":")+1 ) ) >> parity;


  Rank_T *= -1; // switch isospin convention from imsrg to nutbar.

  Op.Lambda2 = Rank_J;
  Op.Rank_T2 = Rank_T;   
  Op.Parity = parity;

  while ( line.find("index") == std::string::npos)   getline(dagfile, line);
  std::vector<int> orbits_in;
  getline(dagfile, line);
  while ( line.size() > 5 and line.find("a")==std::string::npos )  // no specific reason why 5. Just looking for the empty comment line
  {
    int index,n,l,j2,tz2;
    std::istringstream(line.substr(1)) >> index >> n >> l >> j2 >> tz2;
    tz2 *= -1; // we switch isospin conventions for some very good reason I no longer recall
    // find the corresponding j-orbit
    for (size_t i=0;i<jorbits.size(); ++i)
    {
      if (    m_orbits[jorbits[i]].n==n  and m_orbits[jorbits[i]].l2==2*l and m_orbits[jorbits[i]].j2==j2 and m_orbits[jorbits[i]].tz2==tz2)
      {
         orbits_in.push_back(i);
      }
    }
    getline(dagfile, line);
  }

  while ( line.find(" < a || Op || 0 > ") == std::string::npos ) getline(dagfile, line) ;
  Op_ax.zeros( jorbits.size() );

  int ain,a;
  double mat_el;
  
  while ( dagfile >> ain >> mat_el) // this will fail when we hit a comment line !, which is what we want
  {
    a = orbits_in[ain-1]; // fortran indexing...
    int j2a = m_orbits[jorbits[a]].j2;
    int tz2a = m_orbits[jorbits[a]].tz2;
    int la = m_orbits[jorbits[a]].l2/2;
    if ( std::abs(mat_el)>1e-6)
    {
      if ( j2a != Rank_J   )
      {
        std::cout << __func__ << "  WARNING: mat. el. violates triangle contition for rank-" << Rank_J << " operator. -- " << ain <<  " " << mat_el << std::endl;
        continue;
      }
      if ( tz2a != Rank_T  )
      {
        std::cout << __func__ << "  WARNING: mat. el. has wrong isospin projection. dTz = " << Rank_T << " operator. -- " << ain  << " " << mat_el << std::endl;
        continue;
      }
      if ( (la)%2 != (1-parity)/2 )
      {
        std::cout << __func__ << "  WARNING: mat. el. has wrong parity (" << parity << ")-- " << ain  << " " << mat_el << std::endl;
        continue;
      }
    }
    Op_ax(a) = mat_el;
  }
  // done reading the ax part. Now onto the axaxa part

  dagfile.clear();  // forgive the filestream its failbits...
  
  while ( line.find("< ab Jab || Op || c >") == std::string::npos ) 

  {
    getline(dagfile, line) ;
  }
  Op_axaxa.zeros( ket_J.size(), jorbits.size() ) ; // Note this is not a square matrix, we have 2-body states on the left, one-body states on the right.

  int bin,cin,b,c,Jab;

  while ( dagfile >> ain >> bin >> cin >> Jab >> mat_el ) // this will fail at the end of file.
  {
    a = orbits_in[ain-1];
    b = orbits_in[bin-1];
    c = orbits_in[cin-1];
    int j2a = m_orbits[jorbits[a]].j2;
    int j2b = m_orbits[jorbits[b]].j2;
    int j2c = m_orbits[jorbits[c]].j2;
    int tz2a = m_orbits[jorbits[a]].tz2;
    int tz2b = m_orbits[jorbits[b]].tz2;
    int tz2c = m_orbits[jorbits[c]].tz2;
    int la = m_orbits[jorbits[a]].l2/2;
    int lb = m_orbits[jorbits[b]].l2/2;
    int lc = m_orbits[jorbits[c]].l2/2;
    if ( std::abs(mat_el)>1.0e-6)
    {
      if ( (std::abs(2*Jab-j2c) > Rank_J) or (2*Jab+j2c < Rank_J) )
      {
        std::cout << "WARNING: mat. el. violates triangle contition for rank-" << Rank_J << " operator. -- "
                 << ain << " " << bin << " " << cin << " " <<  Jab << " " << mat_el << std::endl;
        continue;
      }
      if ( (tz2a+tz2b-tz2c) != Rank_T  )  // we applied a minus sign to tz2 and Rank_T, so the isospin conventions should be consistent.
      {
        std::cout << "WARNING: mat. el. has wrong isospin projection. dTz = " << Rank_T << " operator. -- " 
                 << ain << " " << bin << " " << cin << " " <<  Jab <<  " " << mat_el << std::endl;
        continue;
      }
      if ( (la+lb+lc)%2 != (1-parity)/2 )
      {
        std::cout << "WARNING: mat. el. has wrong parity -- " 
                 << ain << " " << bin << " " << cin << " " <<  Jab <<  " " << mat_el << std::endl;
        continue;
      }
    }
    Jab *=2;
    if (a>b)
    {
      std::swap(a,b);
      mat_el *= -(1-std::abs(j2a+j2b-Jab)%4);
    }

    size_t ibra=0;
    while( ibra<ket_a.size() and not( (ket_a[ibra]==a) and (ket_b[ibra]==b) and ket_J[ibra]==Jab) ) ibra++;

    if (ibra >= ket_a.size())
    {
      std::cout << "trouble:  " << a << " " << b << " " << c << " " << Jab <<  " " << mat_el << std::endl;
      std::cout << "   ibra = " << ibra << " ket size = " << ket_a.size() << std::endl;
    }

    Op_axaxa(ibra,c) = mat_el;


  }



  return Op;

}


void ReadWrite::WriteLogHeader()
{
  auto& jorbits = settings.jorbits;
  auto& m_orbits = settings.m_orbits;
  auto& ket_a = settings.ket_a;
  auto& ket_b = settings.ket_b;
  auto& ket_J = settings.ket_J;
  logfile << "###  One-body states ####" << std::endl;
  logfile << " index  n    l   j2   tz2  (proton has tz=+1/2)" << std::endl;
  for (size_t i=0; i<jorbits.size(); i++)
  {
    auto& orb = m_orbits[jorbits[i]];
    logfile << std::setw(4) << i << " " << std::setw(4) << orb.n << " " << std::setw(4) << orb.l2/2 << " " << std::setw(4) << orb.j2 << " " << std::setw(4) << orb.tz2 << std::endl;
  }
  logfile << std::endl;
  logfile << "###  Two-body states ####" << std::endl;
  logfile << "   i    a    b    Jab " << std::endl;
  for (size_t i=0; i<ket_a.size(); i++)
  {
    logfile << std::setw(4) << i << " " << std::setw(4) << ket_a[i] << " " << std::setw(4) << ket_b[i] << " " << std::setw(6) << ket_J[i] << std::endl;
  }

}



void ReadWrite::WriteLog_Tensor1b( int indexJi, int ivec, int indexJf, int fvec, TensorOperator& TensorOp, arma::mat& obtd )
{
  logfile << "2Jf nJf  2Ji nJi  2Lambda = " << std::setw(3) << std::setprecision(1) << settings.J2_f[indexJf] << " " << fvec+1
  << "    " << std::setw(3) << std::setprecision(1) << settings.J2_i[indexJi] << " " << ivec+1
  << "    " << std::setw(3) << std::setprecision(1) << TensorOp.Lambda*2  << std::endl;
  logfile << "======= One Body Terms ======" << std::endl;
  logfile << std::setw(4) << "a" << " " << std::setw(4) << "b"
          << std::setw(14) << "obtd(a,b)" << " "
          << std::setw(14) << "<a||Op||b>" << " "
          << std::setw(14) << "obtd * Op" << " " 
          << std::setw(14) << "Sum obtd * Op" << std::endl;
   int nc = TensorOp.OneBody.n_cols;
   int nr = TensorOp.OneBody.n_rows;
   float runningsum = 0;
   for (int c=0;c<nc;c++)
   {
     for (int r=0;r<nr;r++)
     {
       double obme = obtd(r,c) * TensorOp.OneBody(r,c);
       runningsum += obme;
       logfile << std::setw(4) << std::fixed << c << " " << std::setw(4) << std::fixed << r
               << std::setw(14) << std::fixed << std::setprecision(6) << obtd(r,c) << " "
               << std::setw(14) << std::fixed << std::setprecision(6) << TensorOp.OneBody(r,c) << " "
               << std::setw(14) << std::fixed << std::setprecision(6) << obme << " " 
               << std::setw(14) << std::fixed << std::setprecision(6) << runningsum << std::endl;
     }
   }

}




void ReadWrite::WriteLog_Tensor2b( TensorOperator& TensorOp, arma::mat& tbtd )
{
  logfile << "======= Two Body Terms ======" << std::endl;
  logfile << std::setw(4) << "A" << " " << std::setw(4) << "B"
          << std::setw(14) << "tbtd(A,B)" << " "
          << std::setw(14) << "<A||Op||B>" << " "
          << std::setw(14) << "tbtd * Op" << " " 
          << std::setw(14) << "Sum tbtd * Op" << std::endl;
   int nc = TensorOp.TwoBody.n_cols;
   int nr = TensorOp.TwoBody.n_rows;
   float runningsum = 0;
   for (int c=0;c<nc;c++)
   {
     for (int r=0;r<nr;r++)
     {
       double tbme = tbtd(r,c) * TensorOp.TwoBody(r,c);
       runningsum += tbme;
       logfile << std::setw(4) << std::fixed << r << " " << std::setw(4) << std::fixed << c
               << std::setw(14) << std::fixed << std::setprecision(6) << tbtd(r,c) << " "
               << std::setw(14) << std::fixed << std::setprecision(6) << TensorOp.TwoBody(r,c) << " "
               << std::setw(14) << std::fixed << std::setprecision(6) << tbme << " " 
               << std::setw(14) << std::fixed << std::setprecision(6) << runningsum << std::endl;
     }
   }
   logfile << "#" << std::endl;
   logfile << "#" << std::endl;
}





void ReadWrite::WriteLog_Dagger_ax(  int indexJi, int ivec, int indexJf, int fvec, DaggerOperator& DaggerOp, arma::vec& td_ax )
{
  logfile << "2Jf nJf  2Ji nJi  2Lambda = " << std::setw(3) << std::setprecision(1) << settings.J2_f[indexJf] << " " << fvec+1
  << "    " << std::setw(3) << std::setprecision(1) << settings.J2_i[indexJi] << " " << ivec+1
  << "    " << std::setw(3) << std::setprecision(1) << DaggerOp.Lambda2  << std::endl;
  logfile << "======= Dagger ax terms ======" << std::endl;
  logfile << std::setw(4) << "a" << " " 
          << std::setw(14) << "tr.dens." << " "
          << std::setw(14) << "<a||Op||0>" << " "
          << std::setw(14) << "tr.dens. * Op" << " " 
          << std::setw(14) << "Sum" << std::endl;

  int nc = DaggerOp.Op_ax.n_rows;
  float runningsum = 0;
  for (int c=0;c<nc;c++)
  {
    double me_ax = td_ax(c) * DaggerOp.Op_ax(c);
    runningsum += me_ax;
    logfile << std::setw(4) << std::fixed << c << " "
            << std::setw(14) << std::fixed << std::setprecision(6) << td_ax(c) << " "
            << std::setw(14) << std::fixed << std::setprecision(6) << DaggerOp.Op_ax(c) << " "
            << std::setw(14) << std::fixed << std::setprecision(6) << me_ax << " " 
            << std::setw(14) << std::fixed << std::setprecision(6) << runningsum << std::endl;
  }
   logfile << "#" << std::endl;
   logfile << "#" << std::endl;

}






void ReadWrite::WriteLog_Dagger_axaxa( DaggerOperator& DaggerOp, arma::mat& td_axaxa )
{
  logfile << "====== Dagger axaxa terms ====" << std::endl;
  logfile << std::setw(4) << "A" << " " << std::setw(4) << "b"
          << std::setw(14) << "tr.dens(A,b)" << " "
          << std::setw(14) << "<A||Op||b>" << " "
          << std::setw(14) << "trdens * Op" << " " 
          << std::setw(14) << "Sum" << std::endl;

   int nc = DaggerOp.Op_axaxa.n_cols;
   int nr = DaggerOp.Op_axaxa.n_rows;
   float runningsum = 0;
   for (int c=0;c<nc;c++)
   {
     for (int r=0;r<nr;r++)
     {
       double me_axaxa = td_axaxa(r,c) * DaggerOp.Op_axaxa(r,c);
       runningsum += me_axaxa;
       logfile << std::setw(4) << std::fixed << r << " " << std::setw(4) << std::fixed << c
               << std::setw(14) << std::fixed << std::setprecision(6) << td_axaxa(r,c) << " "
               << std::setw(14) << std::fixed << std::setprecision(6) << DaggerOp.Op_axaxa(r,c) << " "
               << std::setw(14) << std::fixed << std::setprecision(6) << me_axaxa << " " 
               << std::setw(14) << std::fixed << std::setprecision(6) << runningsum << std::endl;
     }
   }
   logfile << "#" << std::endl;
   logfile << "#" << std::endl;



}





void ReadWrite::WriteScalarResults( TransitionDensity& trans, std::vector<ScalarOperator>& ScalarOps, std::vector<ScalarNME>& scalarnme )
{
  for (size_t isc=0; isc<ScalarOps.size(); ++isc)
  {
    std::ostringstream oss;
    oss << "nutbar_scalar" << isc << "_" <<  settings.basename_vectors_f.substr( settings.basename_vectors_f.find_last_of("/")+1 ) ;
//    std::string filename = "nutbar_scalar" + oss.str() + "_" + ;
    if (settings.basename_vectors_i != settings.basename_vectors_f)
    {
       oss << "_" << settings.basename_vectors_i.substr( settings.basename_vectors_i.find_last_of("/")+1 );
//       filename +  "_" + settings.basename_vectors_i.substr( settings.basename_vectors_i.find_last_of("/")+1 );
    }
    oss << ".dat";
    std::ofstream scalar_out( oss.str() ); 

    scalar_out << "############################################################################################################" << std::endl;
    scalar_out << "# initial vector basename: " << settings.basename_vectors_i << std::endl;
    scalar_out << "# final   vector basename: " << settings.basename_vectors_f << std::endl;
    scalar_out << "# Operator file: " << settings.scalar_op_files[isc] << std::endl;
    scalar_out << "# Zero body term: " << ScalarOps[isc].ZeroBody << std::endl;
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

           double Ef = trans.nuvec_list_f[indexJf].alpha[njf];
           double Ei = trans.nuvec_list_i[indexJi].alpha[nji];
           double zerobody = (nji==njf) ? ScalarOps[isc].ZeroBody : 0;
           double NME1 = scalarnme[isc].OneBody(indexJf,indexJi)(njf,nji);
           double NME2 = scalarnme[isc].TwoBody(indexJf,indexJi)(njf,nji);
           scalar_out << std::fixed      << std::setw(4) << std::setprecision(1)  << Jf << " " << std::setw(3) << njf+1 << "    "
                      << std::fixed      << std::setw(4) << std::setprecision(1)  << Ji << " " << std::setw(3) << nji+1 << " "
                      << std::fixed      << std::setw(10) << std::setprecision(3) << Ef << "  "
                      << std::fixed      << std::setw(10) << std::setprecision(3) << Ei << "  "
                      << std::scientific << std::setw(14) << std::setprecision(6) << NME1 << "  "
                      << std::scientific << std::setw(14) << std::setprecision(6) << NME2 << "  "
                      << std::scientific << std::setw(14) << std::setprecision(6) << NME1 + NME2 << "  "
                      << std::scientific << std::setw(14) << std::setprecision(6) << zerobody + NME1 + NME2 
                      << std::endl;
         }
        }
      }
    }
  }

}






void ReadWrite::WriteTensorResults( TransitionDensity& trans, TensorOperator& TensorOp, TensorNME& tensornme )
{

  std::ostringstream oss;
  if ( settings.tensor_op_files.size()<1) return;
  oss << "nutbar_tensor" << TensorOp.Lambda << "_" << settings.basename_vectors_f.substr( settings.basename_vectors_f.find_last_of("/")+1 );
  if (settings.basename_vectors_i != settings.basename_vectors_f)
  {
     oss <<  "_" << settings.basename_vectors_i.substr( settings.basename_vectors_i.find_last_of("/")+1 );
  }
  oss << ".dat";
  std::ofstream tensor_out(oss.str());
  
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

         double Ef = trans.nuvec_list_f[indexJf].alpha[njf];
         double Ei = trans.nuvec_list_i[indexJi].alpha[nji];
         double NME1 = tensornme.OneBody(indexJf,indexJi)(njf,nji);
         double NME2 = tensornme.TwoBody(indexJf,indexJi)(njf,nji);
         tensor_out << std::fixed      << std::setw(4) << std::setprecision(1)  << Jf << " " << std::setw(3) << njf+1 << "    "
                    << std::fixed      << std::setw(4) << std::setprecision(1)  << Ji << " " << std::setw(3) << nji+1 << " "
                    << std::fixed      << std::setw(10) << std::setprecision(3) << Ef << "  "
                    << std::fixed      << std::setw(10) << std::setprecision(3) << Ei << "  "
                    << std::scientific << std::setw(14) << std::setprecision(6) << NME1 << "  "
                    << std::scientific << std::setw(14) << std::setprecision(6) << NME2 << "  "
                    << std::scientific << std::setw(14) << std::setprecision(6) << NME1+NME2 << "  "
                    << std::endl;
       }
      }
    }
  }


}





void ReadWrite::WriteDaggerResults( TransitionDensity& trans, std::vector<DaggerOperator>& DaggerOps, std::vector<DaggerNME>& daggernme )
{

  for (size_t idag=0; idag<DaggerOps.size(); ++idag)
  {
    std::ostringstream oss;
    oss << "nutbar_dagger" << idag << "_" <<  settings.basename_vectors_f.substr( settings.basename_vectors_f.find_last_of("/")+1 ) ;
    if (settings.basename_vectors_i != settings.basename_vectors_f)
    {
       oss << "_" << settings.basename_vectors_i.substr( settings.basename_vectors_i.find_last_of("/")+1 );
    }
    oss << ".dat";
    std::ofstream dagger_out( oss.str() );

    dagger_out << "##########################################################################################" << std::endl;
    dagger_out << "# initial vector basename: " << settings.basename_vectors_i << std::endl;
    dagger_out << "# final   vector basename: " << settings.basename_vectors_f << std::endl;
    dagger_out << "# Operator file: " << settings.dagger_op_files[idag] << std::endl;
    dagger_out << "# Jf  nJf     Ji  nJi      Ef          Ei         <Opax>          <Opaxaxa>         <Op> " << std::endl;
    dagger_out << "###########################################################################################" << std::endl;
    
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

           double Ef = trans.nuvec_list_f[indexJf].alpha[njf];
           double Ei = trans.nuvec_list_i[indexJi].alpha[nji];
           double NME1 = daggernme[idag].ax(indexJf,indexJi)(njf,nji);
           double NME2 = daggernme[idag].axaxa(indexJf,indexJi)(njf,nji);
           dagger_out << std::fixed      << std::setw(4) << std::setprecision(1)  << Jf << " " << std::setw(3) << njf+1 << "    "
                      << std::fixed      << std::setw(4) << std::setprecision(1)  << Ji << " " << std::setw(3) << nji+1 << " "
                      << std::fixed      << std::setw(10) << std::setprecision(3) << Ef << "  "
                      << std::fixed      << std::setw(10) << std::setprecision(3) << Ei << "  "
                      << std::scientific << std::setw(14) << std::setprecision(6) << NME1 << "  "
                      << std::scientific << std::setw(14) << std::setprecision(6) << NME2 << "  "
                      << std::scientific << std::setw(14) << std::setprecision(6) << NME1+NME2 << "  "
                      << std::endl;
         }
        }
      }
    }

 }


}
















// Write out the eigenstd::vectors in the Darmstadt MBPT/NCSM format
// so it can be read in by Petr Navratil's TRDENS code
void ReadWrite::WriteEGV(TransitionDensity& trans, std::string fname)
{
  std::vector<MschemeOrbit> mscheme_orbits;

  // Count number of oscillator quanta, basically emax with e = 2n+l
  int Nshell=0;
  for (auto morbit : settings.m_orbits)
  {
    Nshell = std::max(Nshell, 2*morbit.n + morbit.l2/2);
  }
  Nshell+=1;  // if emax = 0, there is 1 shell


  // find all the core orbits
  for (int N=0;N<Nshell;++N) 
  {
    int sumN = (N+1)*(N+2)*(N+3)/3;
    for (int n=0;2*n<=N;++n)
    {
      int l2 = 2*(N-2*n);
      for (int j2=l2+1;j2>=std::max(0,l2-1);j2-=2)
      {
       for (int mj2=-j2;mj2<=j2;mj2+=2)
       {
        if ( sumN<=settings.Zcore )        mscheme_orbits.emplace_back( MschemeOrbit(n,l2,j2,mj2,+1 ) );
        if ( sumN<=settings.Acore-settings.Zcore )  mscheme_orbits.emplace_back( MschemeOrbit(n,l2,j2,mj2,-1 ) );
       }
      }
    }
  }


  // loop over the valence orbits
  int n_core_orbits = mscheme_orbits.size();
  for (auto& morbit : settings.jbasis_i.m_orbits)
  {
     mscheme_orbits.push_back ( morbit );
  }
  // calculate parity of 0hw configurations
  int N=0;
  for (N=0;(N+1)*(N+2)*(N+3)/3 < settings.Z_i;++N){};
  int parity_protons = (N%2);
  for (N=0;(N+1)*(N+2)*(N+3)/3 < settings.A_i-settings.Z_i;++N){};
  int parity_neutrons = (N%2);
  int parity = 1-2*(( parity_protons*(settings.Z_i-settings.Zcore) + parity_neutrons*(settings.A_i-settings.Acore-(settings.Z_i-settings.Zcore)))%2);
  
  std::ofstream output(fname);
  
  std::string interaction_id = "imsrg";
  int Nmax = 0;
  float hw = 20; // this should be irrelevant

  output << std::left << std::setw(10) << settings.Z_i                            << " !  number of protons" << std::endl;
  output << std::left << std::setw(10) << settings.A_i-settings.Z_i                        << " !  number of neutrons" << std::endl;
  output << std::left << std::setw(std::max(10,(int)interaction_id.size()+3)) << interaction_id << " ! interaction id" << std::endl;
  output << std::left << std::setw(10) << hw                             << " ! hbar omega" << std::endl;
  output << std::left << std::setw(10) << Nshell                         << " ! N shells " << std::endl;
  output << std::left << std::setw(10) << mscheme_orbits.size()          << " ! m-scheme orbits " << std::endl;
  output << std::left << std::setw(10) << Nmax                           << " ! Nmax " << std::endl;
  output << std::left << std::setw(10) << trans.amplitudes_i.size()      << " ! number of basis Slater determinants" << std::endl;
  output << std::left << std::setw(10) << std::showpos << parity << std::noshowpos << " ! parity " << std::endl;
  output << std::left << std::setw(10) << trans.MJtot_i                  << " ! total MJ*2" << std::endl;
  output << std::left << std::setw(10) << settings.total_number_levels_i    << " ! Number of eigenstates" << std::endl;
  
  // write out energies, J, T of eigenstates
  for (auto& nuvec : trans.nuvec_list_i)
  {
   int imax = nuvec.no_level;
   for (int ilevel=0;ilevel<imax;++ilevel)
   {
     output << std::right << std::fixed << std::setw(12) << std::setprecision(4) << nuvec.alpha[ilevel] << " " << nuvec.J2/2.0 << " " << 0 << " " << 0.0 << std::endl;
   }
  }
  
  
  output << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  output << "!!!  now list mscheme single-particle basis !!!" << std::endl;
  output << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  // write single-particle basis in mscheme
  
  
  
  int s = 1;
  for (auto& morbit : mscheme_orbits)
  {
    output << std::setw(4) << s << " " << std::setw(3) << morbit.n << " " << std::setw(3) << morbit.l2 << " " << std::setw(3) << morbit.j2 << " " << std::setw(3) << morbit.mj2 << " " << std::setw(3) << morbit.tz2 << std::endl;;
    s++;
  }
  
  output << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  output << "!!!  now list mscheme basis states and amplitudes !!!" << std::endl;
  output << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  
  
  
  // Write m-scheme basis occupations and eigenstd::vector coefficients
  size_t bits_per_word = 8*sizeof(mvec_type);
  for (auto it_amp : trans.amplitudes_i)
  {
    auto& mvec = it_amp.first;
    for (int i=0;i<n_core_orbits;++i) output << i+1 << " ";
    output << " " ;
    for (size_t i=0;i<mscheme_orbits.size()-n_core_orbits;++i)
    {
      if ( mvec[i])
          output << n_core_orbits + i+1 << " ";
    }
    output << "   ";
    for (auto& ampvec : it_amp.second)
    {
      for (auto amp : ampvec )
      {
        output << std::scientific << std::setw(13) << amp << " ";
      }
    }
    output << std::endl;
  
  }

}




void ReadWrite::WriteTRDENS_input( TransitionDensity& trans, std::string fname)
{
  std::ofstream outfile(fname);

  outfile << "T        ! specify which states to take?" << std::endl;
  outfile << "1   " << settings.total_number_levels_i << "   ! ki,nki" << std::endl;
  for (int i=1;i<=settings.total_number_levels_i;++i) outfile << i << " ";
  outfile << std::endl;
  outfile << "1   " << settings.total_number_levels_i << "   ! kf,nkf" << std::endl;
  for (int i=1;i<=settings.total_number_levels_i;++i) outfile << i << " ";
  outfile << std::endl;
  outfile << *std::max_element(std::begin(trans.Jlist_i),std::end(trans.Jlist_i)) << "      ! jtotal2max" << std::endl;
  outfile << "0                   ! irestart" << std::endl;
  outfile << "2                   ! majortot" << std::endl;
  outfile << "F                   ! irem_cal" << std::endl;
  outfile << "F                   ! formfcal" << std::endl;
  outfile << "F                   ! radial  " << std::endl;
  outfile << "F                   ! momdist " << std::endl;
  outfile << "T                   ! twobdcal" << std::endl;
  outfile << "1                   ! ipn     " << std::endl;
  outfile << "IMSRG.int_iso                 " << std::endl;
  outfile << "IMSRG_E2_1b.op_iso            " << std::endl;
  outfile << "IMSRG_E2_2b.op_iso            " << std::endl;
  outfile << "F                   ! cluster " << std::endl;
  outfile << "F                   ! antoine " << std::endl;
  outfile << "24.d0               ! hbar*Omega for antoine file" << std::endl;
  outfile << "1                   ! #init states in antoine file" << std::endl;
  outfile << "0  0 " << std::endl;
  outfile << "1                   ! #fin states in antoine file" << std::endl;
  outfile << "0  0 " << std::endl;
  outfile << "T                   ! mbpt_ncsm " << std::endl;
  outfile << "F                   ! redstick" << std::endl;
  outfile << "F                   ! mfd_james" << std::endl;
  outfile << "F                   ! threebdcal" << std::endl;
  outfile << "F                   ! NCSMC_kernels" << std::endl;
  outfile << "3                   ! num_of_interaction_files" << std::endl;
  outfile << "../vrelnp_rgm.int_H2srg-n3lo2.2_2014 " << std::endl;
  outfile << "../vrelpp_rgm.int_H2srg-n3lo2.2_2014 " << std::endl;
  outfile << "../vrelnn_rgm.int_H2srg-n3lo2.2_2014 " << std::endl;
  outfile << "F                   ! V3Nint " << std::endl;
  outfile << "14 14 14            ! N1_max,N12_max,N123_max " << std::endl;
  outfile << "chi2b3b400cD-02cE0098_srg0800ho40C_eMax14_EMax14_hwHO016.me3j_bin" << std::endl;


  outfile.close();


}







void ReadWrite::PrintOptions()
{
  std::cout << "options: [ ";
  for (std::string opt : settings.options) std::cout << opt << " ";
  std::cout << " ]" << std::endl;
}

void ReadWrite::PrintOperatorNames()
{
  std::cout << "operators: [ " ;
  for (auto ops : settings.scalar_op_files) std::cout << ops << " ";
  for (auto ops : settings.tensor_op_files) std::cout << ops << " ";
  for (auto ops : settings.dagger_op_files) std::cout << ops << " ";
  std::cout << "] " << std::endl;
}




