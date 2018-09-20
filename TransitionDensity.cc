#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <armadillo>
#include <unordered_map>
#include <vector>
#include <omp.h>

#include "TransitionDensity.hh"
#include "JMState.hh"

#define SQRT2 1.4142135623730950488

//#define VERBOSE true


// a and b are presumably skipped here because a and b are used to label proton/neutron
std::vector<char> TransitionDensity::an_code = {  '0','1','2','3','4','5','6','7','8','9',
                                             '_','-','c','d','e','f','g','h','i','j',  
                                             'k','l','m','n','o','p','q','r','s','t',  
                                             'u','v','w','x','y','z','A','B','C','D',  
                                             'E','F','G','H','I','J','K','L','M','N',  
                                             'O','P','Q','R','S','T','U','V','W','X',  
                                             '0','1','2','3','4','5','6','7','8','9',  
                                             '_','-','c'};

std::vector<std::string> TransitionDensity::periodic_table = {
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


TransitionDensity::TransitionDensity()
: total_number_levels_i(0),total_number_levels_f(0),densfile_name("none"),Acore(0),Zcore(0),same_basename_f_i(false)
{}

TransitionDensity::TransitionDensity(std::vector<int> jlist)
:total_number_levels_i(0),total_number_levels_f(0),
 Jlist_i(jlist), Jlist_f(jlist),densfile_name("none"),Acore(0),Zcore(0),same_basename_f_i(false)
{}

TransitionDensity::TransitionDensity(std::vector<int> jlist_i, std::vector<int> jlist_f)
:total_number_levels_i(0),total_number_levels_f(0),
 Jlist_i(jlist_i), Jlist_f(jlist_f),densfile_name("none"),Acore(0),Zcore(0),same_basename_f_i(false)
{}








void TransitionDensity::ReadFiles( )
{

  double t_start;
  std::ifstream testread; // for checking if files exist
  std::ostringstream ostr;
  std::vector<std::string> Afiles_i,Bfiles_i,Afiles_f,Bfiles_f;
  GetAZFromFileName(); // Make a default guess of the core
  ReadSPfile(); // If the sp file is there, use that to get the core
  std::cout << "Acore,Zcore = " << Acore << " " << Zcore << std::endl;

  int nvalence_protons_i = Z_i - Zcore;
  int nvalence_neutrons_i = A_i-Z_i - (Acore-Zcore);
  int nvalence_protons_f = Z_f - Zcore;
  int nvalence_neutrons_f = A_f-Z_f - (Acore-Zcore);


//  MJtot = (*min_element(begin(Jlist),end(Jlist)))%2;
  MJtot_i = (*std::min_element(begin(Jlist_i),end(Jlist_i)))%2;
  MJtot_f = (*std::min_element(begin(Jlist_f),end(Jlist_f)))%2;

  for (auto j : Jlist_i )
  {
    if (j%2 != A_i%2)
    {
      std::cout << "Warning A=" << A_i << " and J*2 = " << j << ".  This isn't good" << std::endl;
    }
  }
  for (auto j : Jlist_f )
  {
    if (j%2 != A_f%2)
    {
      std::cout << "Warning A=" << A_f << " and J*2 = " << j << ".  This isn't good" << std::endl;
    }
  }
  
#ifdef VERBOSE
  std::cout << "TransitionDensity::ReadFiles -- finding all prj and nba files" << std::endl;
#endif
  // Find all the prj and nba files
  int iJ;
  for (int icode : { nvalence_protons_i, -nvalence_protons_i, nvalence_neutrons_i, -nvalence_neutrons_i  } )
  {
    iJ=0;
    while (true)
    {
      ostr.str("");
      ostr.str().clear();
      ostr << basename_i << "a" << an_code[iJ] << "0" << an_code[36 + icode/2];
      testread.open(ostr.str()+".nba");
      if ( testread.fail() ) break;
      if ( find(Afiles_i.begin(),Afiles_i.end(), ostr.str() ) == Afiles_i.end())
        Afiles_i.push_back( ostr.str() );
      testread.close();
      iJ++;
    }
    iJ=0;
    while (true)
    {
      ostr.str("");
      ostr.str().clear();
      ostr << basename_i << "b" << an_code[iJ] << "0" << an_code[36 + icode/2];
      testread.open(ostr.str()+".nba");
      if ( testread.fail() ) break;
      if ( find(Bfiles_i.begin(),Bfiles_i.end(), ostr.str() ) == Bfiles_i.end())
      Bfiles_i.push_back( ostr.str() );
      testread.close();
      iJ++;
    }
  }

#ifdef VERBOSE
  std::cout << "TransitionDensity::ReadFiles -- finding all final prj and nba files" << std::endl;
#endif
  // same thing but for final state files
  for (int icode : { nvalence_protons_f, -nvalence_protons_f, nvalence_neutrons_f, -nvalence_neutrons_f  } )
  {
    iJ=0;
    while (true)
    {
      ostr.str("");
      ostr.str().clear();
//      ostr << basename_f << "a" << iJ << "0" << an_code[36 + icode/2];
      ostr << basename_f << "a" << an_code[iJ] << "0" << an_code[36 + icode/2];
      testread.open(ostr.str()+".nba");
      if ( testread.fail() ) break;
      if ( find(Afiles_f.begin(),Afiles_f.end(), ostr.str() ) == Afiles_f.end())
        Afiles_f.push_back( ostr.str() );
      testread.close();
      iJ++;
    }
    iJ=0;
    while (true)
    {
      ostr.str("");
      ostr.str().clear();
//      ostr << basename_f << "b" << iJ << "0" << an_code[36 + icode/2];
      ostr << basename_f << "b" << an_code[iJ] << "0" << an_code[36 + icode/2];
      testread.open(ostr.str()+".nba");
      if ( testread.fail() ) break;
      if ( find(Bfiles_f.begin(),Bfiles_f.end(), ostr.str() ) == Bfiles_f.end())
      Bfiles_f.push_back( ostr.str() );
      testread.close();
      iJ++;
    }
  }


  t_start = omp_get_wtime();
#ifdef VERBOSE
  std::cout << "TransitionDensity::ReadFiles -- Starting loop over Jlist_i" << std::endl;
#endif

  jbasis_i = JBasis( sps_file_name, Afiles_i, Bfiles_i, Jlist_i ); // Newly added...
  for (int Jtot : Jlist_i )
  {
    #ifdef VERBOSE
      std::cout << "TransitionDensity::ReadFiles -- adding jbasis with " << sps_file_name << "  (";
      for ( auto a : Afiles_i ) std::cout << " " << a;
      std::cout << ")  (";
      for ( auto b : Bfiles_i ) std::cout << " " << b;
      std::cout << ")  " << Jtot << " " << MJtot_i << std::endl;;
    #endif
//    jbasis_list_i.emplace_back( JBasis( sps_file_name, Afiles_i, Bfiles_i, Jtot, MJtot_i));
    #ifdef VERBOSE
      std::cout << "TransitionDensity::ReadFiles -- done with emplace_back" << std::endl;
    #endif
  
  // Guess the name of the xvc file
    ostr.str("");
    ostr.clear();
    ostr << basename_i << an_code[Jtot/2] << an_code[nvalence_protons_i] << an_code[nvalence_neutrons_i] << ".xvc";
    std::string vecfile = ostr.str();
    testread.open(vecfile);
    if ( testread.fail() )
    {
      ostr.str("");
      ostr.clear();
      ostr << basename_i << an_code[Jtot/2] << an_code[nvalence_neutrons_i] << an_code[nvalence_protons_i] << ".xvc";
      vecfile = ostr.str();
      testread.close();
      testread.open(vecfile);
      if ( testread.fail() )
      {
        std::cout << "ERROR! I cant figure out what the *.xvc file should be. Exiting." << std::endl;
        std::cout << "( " << ostr.str() << " ) didnt work. abfile_base = " << basename_i << std::endl;
        return ;
      }
    }
    testread.close();
    #ifdef VERBOSE
    std::cout << "TransitionDensity::ReadFiles -- after testread.close, adding Jtot = " << Jtot << ", vecfile = " << vecfile << std::endl;
    #endif
    nuvec_list_i.emplace_back( NuVec(Jtot) );
    #ifdef VERBOSE
    std::cout << "Begin read " << vecfile << std::endl;
    #endif
    nuvec_list_i.back().ReadFile(vecfile);
    #ifdef VERBOSE
    std::cout << " done." << std::endl;
    #endif
  }

  profiler.timer["Read_initial_std::vectors"] += omp_get_wtime() - t_start;
  t_start = omp_get_wtime();

#ifdef VERBOSE
  std::cout << "TransitionDensity::ReadFiles -- Starting loop over Jlist_f" << std::endl;
#endif

  bool initial_final_same = true;
  for (auto& afile : Afiles_f)
  {
    if ( find( begin(Afiles_i), end(Afiles_i), afile) == end(Afiles_i) ) initial_final_same = false;
  }
  for (auto& bfile : Bfiles_f)
  {
    if ( find( begin(Bfiles_i), end(Bfiles_i), bfile) == end(Bfiles_i) ) initial_final_same = false;
  }

  
  jbasis_f = JBasis( sps_file_name, Afiles_f, Bfiles_f, Jlist_f ); // Newly added...
  for (int Jtot : Jlist_f )
  {
    auto ji_iter = find( begin(Jlist_i), end(Jlist_i), Jtot);
//    if (initial_final_same and ji_iter != end(Jlist_i) )
//    {
//      jbasis_list_f.emplace_back( jbasis_list_i[ ji_iter-begin(Jlist_i) ] );
//    }
//    else
//    {
//      jbasis_list_f.emplace_back( JBasis( sps_file_name, Afiles_f, Bfiles_f, Jtot, MJtot_f));
//    }
  
  // Guess the name of the xvc file
    ostr.str("");
    ostr.clear();
    ostr << basename_f << an_code[Jtot/2] << an_code[nvalence_protons_f] << an_code[nvalence_neutrons_f] << ".xvc";
    std::string vecfile = ostr.str();
    testread.open(vecfile);
    if ( testread.fail() )
    {
      ostr.str("");
      ostr.clear();
      ostr << basename_f << an_code[Jtot/2] << an_code[nvalence_neutrons_f] << an_code[nvalence_protons_f] << ".xvc";
      vecfile = ostr.str();
      testread.close();
      testread.open(vecfile);
      if ( testread.fail() )
      {
        std::cout << "ERROR! I cant figure out what the *.xvc file should be. Exiting." << std::endl;
        std::cout << "( " << ostr.str() << " ) didnt work. abfile_base = " << basename_f << std::endl;
        return ;
      }
    }
    testread.close();
    if (initial_final_same and ji_iter != end(Jlist_i) )
    {
      nuvec_list_f.emplace_back( nuvec_list_i[ ji_iter-begin(Jlist_i) ] );
    }
    else
    {
      nuvec_list_f.emplace_back( NuVec(Jtot) );
      nuvec_list_f.back().ReadFile(vecfile);
    }
  }

  profiler.timer["Read_final_std::vectors"] += omp_get_wtime() - t_start;


//  m_orbits = jbasis_list_i[0].nubasis_a.m_orbits;
//  m_orbits = jbasis_list_i[0].nubasis.m_orbits;
//  m_orbits = jbasis_list_i[0].m_orbits;
  m_orbits = jbasis_i.m_orbits;
  jorbits.clear();
  for (size_t i=0;i<m_orbits.size();++i )
  {
    if (m_orbits[i].mj2 == -m_orbits[i].j2) jorbits.push_back(i);
  }
  

  Nshell=0;
  for (auto morbit : m_orbits)
  {
    Nshell = std::max(Nshell, 2*morbit.n + morbit.l2/2);
  }
  Nshell++;

  for (size_t ivec=0; ivec<nuvec_list_i.size(); ++ivec)
  {
   int imax = nuvec_list_i[ivec].no_level;
   if ( max_states_per_J_i.find(nuvec_list_i[ivec].J2) != max_states_per_J_i.end() ) imax = std::min(imax,max_states_per_J_i[nuvec_list_i[ivec].J2]);
    blank_vector_i.push_back(std::vector<float>(imax, 0.));
    total_number_levels_i += blank_vector_i.back().size();
  }
  for (size_t ivec=0; ivec<nuvec_list_f.size(); ++ivec)
  {
   int imax = nuvec_list_f[ivec].no_level;
   if ( max_states_per_J_f.find(nuvec_list_f[ivec].J2) != max_states_per_J_f.end() ) imax = std::min(imax,max_states_per_J_f[nuvec_list_f[ivec].J2]);
    blank_vector_f.push_back(std::vector<float>(imax, 0.));
    total_number_levels_f += blank_vector_f.back().size();
  }

}

/// Calculate the Mscheme wave functions for each eigenstate.
/// If the final and initial states are the same, then don't calculate it twice.
void TransitionDensity::CalculateMschemeAmplitudes()
{
  std::cout << "start CalculateMschemeAmplitudes" << std::endl;
//  CalculateMschemeAmplitudes_fi( nuvec_list_i, jbasis_list_i, max_states_per_J_i, blank_vector_i, amplitudes_i);
  CalculateMschemeAmplitudes_fi( nuvec_list_i, jbasis_i, max_states_per_J_i, blank_vector_i, amplitudes_i);
  bool same_f_i = true;
  if (basename_i != basename_f or Jlist_i.size() != Jlist_f.size()) same_f_i = false;
  if (same_f_i)
  {
    for (size_t iJ=0; iJ<Jlist_i.size(); ++iJ)
    {
      if (Jlist_i[iJ] != Jlist_f[iJ])
      {
        same_f_i = false;
        break;
      }
      if (max_states_per_J_i[Jlist_i[iJ]] < max_states_per_J_f[Jlist_f[iJ]])
      {
        same_f_i = false;
        break;
      }
    }
  }
  if (same_f_i) amplitudes_f = amplitudes_i;
  else
//     CalculateMschemeAmplitudes_fi( nuvec_list_f, jbasis_list_f, max_states_per_J_f, blank_vector_f, amplitudes_f);
     CalculateMschemeAmplitudes_fi( nuvec_list_f, jbasis_f, max_states_per_J_f, blank_vector_f, amplitudes_f);
}


// vec   : labels a sub-block of good J,pi
// state : basis state with good J
// level : eigenstate of the Hamiltonian
//
// if we call a level |J> and a state |j> and an m-scheme determinant |m>,
// this function evaluates   |J> = sum_{j,m} |m><m|j><j|J>.
//
//void TransitionDensity::CalculateMschemeAmplitudes_fi(std::vector<NuVec>& nuvec_list, std::vector<JBasis>& jbasis_list, std::unordered_map<int,int>& max_states_per_J, std::vector<std::vector<float>>& blank_vector, std::unordered_map< key_type, std::vector<std::vector<float>> >& amplitudes)
void TransitionDensity::CalculateMschemeAmplitudes_fi(std::vector<NuVec>& nuvec_list, JBasis& jbasis, std::unordered_map<int,int>& max_states_per_J, std::vector<std::vector<float>>& blank_vector, std::unordered_map< key_type, std::vector<std::vector<float>> >& amplitudes)
{
  double t_start = omp_get_wtime();
  int nthreads = omp_get_max_threads();
  std::cout << "max_states_per_J:  ";
  for (auto m : max_states_per_J) std::cout << "J = " << m.first/2.0 << " => " << m.second << "states;  ";
  std::cout << std::endl;
  std::cout << "Jvals :  ";
  for (auto& nv : nuvec_list) std::cout << nv.J2/2 << " ";
  std::cout << std::endl;
  // Choose the MJ2 value equal to the minimum J in the basis. I do this because m-scheme dimension grows with decreasing M.
  int MJ2 = (*std::min_element(begin(jbasis.basis_states),end(jbasis.basis_states),jbasis.basis_states.value_comp())).first % 2 ;
  for (size_t ivec=0; ivec<nuvec_list.size(); ++ivec)
  {
   const auto& nuvec = nuvec_list[ivec];
   int J2 = nuvec.J2;
//   const auto& jbasis = jbasis_list[ivec];
   std::cout << "Naive m-scheme dimension = " << jbasis.GetNaiveMschemeDimension(J2) << "  ";
   int imax = nuvec.no_level;
  #ifdef VERBOSE
   std::cout << "ivec = " << ivec << std::endl;
   std::cout << "no_state = " << nuvec.no_state << std::endl;
   std::cout << "no_level = " << nuvec.no_level << "  max_states_per_J = " << max_states_per_J[nuvec.J2] << std::endl;
  #endif
   if ( max_states_per_J.find(nuvec.J2) != max_states_per_J.end() ) imax = std::min(imax,max_states_per_J[nuvec.J2]);
   std::cout << "imax = " << imax << std::endl;


//   if (nuvec.no_state > jbasis.basis_states.size() )
   if (nuvec.no_state > jbasis.basis_states[J2].size() )
   {
     std::cout << "ERROR -- TransitionDensity::CalculateMschemeAmplitudes_fi -- nuvec.no_state = " << nuvec.no_state << ",  basis_states.size() = " << jbasis.basis_states[J2].size()
          << "   ivec = " << ivec << "  J = " << nuvec.J2/2 << std::endl;
     return;
   }
  
   double t_start_inner = omp_get_wtime();
   std::vector<std::unordered_map< key_type, std::vector<float> >> local_amplitudes( nthreads );
   // give each thread a local copy of the amplitudes
   // loop over J-coupled basis states
   // I make this loop parallel because there are typically lots of basis states,
   // while there often may only be one ilevel and one ivec.
   // There is probably a more clever way to do this, since with a large calculation
   // using lots of threads will chew up lots of memory.
   #pragma omp parallel for schedule(dynamic,1)
   for (int istate=0;istate<nuvec.no_state;++istate)
   {
     int thread_num = omp_get_thread_num();
//     double t_getstate = omp_get_wtime();
//     const JMState jmst = jbasis.GetBasisState(istate);
     const JMState jmst = jbasis.GetBasisState(istate,J2,MJ2);  // can we be smarter about picking an initial projection?
//     profiler.timer["GetBasisState"] += omp_get_wtime() - t_getstate;
     for (auto& it_mstate : jmst.m_coefs)
     {
       auto& key = it_mstate.first;
       const float& m_coef = it_mstate.second;

       if( local_amplitudes[thread_num].find(key) == local_amplitudes[thread_num].end() ) local_amplitudes[thread_num][key] = std::vector<float>(imax,0.);
       for (size_t ilevel=0;ilevel<imax;++ilevel) local_amplitudes[thread_num][key][ilevel] += m_coef * nuvec.coefT[ilevel][istate];
     }
   }
   profiler.timer["amplitudes_calc"] += omp_get_wtime() - t_start_inner;
   t_start_inner = omp_get_wtime();

   // Now accumulate amplitudes from all threads
   {
     for (int ith=0;ith<nthreads;++ith)
     {
       for (auto& it_amp : local_amplitudes[ith] )
       {
         auto& key = it_amp.first;
         if( amplitudes.find(key) == amplitudes.end() ) amplitudes[key] = blank_vector;
         for (size_t ilevel=0;ilevel<imax;++ilevel) amplitudes[key][ivec][ilevel] +=  it_amp.second[ilevel];
       }
     }
   }
   profiler.timer["amplitudes_accumulate"] += omp_get_wtime() - t_start_inner;
  }

  profiler.timer["CalculateMschemeAmplitudes"] += omp_get_wtime() - t_start;
}




/// Calculate a single one body transition density
double TransitionDensity::OBTD(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a, int m_index_b, int Lambda2 )
{

  int J2i = Jlist_i[J_index_i];
  int J2f = Jlist_f[J_index_f];

  if ((J2i+Lambda2 < J2f) or (std::abs(J2i-Lambda2)>J2f)) return 0;
 
 
  int Mi = J2i%2;
  int Mf = J2f%2;
  int mu = Mf-Mi;

  int j2_a = m_orbits[m_index_a].j2;
  int j2_b = m_orbits[m_index_b].j2;

  m_index_a += (j2_a - m_orbits[m_index_a].mj2)/2;
  m_index_b += (j2_b - m_orbits[m_index_b].mj2)/2;

  
  // find m-scheme orbits so that m_a = m_b, which will work for mu=0
  double obd = 0;
  int ma_min = std::max(-j2_a,mu-j2_b);
  int ma_max = std::min(j2_a,mu+j2_b);

  std::vector<key_type> keys_i;
  std::vector<double> amp_vec_i;
  for (auto& it_amp : amplitudes_i )
  {
     double amp_i = it_amp.second[J_index_i][eigvec_i];
     if (std::abs(amp_i)<1e-7) continue;
     auto& key = it_amp.first;
     keys_i.push_back( key );
     amp_vec_i.push_back(amp_i);
  }



  double clebsch_fi = CG(J2i,Mi,Lambda2,mu,J2f,Mf);
  if (std::abs(clebsch_fi)<1e-9)
  {
     clebsch_fi = CG(J2i,Mi+2,Lambda2,mu-2,J2f,Mf);
     if (std::abs(clebsch_fi)<1e-9)
     {
        std::cout << " ERROR:    Still got zero Clebsch" << std::endl;
        return 0;
     }
     else
     {
       Jplus(keys_i, amp_vec_i, J2i, Mi);
       Mi += 2;
       mu -=2;
     }
  }


  for ( int ma=ma_min;ma<=ma_max;ma+=2)
  {
    int mb = ma - mu;
    int ia = m_index_a - ( j2_a -ma )/2;
    int ib = m_index_b - ( j2_b -mb )/2;


    // convention: tilded destruction operator b~(m) = (-1)**(jb + mb) b(-m)
    //                                         b(m)  = (-1)**(jb -mb) b~(-m)
    int phase_b = (1-(j2_b-mb)%4);
    double clebsch = CG(j2_a,ma,j2_b,-mb,Lambda2,mu) ;

    for (size_t iamp=0; iamp<amp_vec_i.size(); ++iamp)
    {
      auto& key = keys_i[iamp];
      auto& amp_i = amp_vec_i[iamp];

      if ( not key[ib] ) continue;
      auto new_key = key;
      new_key.set(ib,0);
      if ( new_key[ia] ) continue;
      new_key.set(ia,1);
      if (amplitudes_f.find(new_key) == amplitudes_f.end() ) continue;


      int phase_ladder = 1;
      for (int iphase=std::min(ia,ib)+1;iphase<std::max(ia,ib);++iphase) if( key[iphase]) phase_ladder *=-1;
      double amp_f = amplitudes_f[new_key][J_index_f][eigvec_f];

      obd += clebsch * amp_i * amp_f * phase_ladder * phase_b;
    }
  }
  
  obd *= sqrt((J2f+1.)/(Lambda2+1)) / clebsch_fi;

  return obd;

}



/// Calculate a single two body transition density
double TransitionDensity::TBTD(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a, int m_index_b, int m_index_c, int m_index_d, int J2ab, int J2cd, int Lambda2 )
{

  int J2i = Jlist_i[J_index_i];
  int J2f = Jlist_f[J_index_f];
//  int mu = Lambda2%2;
  int Mi = J2i%2;
  int Mf = J2f%2;
  int mu = Mf-Mi;
  
  int j2_a = m_orbits[m_index_a].j2;
  int j2_b = m_orbits[m_index_b].j2;
  int j2_c = m_orbits[m_index_c].j2;
  int j2_d = m_orbits[m_index_d].j2;

  // start out with the maximally-projected m state
  m_index_a += (j2_a - m_orbits[m_index_a].mj2)/2;
  m_index_b += (j2_b - m_orbits[m_index_b].mj2)/2;
  m_index_c += (j2_c - m_orbits[m_index_c].mj2)/2;
  m_index_d += (j2_d - m_orbits[m_index_d].mj2)/2;


  // check some triangle conditions
  if ((J2f+Lambda2 < J2i) or (std::abs(J2f-Lambda2)>J2i)) return 0;
  if ((j2_a+j2_b<J2ab) or (std::abs(j2_a-j2_b)>J2ab)) return 0;
  if ((j2_c+j2_d<J2cd) or (std::abs(j2_c-j2_d)>J2cd)) return 0;
  if ((J2ab+J2cd < Lambda2) or (std::abs(J2ab-J2cd)>Lambda2)) return 0;

  
//  std::vector<std::vector<mvec_type>> keys_i;
  std::vector<key_type> keys_i;
  std::vector<double> amp_vec_i;
  for (auto& it_amp : amplitudes_i )
  {
     double amp_i = it_amp.second[J_index_i][eigvec_i];
     if (std::abs(amp_i)<1e-7) continue;
     auto& key = it_amp.first;
     keys_i.push_back( key );
     amp_vec_i.push_back(amp_i);
  }

  double clebsch_fi = CG(J2i,Mi,Lambda2,mu,J2f,Mf);
  if (std::abs(clebsch_fi)<1e-9)
  {
     clebsch_fi = CG(J2i,Mi+2,Lambda2,mu-2,J2f,Mf);
     if (std::abs(clebsch_fi)<1e-9)
     {
        std::cout << " ERROR:    Still got zero Clebsch" << std::endl;
        std::cout << "           J2i=" << J2i << " M2i=" << Mi+2 << " Lambda2=" << Lambda2 << " mu=" << mu-2 << " J2f=" << J2f << " M2f=" << Mf << std::endl;
        return 0;
     }
     else
     {
       Jplus(keys_i, amp_vec_i, J2i, Mi);
       Mi+=2;
       mu-=2;
     }
  }


  // Restricted sum a<=b c<=d gives 1/[(1+delta_ab)(1+delta_cd)], while
  // using normalized TBMEs gives sqrt[ (1+delta_ab)(1+delta_cd) ].
  // The additional factor of 2 comes from being able to limit ma < mb or mc < md
  double norm = 1;
  if (m_index_a==m_index_b) norm *= SQRT2;  
  if (m_index_c==m_index_d) norm *= SQRT2;  

  double tbd = 0;
  int Mab_min = std::max(-J2ab,mu-J2cd);
  int Mab_max = std::min( J2ab,mu+J2cd);
  for (int Mab=Mab_min;Mab<=Mab_max;Mab+=2)
  {
    int Mcd = mu - Mab;
    int phasecd = (1- std::abs(J2cd - Mcd)%4); // phase from getting rid of the tildes
    double clebsch_abcd = CG(J2ab,Mab,J2cd,Mcd,Lambda2,mu);
    if ( std::abs(clebsch_abcd)<1e-7) continue;
    int ma_min = std::max(-j2_a, Mab-j2_b);
    int ma_max = std::min(j2_a, Mab+j2_b);
    if (m_index_a==m_index_b) ma_max = std::min(ma_max, Mab/2);

    for (int ma=ma_min;ma<=ma_max;ma+=2)
    {
      int mb=Mab-ma;
      double clebsch_ab = CG(j2_a,ma, j2_b,mb, J2ab,Mab);
      if ( std::abs(clebsch_ab)<1e-7) continue;
      int ia = m_index_a - ( j2_a -ma )/2;
      int ib = m_index_b - ( j2_b -mb )/2;
      if (ib==ia) continue;
      int mc_min = std::max(-j2_c, Mab-j2_d);
      int mc_max = std::min(j2_c, Mab+j2_d);
      if (m_index_c == m_index_d) mc_max = std::min(j2_c, Mab/2);

      for (int mc=mc_min;mc<=mc_max;mc+=2)
      {
        int md = -Mcd-mc;
        double clebsch_cd = CG(j2_c,mc, j2_d, md, J2cd,-Mcd); // Mcd = -Mab, and another (-) comes from getting rid of the tildes
        if ( std::abs(clebsch_cd)<1e-7) continue;
        int ic = m_index_c - ( j2_c -mc )/2;
        int id = m_index_d - ( j2_d -md )/2;
        if (ic==id) continue;


        for ( size_t iamp=0; iamp<amp_vec_i.size(); ++iamp )
        {
          auto& key = keys_i[iamp];
          auto& amp_i = amp_vec_i[iamp];

          if (not (key[ic] && key[id]) ) continue;
          auto new_key = key;
          new_key.set(ic,0).set(id,0); // remove particles from c and then from d  (d-c-)
          if ( new_key[ia] || new_key[ib]) continue;
          new_key.set(ib,1).set(ia,1); // add particles to b and then to a  (a+b+)

          auto iter_newkey = amplitudes_f.find(new_key);
          if (iter_newkey == amplitudes_f.end() 
            or iter_newkey->second.size() < J_index_f
            or iter_newkey->second[J_index_f].size() < eigvec_f) continue;

          double amp_f = iter_newkey->second[J_index_f][eigvec_f];

          // pick up a phase from commuting the ladder operators
          int phase_ladder = (ia>ib xor id<ic) ? -1 : 1;
          for (int iphase = std::min(ic,id)+1;iphase<std::max(ic,id);++iphase)  if( key[iphase]) phase_ladder *=-1;
          for (int iphase = std::min(ia,ib)+1;iphase<std::max(ia,ib);++iphase)  if( new_key[iphase]) phase_ladder *=-1;

          tbd += amp_i * amp_f * clebsch_abcd * clebsch_ab * clebsch_cd * phasecd * phase_ladder;
        }
      }
    }
  }

  tbd *= sqrt((J2f+1.)/(Lambda2+1.)) / clebsch_fi * norm;
  return tbd;

}



double TransitionDensity::TD_ax(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a )
{

  int J2i = Jlist_i[J_index_i];
  int J2f = Jlist_f[J_index_f];

  int j2_a = m_orbits[m_index_a].j2;
  m_index_a += (j2_a - m_orbits[m_index_a].mj2)/2;

  if ((J2i+j2_a < J2f) or (std::abs(J2i-j2_a)>J2f)) return 0;
 
  int Mi = J2i%2;
  int Mf = J2f%2;
//  int mu = Mf-Mi;
  int ma = Mf-Mi;


  
  double td = 0;
//  int ma_min = std::max(-j2_a,mu-j2_b);
//  int ma_max = std::min(j2_a,mu+j2_b);

  std::vector<key_type> keys_i;
  std::vector<double> amp_vec_i;
  for (auto& it_amp : amplitudes_i )
  {
     double amp_i = it_amp.second[J_index_i][eigvec_i];
     if (std::abs(amp_i)<1e-7) continue;
     auto& key = it_amp.first;
     keys_i.push_back( key );
     amp_vec_i.push_back(amp_i);
  }



//  double clebsch_fi = CG(J2i,Mi,Lambda2,mu,J2f,Mf);
  double clebsch_fi = CG(J2i,Mi,j2_a,ma,J2f,Mf);
  if (std::abs(clebsch_fi)<1e-9)
  {
     clebsch_fi = CG(J2i,Mi+2,j2_a,ma-2,J2f,Mf);
     if (std::abs(clebsch_fi)<1e-9)
     {
        std::cout << " ERROR:    Still got zero Clebsch" << std::endl;
        return 0;
     }
     else
     {
       Jplus(keys_i, amp_vec_i, J2i, Mi);
       Mi += 2;
       ma -=2;
     }
  }


//  for ( int ma=ma_min;ma<=ma_max;ma+=2)
//  {
//    int mb = ma - mu;
    int ia = m_index_a - ( j2_a -ma )/2;
//    int ib = m_index_b - ( j2_b -mb )/2;


    // convention: tilded destruction operator b~(m) = (-1)**(jb + mb) b(-m)
    //                                         b(m)  = (-1)**(jb -mb) b~(-m)
//    int phase_b = (1-(j2_b-mb)%4);
//    double clebsch = CG(j2_a,ma,j2_b,-mb,Lambda2,mu) ;

    // loop through initial state amplitudes
    for (size_t iamp=0; iamp<amp_vec_i.size(); ++iamp)
    {
      auto& key = keys_i[iamp];
      auto& amp_i = amp_vec_i[iamp];

      if (key[ia]) continue; // if orbit a is occupied in the initial configuration, we get zero

//      if ( not key[ib] ) continue;
      auto new_key = key;
//      new_key.set(ia,1);
//      if ( new_key[ia] ) continue;
      new_key.set(ia,1);
      if (amplitudes_f.find(new_key) == amplitudes_f.end() ) continue;


      int phase_ladder = 1 - 2*(ia%2); // pick up a phase from commuting the a+ over to its place

//      for (int iphase=std::min(ia,ib)+1;iphase<std::max(ia,ib);++iphase) if( key[iphase]) phase_ladder *=-1;

      double amp_f = amplitudes_f[new_key][J_index_f][eigvec_f];

      td += amp_i * amp_f * phase_ladder;
    }
//  }
  
  td *= - sqrt(J2f+1.) / clebsch_fi;  // minus sign from Wigner-Eckart convention of a phase (-1)^{2*lambda}, and lambda=ja is half-integer in this case.


  return td;
}











double TransitionDensity::TD_axaxa(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a, int m_index_b, int m_index_c, int J2ab )
{
  double td = 0;

  return td;
}




arma::mat TransitionDensity::CalcOBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2)
{

  double t_start = omp_get_wtime();

  int Ji = Jlist_i[J_index_i];
  int Jf = Jlist_f[J_index_f];
  size_t njorb = jorbits.size();
  arma::mat obtd(njorb,njorb,arma::fill::zeros);
//  std::ofstream densout( densfile_name, std::ios::app );
//  if ( densfile_name != "none")
//  {
//     densout << std::endl;
//     densout << "Jf nJf  Ji nJi  Lambda = " << std::setw(3) << std::setprecision(1) << Jf*0.5 << " " << eigvec_f+1
//             << "    " << std::setw(3) << std::setprecision(1) << Ji*0.5 << " " << eigvec_i+1
//             << "    " << std::setw(3) << std::setprecision(1) << Lambda2*0.5  << std::endl;
//     densout << "-------------- OBTD ---------------------" << std::endl;
//  }

  for (size_t i=0; i<njorb; ++i)
  {
    int j2i = m_orbits[jorbits[i]].j2;
    int jmin = 0;
    if (basename_i==basename_f and Ji==Jf and eigvec_i==eigvec_f) jmin = i;
    for (size_t j=jmin; j<njorb; ++j)
    {
      obtd(i,j) = OBTD( J_index_i, eigvec_i, J_index_f, eigvec_f, jorbits[i], jorbits[j], Lambda2);
      if (Lambda2 == 0)  obtd(i,j) /= sqrt( Ji+1 );

//      if (densfile_name != "none" and std::abs(obtd(i,j))>1e-7)
//      {
//         densout << std::setw(3) << i << " " << std::setw(3) << j << " "  << std::setw(12) << std::fixed << std::setprecision(8) << obtd(i,j)  << std::endl;
//      }
      if (basename_i==basename_f and Ji==Jf and eigvec_i==eigvec_f)
      {
        int j2j = m_orbits[jorbits[j]].j2;
        obtd(j,i) = (1-std::abs(j2j-j2i)%4) * obtd(i,j);
      }
    }
  }

//  densout.close();
  profiler.timer["CalcOBTD"] += omp_get_wtime() - t_start;
  return obtd;

}






arma::mat TransitionDensity::CalcTBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2)
{

  double t_start = omp_get_wtime();
  // generate all the two body states that are needed
   SetupKets();

  arma::mat tbtd(ket_J.size(), ket_J.size(), arma::fill::zeros);
  int Ji = Jlist_i[J_index_i];
  int Jf = Jlist_f[J_index_f];
  if (std::abs(Ji-Jf)>Lambda2 or Ji+Jf<Lambda2) return tbtd;


  #pragma omp parallel for schedule(dynamic,1)
  for (size_t ibra=0;ibra<ket_J.size();++ibra)
  {
    int a = ket_a[ibra];
    int b = ket_b[ibra];
    int J2ab = ket_J[ibra];
    size_t iket_min = ((basename_i==basename_f) and (Ji == Jf) and (eigvec_i==eigvec_f)) ? ibra : 0;
    for (size_t iket=iket_min;iket<ket_J.size();++iket)
    {
      int c = ket_a[iket];
      int d = ket_b[iket];
      int J2cd = ket_J[iket];
      tbtd(ibra,iket) = TBTD(  J_index_i,  eigvec_i,  J_index_f,  eigvec_f,
                               jorbits[a], jorbits[b],  jorbits[c], jorbits[d],
                                                          J2ab,  J2cd,  Lambda2 );

      // if final == initial, only calculate half of the matrix
      if  ((basename_i==basename_f) and (Ji == Jf) and (eigvec_i==eigvec_f))
      {
        tbtd(iket,ibra) = tbtd(ibra,iket) * (1-std::abs(J2ab-J2cd)%4); 
      }
    }
  }

  if (Lambda2==0) tbtd /= sqrt( Ji+1.);
//  if ( densfile_name != "none" )
//  {
//    std::ofstream densout(densfile_name, std::ios::app);
//    densout << std::endl;
//    densout << "-------------- TBTD ---------------------" << std::endl;
//    for (size_t i=0;i<tbtd.n_rows;++i)
//    {
//      for (size_t j=0;j<tbtd.n_cols;++j)
//      {
//         if (std::abs(tbtd(i,j))>1e-7)
//         densout << std::setw(3) << i << " " << std::setw(3) << j << " "  << std::setw(12) << std::fixed << std::setprecision(8) << tbtd(i,j) << std::endl;
//      }
//    }
//    densout.close();
//  }

  profiler.timer["CalcTBTD"] += omp_get_wtime() - t_start;
  return tbtd;
}



arma::vec TransitionDensity::CalcTransitionDensity_ax( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f)
{

  double t_start = omp_get_wtime();

  int Ji = Jlist_i[J_index_i];
  int Jf = Jlist_f[J_index_f];
  size_t njorb = jorbits.size();
  arma::vec td(njorb,arma::fill::zeros);
  std::ofstream densout( densfile_name, std::ios::app );
  if ( densfile_name != "none")
  {
     densout << std::endl;
     densout << "Jf nJf  Ji nJi = " << std::setw(3) << std::setprecision(1) << Jf*0.5 << " " << eigvec_f+1
             << "    " << std::setw(3) << std::setprecision(1) << Ji*0.5 << " " << eigvec_i+1
             << std::endl;
     densout << "-------------- <f||a+||i> --------------------" << std::endl;
  }

  for (size_t i=0; i<njorb; ++i)
  {
    int j2i = m_orbits[jorbits[i]].j2;
//    int jmin = 0;
//    if (basename_i==basename_f and Ji==Jf and eigvec_i==eigvec_f) jmin = i;
//    for (size_t j=jmin; j<njorb; ++j)
//    {
      td(i) = TD_ax( J_index_i, eigvec_i, J_index_f, eigvec_f, jorbits[i] );
//      if (Lambda2 == 0)  obtd(i,j) /= sqrt( Ji+1 );

      if (densfile_name != "none" and std::abs(td(i))>1e-7)
      {
         densout << std::setw(3) << i << " "  << std::setw(12) << std::fixed << std::setprecision(8) << td(i)  << std::endl;
      }

  }

  densout.close();
  profiler.timer["CalcTD_ax"] += omp_get_wtime() - t_start;
  return td;

}




arma::mat TransitionDensity::CalcTransitionDensity_axaxa( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f)
{
  arma::mat td;

  return td;

}





arma::mat TransitionDensity::ReadOBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, std::string fname)
{
  bool found_it = false;
  int Ji = Jlist_i[J_index_i];
  int Jf = Jlist_f[J_index_f];
  size_t njorb = jorbits.size();
  std::ifstream densfile(fname);
  if (not densfile.good())
  {
   std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
   std::cout << "!! TransitionDensity::ReadOBTD -- trouble reading " << fname  << std::endl;
   std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
   exit(EXIT_FAILURE);
  }
  std::string line;
  std::ostringstream line_to_find;
  line_to_find <<  "Jf nJf  Ji nJi  Lambda = "   << std::setw(3) << std::setprecision(1) << Jf*0.5 << " " << eigvec_f+1
          << "    " << std::setw(3) << std::setprecision(1) << Ji*0.5 << " " << eigvec_i+1
          << "    " << std::setw(3) << std::setprecision(1) << Lambda2*0.5;


  while ( getline( densfile, line) )
  {
     if( line == line_to_find.str() )
     {
       found_it = true;
       break;
     }
  }
  while ( getline( densfile, line) )
  {
     if( line == "-------------- OBTD ---------------------" ) break;
  }
  arma::mat obtd(njorb,njorb,arma::fill::zeros);
  int i,j;
  double obd;
  while ( line.size() > 2 )
  {
    getline( densfile, line );
    std::istringstream(line) >> i >> j >> obd;
    obtd(i,j) = obd;
    if (basename_i==basename_f and Ji==Jf and eigvec_i==eigvec_f)
    {
      int j2i = m_orbits[jorbits[i]].j2;
      int j2j = m_orbits[jorbits[j]].j2;
      obtd(j,i) = (1-std::abs(j2j-j2i)%4) * obtd(i,j);
    }

  }
  densfile.close();
  if (not found_it)
  {
    std::cout << "WARNING!! Didn't find " << line_to_find.str() << "  in nutbar_densities.dat " << std::endl;
  }
  return obtd;
}


arma::mat TransitionDensity::ReadTBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, std::string fname)
{
  bool found_it = false;
  int Ji = Jlist_i[J_index_i];
  int Jf = Jlist_f[J_index_f];
  std::ifstream densfile(fname);
  if (not densfile.good())
  {
   std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
   std::cout << "!! TransitionDensity::ReadTBTD -- trouble reading " << fname  << std::endl;
   std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
   exit(EXIT_FAILURE);
  }
  std::string line;
  std::ostringstream line_to_find;
  line_to_find << "Jf nJf  Ji nJi  Lambda = "  << std::setw(3) << std::setprecision(1) << Jf*0.5 << " " << eigvec_f+1
          << "    " << std::setw(3) << std::setprecision(1) << Ji*0.5 << " " << eigvec_i+1
          << "    " << std::setw(3) << std::setprecision(1) << Lambda2*0.5;


  while ( getline( densfile, line) )
  {
     if( line == line_to_find.str() )
     {
       found_it = true;
       break;
     }
  }
  while ( getline( densfile, line) )
  {
     if( line == "-------------- TBTD ---------------------" ) break;
  }
  arma::mat tbtd(ket_J.size(), ket_J.size(), arma::fill::zeros);
  int ibra,iket;
  double tbd;
  while ( line.size() > 2 )
  {
    getline( densfile, line );
    std::istringstream(line) >> ibra >> iket >> tbd;

    int J2ab = ket_J[ibra];
    int J2cd = ket_J[iket];
    tbtd(ibra,iket) = tbd;

    if  ((basename_i==basename_f) and (Ji == Jf) and (eigvec_i==eigvec_f))
    {
      tbtd(iket,ibra) = tbtd(ibra,iket) * (1-std::abs(J2ab-J2cd)%4); 
    }
  }
  densfile.close();
  if (not found_it)
  {
    std::cout << "WARNING!! Didn't find " << line_to_find.str() << "  in nutbar_densities.dat " << std::endl;
  }
  return tbtd;
}

/*
arma::mat TransitionDensity::GetOneBodyTransitionOperator( std::string filename, int& Rank_J, int& Rank_T, int& parity )
{

  std::ifstream opfile(filename);

  if (not opfile.good())
  {
    std::cout << "Trouble reading file " << filename << std::endl;
    Rank_J = -1;
    return arma::mat();
  }

  std::string line;

  while ( line.find("Rank_J") == std::string::npos)   getline(opfile, line);
  std::istringstream( line.substr( line.rfind(":")+1 ) ) >> Rank_J;
  while ( line.find("Rank_T") == std::string::npos)   getline(opfile, line);
  std::istringstream( line.substr( line.rfind(":")+1 ) ) >> Rank_T;
  while ( line.find("Parity") == std::string::npos)   getline(opfile, line);
  std::istringstream( line.substr( line.rfind(":")+1 ) ) >> parity;


  while ( line.find("index") == std::string::npos)   getline(opfile, line);
  std::vector<int> orbits_in;
  getline(opfile, line);
  while ( line.size() > 5 and line.find("a")==std::string::npos )  // no specific reason why 5. Just looking for the empty comment line
  {
    int index,n,l,j2,tz2;
    std::istringstream(line.substr(1)) >> index >> n >> l >> j2 >> tz2;
    // find the corresponding j-orbit
    for (size_t i=0;i<jorbits.size(); ++i)
    { // we use -tz2 because we switch isospin conventions
      if (    m_orbits[jorbits[i]].n==n  and m_orbits[jorbits[i]].l2==2*l and m_orbits[jorbits[i]].j2==j2 and m_orbits[jorbits[i]].tz2==-tz2)
      {
         orbits_in.push_back(i);
      }
    }

    getline(opfile, line);
  }


  arma::mat Op1b(jorbits.size(),jorbits.size(),arma::fill::zeros);

  getline(opfile, line); // skip final header
  int ain,bin,a,b;
  double Op_ab;
  while ( opfile >> ain >> bin >> Op_ab)
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
      std::cout << "WARNING: mat. el. violates triangle contition for rank-" << Rank_J << " operator. -- " << ain << " " << bin << " " << Op_ab << std::endl;
      continue;
    }
    if ( (std::abs(tz2a-tz2b) != 2*Rank_T) and std::abs(Op_ab)>1e-6   )
    {
      std::cout << "WARNING: mat. el. has wrong isospin projection. dTz = " << Rank_T << " operator. -- " << ain << " " << bin << " " << Op_ab << std::endl;
      continue;
    }
    if ( (la+lb)%2 != (1-parity)/2 and std::abs(Op_ab)>1e-6  )
    {
      std::cout << "WARNING: mat. el. has wrong parity (" << parity << ")-- " << ain << " " << bin << " " << Op_ab << std::endl;
      continue;
    }
    Op1b(a,b) = Op_ab;
    Op1b(b,a) = (1 - std::abs(j2a-j2b)%4) * Op_ab; // phase factor (-1)^(ja-jb)
  }

  return Op1b;
}
*/

/*
arma::mat TransitionDensity::GetTwoBodyTransitionOperator( std::string filename , int& Rank_J, int& Rank_T, int& parity)
{

  std::ifstream opfile(filename);
  if (not opfile.good())
  {
    std::cout << "Trouble reading file " << filename << std::endl;
    Rank_J = -2;
    return arma::mat();
  }
  std::string line;

  while ( line.find("Rank_J") == std::string::npos)   getline(opfile, line);
  std::istringstream( line.substr( line.rfind(":")+1 ) ) >> Rank_J;
  while ( line.find("Rank_T") == std::string::npos)   getline(opfile, line);
  std::istringstream( line.substr( line.rfind(":")+1 ) ) >> Rank_T;
  while ( line.find("Parity") == std::string::npos)   getline(opfile, line);
  std::istringstream( line.substr( line.rfind(":")+1 ) ) >> parity;


  // Read in the single particle basis used in the file
  while ( line.find("index") == std::string::npos)   getline(opfile, line);
  std::vector<int> orbits_in;
  getline(opfile, line);
  while ( line.size() > 5 and line.find("a")==std::string::npos )  // no specific reason why 5. Just looking for the empty comment line
  {
    int index,n,l,j2,tz2;
    std::istringstream(line.substr(1)) >> index >> n >> l >> j2 >> tz2;
    // find the corresponding j-orbit
    for (size_t i=0;i<jorbits.size(); ++i)
    {// we use -tz2 because we switch isospin conventions
      if (    m_orbits[jorbits[i]].n==n  and m_orbits[jorbits[i]].l2==2*l and m_orbits[jorbits[i]].j2==j2 and m_orbits[jorbits[i]].tz2==-tz2)
      {
         orbits_in.push_back(i);
      }
    }
    getline(opfile, line);
  }


  SetupKets();


  arma::mat Op2b(ket_J.size(), ket_J.size(), arma::fill::zeros);

  if ( line.find("Jab")==std::string::npos )
    getline(opfile, line); // skip final header
  int ain,bin,cin,din,a,b,c,d,Jab,Jcd;
  double Op_abcd;
  while ( opfile >> ain >> bin >> cin >> din >> Jab >> Jcd >> Op_abcd )
  {
    a = orbits_in[ain-1]; // fortran indexing...
    b = orbits_in[bin-1];
    c = orbits_in[cin-1]; // fortran indexing...
    d = orbits_in[din-1];
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
      Op_abcd *= -(1-std::abs(m_orbits[jorbits[a]].j2 + m_orbits[jorbits[b]].j2 - Jab)%4);
    }
    if (c>d)
    {
      std::swap(c,d);
      Op_abcd *= -(1-std::abs(m_orbits[jorbits[c]].j2 + m_orbits[jorbits[d]].j2 - Jcd)%4);
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

  return Op2b;
}
*/

/*
void TransitionDensity::GetScalarTransitionOperator( std::string filename, double& Op0b, arma::mat& Op1b, arma::mat& Op2b)
{


  // generate all the two body states that are needed
  SetupKets();

  Op1b.set_size( jorbits.size(), jorbits.size() );
  Op2b.set_size( ket_a.size(), ket_a.size() );
  Op1b.zeros();
  Op2b.zeros();
  
  std::ifstream opfile(filename);
  if (! opfile.good() )
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "!!! ERROR: Trouble reading " << filename << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    return;
  }
  std::string line = "!";

  int dummy;

  while (line.find("!") != std::string::npos)
  {
    getline(opfile,line);
    auto zb_ptr = line.find("Zero body term:");
    if ( zb_ptr != std::string::npos)
       std::istringstream( line.substr(zb_ptr + 16) ) >> Op0b;
  }

  std::istringstream iss( line );
  iss >> dummy;
  for (size_t i=0; i<jorbits.size(); ++i)
  {
     iss >> Op1b(i,i);
     Op1b(i,i) *= sqrt( m_orbits[jorbits[i]].j2 + 1); // convert to a reduced matrix element
  }

//  std::cout << Op1b << std::endl;

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

}
*/


arma::vec TransitionDensity::GetDaggerOperator_ax( std::string filename )
{
  arma::vec Op_ax;

  return Op_ax;
}

arma::mat TransitionDensity::GetDaggerOperator_axaxa( std::string filename )
{
  arma::mat Op_axaxa;

  return Op_axaxa;
}







void TransitionDensity::SetupKets()
{
  if (ket_a.size()>0) return;
  for (size_t a=0;a<jorbits.size();++a)
  {
    int ja = m_orbits[jorbits[a]].j2;
    for (size_t b=a; b<jorbits.size();++b)
    {      
      int jb = m_orbits[jorbits[b]].j2;
      int Jmin = std::abs(ja-jb);
      int Jmax = ja+jb;
      for (int J2=Jmin;J2<=Jmax;J2+=2)
      {
        if (a==b and (J2%4)>0) continue;
        ket_a.push_back(a);
        ket_b.push_back(b);
        ket_J.push_back(J2);
      }
    }
  }
}





void TransitionDensity::Jplus(std::vector<key_type>& mvecs_in, std::vector<double>& amp_in, int J2, int M2)
{
  if (M2+2 > J2)
  {
     mvecs_in.clear();
     amp_in.clear();
     return;
  }
  std::unordered_map<key_type,double> amps_out;

  for (size_t i=0;i<mvecs_in.size();++i)
  {
    auto& mvec_in = mvecs_in[i];
    for (size_t i_m=0; i_m<m_orbits.size();++i_m)
    {
      if ( (not mvec_in[i_m]) or (mvec_in[i_m+1])  ) continue; // Pauli principle
      
      int j2  = m_orbits[i_m].j2;
      int mj2 = m_orbits[i_m].mj2;
      if (mj2==j2) continue;
      key_type temp_mvec_out = mvec_in;
      temp_mvec_out.set(i_m,0).set(i_m+1,1);
      amps_out[temp_mvec_out] += sqrt( j2*(j2+2)-mj2*(mj2+2) )*0.5 * amp_in[i];
    }
  }
  
  mvecs_in.clear();
  amp_in.clear();
  for (auto& it_amp : amps_out)
  {
    if ( std::abs(it_amp.second)>1e-6)
    {
      mvecs_in.push_back(it_amp.first);
      amp_in.push_back(it_amp.second);
    }
  }
  // don't forget to normalize...
  double norm = sqrt( inner_product(begin(amp_in),end(amp_in),begin(amp_in),0.0) );
  for (size_t i=0;i<amp_in.size();++i) amp_in[i] /= norm;
}





// Write out the eigenstd::vectors in the Darmstadt MBPT/NCSM format
// so it can be read in by Petr Navratil's TRDENS code
/*
void TransitionDensity::WriteEGV(std::string fname)
{
  std::vector<MschemeOrbit> mscheme_orbits;

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
        if ( sumN<=Zcore )        mscheme_orbits.emplace_back( MschemeOrbit(n,l2,j2,mj2,+1 ) );
        if ( sumN<=Acore-Zcore )  mscheme_orbits.emplace_back( MschemeOrbit(n,l2,j2,mj2,-1 ) );
       }
      }
    }
  }


  // loop over the valence orbits
  int n_core_orbits = mscheme_orbits.size();
//  for (auto& morbit : jbasis_list_i[0].nubasis_a.m_orbits)
//  for (auto& morbit : jbasis_list_i[0].nubasis.m_orbits)
//  for (auto& morbit : jbasis_list_i[0].m_orbits)
  for (auto& morbit : jbasis_i.m_orbits)
  {
     mscheme_orbits.push_back ( morbit );
  }
  // calculate parity of 0hw configurations
  int N=0;
  for (N=0;(N+1)*(N+2)*(N+3)/3 < Z_i;++N){};
  int parity_protons = (N%2);
  for (N=0;(N+1)*(N+2)*(N+3)/3 < A_i-Z_i;++N){};
  int parity_neutrons = (N%2);
  int parity = 1-2*(( parity_protons*(Z_i-Zcore) + parity_neutrons*(A_i-Acore-(Z_i-Zcore)))%2);
  
  std::ofstream output(fname);
  
  std::string interaction_id = "imsrg";
  int Nmax = 0;
  float hw = 20; // this should be irrelevant

  output << std::left << std::setw(10) << Z_i                            << " !  number of protons" << std::endl;
  output << std::left << std::setw(10) << A_i-Z_i                        << " !  number of neutrons" << std::endl;
  output << std::left << std::setw(std::max(10,(int)interaction_id.size()+3)) << interaction_id << " ! interaction id" << std::endl;
  output << std::left << std::setw(10) << hw                             << " ! hbar omega" << std::endl;
  output << std::left << std::setw(10) << Nshell                         << " ! N shells " << std::endl;
  output << std::left << std::setw(10) << mscheme_orbits.size()          << " ! m-scheme orbits " << std::endl;
  output << std::left << std::setw(10) << Nmax                           << " ! Nmax " << std::endl;
  output << std::left << std::setw(10) << amplitudes_i.size()            << " ! number of basis Slater determinants" << std::endl;
  output << std::left << std::setw(10) << std::showpos << parity << std::noshowpos << " ! parity " << std::endl;
  output << std::left << std::setw(10) << MJtot_i                        << " ! total MJ*2" << std::endl;
  output << std::left << std::setw(10) << total_number_levels_i          << " ! Number of eigenstates" << std::endl;
  
  // write out energies, J, T of eigenstates
  for (auto& nuvec : nuvec_list_i)
  {
   int imax = nuvec.no_level;
   if ( max_states_per_J_i.find(nuvec.J2) != max_states_per_J_i.end() ) imax = std::min(imax,max_states_per_J_i[nuvec.J2]);
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
  for (auto it_amp : amplitudes_i)
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
*/

/*

void TransitionDensity::WriteTRDENS_input(std::string fname)
{
  std::ofstream outfile(fname);

  outfile << "T        ! specify which states to take?" << std::endl;
  outfile << "1   " << total_number_levels_i << "   ! ki,nki" << std::endl;
  for (int i=1;i<=total_number_levels_i;++i) outfile << i << " ";
  outfile << std::endl;
  outfile << "1   " << total_number_levels_i << "   ! kf,nkf" << std::endl;
  for (int i=1;i<=total_number_levels_i;++i) outfile << i << " ";
  outfile << std::endl;
  outfile << *std::max_element(begin(Jlist_i),end(Jlist_i)) << "      ! jtotal2max" << std::endl;
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

*/


void TransitionDensity::GetAZFromFileName(  )
{
//  std::string trimmed_basename_i = (basename_i.find("/")==std::string::npos) ? basename_i : basename_i.substr( basename_i.find_last_of("/")+1 );
//  std::string trimmed_basename_f = (basename_f.find("/")==std::string::npos) ? basename_f : basename_f.substr( basename_f.find_last_of("/")+1 );
//  std::string element_i = trimmed_basename_i.substr( 0,2);
//  std::string element_f = trimmed_basename_f.substr( 0,2);
//
//  auto el_position_i = find( periodic_table.begin(),periodic_table.end(), element_i);
//  auto el_position_f = find( periodic_table.begin(),periodic_table.end(), element_f);
//  if (el_position_i == periodic_table.end())
//  {
//   std::cout << "ERROR! : could not find " << element_i << " in periodic table" << std::endl;
//   return;
//  }
//  if (el_position_f == periodic_table.end())
//  {
//   std::cout << "ERROR! : could not find " << element_i << " in periodic table" << std::endl;
//   return;
//  }
//  Z_i = el_position_i - periodic_table.begin();
//  Z_f = el_position_f - periodic_table.begin();
//  std::istringstream( trimmed_basename_i.substr(2,2) ) >> A_i;
//  std::istringstream( trimmed_basename_f.substr(2,2) ) >> A_f;
//
//
//  // now guess at what the core should be, assuming a full major oscillator shell valence space
//  for (int N=0;(N+1)*(N+2)*(N+3)/3<=Z_i;++N) Zcore = (N+1)*(N+2)*(N+3)/3; 
//  for (int N=0;(N+1)*(N+2)*(N+3)/3<=(A_i-Z_i);++N) Acore = Zcore + (N+1)*(N+2)*(N+3)/3; 
//
//  if ( A_i < Acore ) A_i +=100;
//  if ( A_f < Acore ) A_f +=100;
  
}



void TransitionDensity::ReadSPfile()
{
  std::string sp_file_name = sps_file_name.substr(0,sps_file_name.find_last_of(".")) + ".sp";
  std::ifstream infile(sp_file_name);
  if (not infile.good())
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "TransitionDensity::ReadSPfile -- error reading " << sp_file_name << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    exit(EXIT_FAILURE);
  }
  infile.ignore(256,'\n');
  infile.ignore(256,'\n');
  infile >> Acore >> Zcore;

}


void TransitionDensity::SetDensFile( std::string fname )
{
  densfile_name = fname;
  std::ofstream densout(densfile_name);
  densout << "# One and two body transition densities" << std::endl << "#" << std::endl;
  densout << "# One-body basis:" << std::endl;
  densout << "# i   n   l   2j  2tz (proton=+1)" << std::endl;
  for (size_t i=0;i<jorbits.size();++i)
  {
    auto morb = m_orbits[jorbits[i]];
    densout << std::setw(3) << i << " " << std::setw(3) << morb.n << " " << std::setw(3) << morb.l2/2 << " "
            << std::setw(3) << morb.j2 << " " << std::setw(3) << morb.tz2 << std::endl;
  }

  densout << "# Two-body basis: " << std::endl;
  densout << "# i   a   b   J" << std::endl;
  for (size_t i=0;i<ket_a.size();++i)
  {
    densout << std::setw(3) << i << " " << std::setw(3) << ket_a[i] << " " << std::setw(3) << ket_b[i] << " " << std::setw(3) << ket_J[i]/2 << std::endl;
  }

  densout.close(); 

}


std::vector<float> operator*(const float lhs, const std::vector<float>& rhs)
{
  std::vector<float> vout = rhs;
  for (size_t i=0;i<vout.size();++i) vout[i] *= lhs;
  return vout;
}

std::vector<float> operator*(const std::vector<float>& lhs, const std::vector<float>& rhs)
{
  std::vector<float> vout = rhs;
  for (size_t i=0;i<vout.size();++i) vout[i] *= lhs[i];
  return vout;
}

std::vector<float>& operator+=(std::vector<float>& lhs, const std::vector<float>& rhs)
{
  for (size_t i=0;i<lhs.size();++i) lhs[i] += rhs[i];
  return lhs;
}







