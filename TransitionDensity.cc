#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <armadillo>
#include <unordered_map>
#include <omp.h>

#include "TransitionDensity.hh"
#include "JMState.hh"

#define SQRT2 1.4142135623730950488

//#define VERBOSE true

using namespace std;

// a and b are presumably skipped here because a and b are used to label proton/neutron
vector<char> TransitionDensity::an_code = {  '0','1','2','3','4','5','6','7','8','9',
                                             '_','-','c','d','e','f','g','h','i','j',  
                                             'k','l','m','n','o','p','q','r','s','t',  
                                             'u','v','w','x','y','z','A','B','C','D',  
                                             'E','F','G','H','I','J','K','L','M','N',  
                                             'O','P','Q','R','S','T','U','V','W','X',  
                                             '0','1','2','3','4','5','6','7','8','9',  
                                             '_','-','c'};

vector<string> TransitionDensity::periodic_table = {
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
: total_number_levels_i(0),total_number_levels_f(0),densfile_name("none"),Acore(0),Zcore(0)
{}

TransitionDensity::TransitionDensity(vector<int> jlist)
:total_number_levels_i(0),total_number_levels_f(0),
 Jlist_i(jlist), Jlist_f(jlist),densfile_name("none"),Acore(0),Zcore(0)
{}

TransitionDensity::TransitionDensity(vector<int> jlist_i, vector<int> jlist_f)
:total_number_levels_i(0),total_number_levels_f(0),
 Jlist_i(jlist_i), Jlist_f(jlist_f),densfile_name("none"),Acore(0),Zcore(0)
{}








void TransitionDensity::ReadFiles( )
{

  double t_start;
  ifstream testread; // for checking if files exist
  ostringstream ostr;
  vector<string> Afiles_i,Bfiles_i,Afiles_f,Bfiles_f;
  GetAZFromFileName(); // Make a default guess of the core
  ReadSPfile(); // If the sp file is there, use that to get the core
  cout << "Acore,Zcore = " << Acore << " " << Zcore << endl;

  int nvalence_protons_i = Z_i - Zcore;
  int nvalence_neutrons_i = A_i-Z_i - (Acore-Zcore);
  int nvalence_protons_f = Z_f - Zcore;
  int nvalence_neutrons_f = A_f-Z_f - (Acore-Zcore);


//  MJtot = (*min_element(begin(Jlist),end(Jlist)))%2;
  MJtot_i = (*min_element(begin(Jlist_i),end(Jlist_i)))%2;
  MJtot_f = (*min_element(begin(Jlist_f),end(Jlist_f)))%2;

  for (auto j : Jlist_i )
  {
    if (j%2 != A_i%2)
    {
      cout << "Warning A=" << A_i << " and J*2 = " << j << ".  This isn't good" << endl;
    }
  }
  for (auto j : Jlist_f )
  {
    if (j%2 != A_f%2)
    {
      cout << "Warning A=" << A_f << " and J*2 = " << j << ".  This isn't good" << endl;
    }
  }
  
#ifdef VERBOSE
  cout << "TransitionDensity::ReadFiles -- finding all prj and nba files" << endl;
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
//      ostr << basename_i << "a" << iJ << "0" << an_code[36 + icode/2];
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
//      ostr << basename_i << "b" << iJ << "0" << an_code[36 + icode/2];
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
  cout << "TransitionDensity::ReadFiles -- finding all final prj and nba files" << endl;
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
  cout << "TransitionDensity::ReadFiles -- Starting loop over Jlist_i" << endl;
#endif

  for (int Jtot : Jlist_i )
  {
    #ifdef VERBOSE
      cout << "TransitionDensity::ReadFiles -- adding jbasis with " << sps_file_name << "  (";
      for ( auto a : Afiles_i ) cout << " " << a;
      cout << ")  (";
      for ( auto b : Bfiles_i ) cout << " " << b;
      cout << ")  " << Jtot << " " << MJtot_i << endl;;
    #endif
    jbasis_list_i.emplace_back( JBasis( sps_file_name, Afiles_i, Bfiles_i, Jtot, MJtot_i));
    #ifdef VERBOSE
      cout << "TransitionDensity::ReadFiles -- done with emplace_back" << endl;
    #endif
  
  // Guess the name of the xvc file
    ostr.str("");
    ostr.clear();
    ostr << basename_i << an_code[Jtot/2] << an_code[nvalence_protons_i] << an_code[nvalence_neutrons_i] << ".xvc";
    string vecfile = ostr.str();
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
        cout << "ERROR! I cant figure out what the *.xvc file should be. Exiting." << endl;
        cout << "( " << ostr.str() << " ) didnt work. abfile_base = " << basename_i << endl;
        return ;
      }
    }
    testread.close();
    #ifdef VERBOSE
    cout << "TransitionDensity::ReadFiles -- after testread.close, adding Jtot = " << Jtot << ", vecfile = " << vecfile << endl;
    #endif
    nuvec_list_i.emplace_back( NuVec(Jtot) );
    #ifdef VERBOSE
    cout << "Begin read " << vecfile << endl;
    #endif
    nuvec_list_i.back().ReadFile(vecfile);
    #ifdef VERBOSE
    cout << " done." << endl;
    #endif
  }

  profiler.timer["Read_initial_vectors"] += omp_get_wtime() - t_start;
  t_start = omp_get_wtime();

#ifdef VERBOSE
  cout << "TransitionDensity::ReadFiles -- Starting loop over Jlist_f" << endl;
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

  for (int Jtot : Jlist_f )
  {
    auto ji_iter = find( begin(Jlist_i), end(Jlist_i), Jtot);
    if (initial_final_same and ji_iter != end(Jlist_i) )
    {
      jbasis_list_f.emplace_back( jbasis_list_i[ ji_iter-begin(Jlist_i) ] );
    }
    else
    {
      jbasis_list_f.emplace_back( JBasis( sps_file_name, Afiles_f, Bfiles_f, Jtot, MJtot_f));
    }
  
  // Guess the name of the xvc file
    ostr.str("");
    ostr.clear();
    ostr << basename_f << an_code[Jtot/2] << an_code[nvalence_protons_f] << an_code[nvalence_neutrons_f] << ".xvc";
    string vecfile = ostr.str();
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
        cout << "ERROR! I cant figure out what the *.xvc file should be. Exiting." << endl;
        cout << "( " << ostr.str() << " ) didnt work. abfile_base = " << basename_f << endl;
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

  profiler.timer["Read_final_vectors"] += omp_get_wtime() - t_start;


//  m_orbits = jbasis_list_i[0].nubasis_a.m_orbits;
//  m_orbits = jbasis_list_i[0].nubasis.m_orbits;
  m_orbits = jbasis_list_i[0].m_orbits;
  jorbits.clear();
  for (size_t i=0;i<m_orbits.size();++i )
  {
    if (m_orbits[i].mj2 == -m_orbits[i].j2) jorbits.push_back(i);
  }
  

  Nshell=0;
  for (auto morbit : m_orbits)
  {
    Nshell = max(Nshell, 2*morbit.n + morbit.l2/2);
  }
  Nshell++;

  for (size_t ivec=0; ivec<nuvec_list_i.size(); ++ivec)
  {
   int imax = nuvec_list_i[ivec].no_level;
   if ( max_states_per_J_i.find(nuvec_list_i[ivec].J2) != max_states_per_J_i.end() ) imax = min(imax,max_states_per_J_i[nuvec_list_i[ivec].J2]);
    blank_vector_i.push_back(vector<float>(imax, 0.));
    total_number_levels_i += blank_vector_i.back().size();
  }
  for (size_t ivec=0; ivec<nuvec_list_f.size(); ++ivec)
  {
   int imax = nuvec_list_f[ivec].no_level;
   if ( max_states_per_J_f.find(nuvec_list_f[ivec].J2) != max_states_per_J_f.end() ) imax = min(imax,max_states_per_J_f[nuvec_list_f[ivec].J2]);
    blank_vector_f.push_back(vector<float>(imax, 0.));
    total_number_levels_f += blank_vector_f.back().size();
  }

  cout << "done reading files" << endl;
}

/// Calculate the Mscheme wave functions for each eigenstate.
/// If the final and initial states are the same, then don't calculate it twice.
void TransitionDensity::CalculateMschemeAmplitudes()
{
  cout << "start CalculateMschemeAmplitudes" << endl;
  CalculateMschemeAmplitudes_fi( nuvec_list_i, jbasis_list_i, max_states_per_J_i, blank_vector_i, amplitudes_i);
  bool same_f_i = true;
  if (basename_i != basename_f or Jlist_i.size() != Jlist_f.size()) same_f_i = false;
  if (same_f_i)
  {
    for (size_t iJ=0; iJ<=Jlist_i.size(); ++iJ)
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
     CalculateMschemeAmplitudes_fi( nuvec_list_f, jbasis_list_f, max_states_per_J_f, blank_vector_f, amplitudes_f);
  cout << "done" << endl;
}


// vec   : labels a sub-block of good J,pi
// state : basis state with good J
// level : eigenstate of the Hamiltonian
//
// if we call a level |J> and a state |j> and an m-scheme determinant |m>,
// this function evaluates   |J> = sum_{j,m} |m><m|j><j|J>.
//
void TransitionDensity::CalculateMschemeAmplitudes_fi(vector<NuVec>& nuvec_list, vector<JBasis>& jbasis_list, unordered_map<int,int>& max_states_per_J, vector<vector<float>>& blank_vector, unordered_map< key_type, vector<vector<float>> >& amplitudes)
{
  double t_start = omp_get_wtime();
  int nthreads = omp_get_max_threads();
  cout << "max_states_per_J:  ";
  for (auto m : max_states_per_J) cout << m.first << ", " << m.second << ";  ";
  cout << endl;
  cout << "Jvals :  ";
  for (auto& jb : jbasis_list) cout << jb.J2/2 << " ";
  cout << endl;
  for (size_t ivec=0; ivec<nuvec_list.size(); ++ivec)
  {
   const auto& nuvec = nuvec_list[ivec];
   const auto& jbasis = jbasis_list[ivec];
   int imax = nuvec.no_level;
//  #ifdef VERBOSE
   cout << "ivec = " << ivec << endl;
   cout << "no_state = " << nuvec.no_state << endl;
   cout << "no_level = " << nuvec.no_level << "  max_states_per_J = " << max_states_per_J[nuvec.J2] << endl;
//  #endif
   if ( max_states_per_J.find(nuvec.J2) != max_states_per_J.end() ) imax = min(imax,max_states_per_J[nuvec.J2]);
//   cout << "imax = " << imax << endl;


   if (nuvec.no_state > jbasis.basis_states.size() )
   {
     cout << "ERROR -- TransitionDensity::CalculateMschemeAmplitudes_fi -- nuvec.no_state = " << nuvec.no_state << ",  basis_states.size() = " << jbasis.basis_states.size()
          << "   ivec = " << ivec << "  J = " << jbasis.J2/2 << endl;
     return;
   }
  
   double t_start_inner = omp_get_wtime();
   vector<unordered_map< key_type, vector<float> >> local_amplitudes( nthreads );
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
     const JMState jmst = jbasis.GetBasisState(istate);
//     profiler.timer["GetBasisState"] += omp_get_wtime() - t_getstate;
     for (auto& it_mstate : jmst.m_coefs)
     {
       auto& key = it_mstate.first;
       const float& m_coef = it_mstate.second;

       if( local_amplitudes[thread_num].find(key) == local_amplitudes[thread_num].end() ) local_amplitudes[thread_num][key] = vector<float>(imax,0.);
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

  if ((J2i+Lambda2 < J2f) or (abs(J2i-Lambda2)>J2f)) return 0;
 
//  cout << "  In OBTD " << J2i << " " << J2f << " " << J_index_i << " " << J_index_f << "  " << eigvec_i << " " << eigvec_f << endl;
//  cout << "    with Lambda2 = " << Lambda2 << endl;
 
  int mu = Lambda2%2;
  int Mi = J2i%2;
  int Mf = J2f%2;

  int j2_a = m_orbits[m_index_a].j2;
  int j2_b = m_orbits[m_index_b].j2;
//  cout << "    ja,jb = " << j2_a << " " << j2_b << endl;

//  cout << "    a,b = " << m_index_a << " " << m_index_b << endl;
//  cout <<  "     " << m_orbits[m_index_a].n << " " << m_orbits[m_index_a].l2/2 << " " << j2_a << " " << m_orbits[m_index_a].tz2
//       <<  "     " << m_orbits[m_index_b].n << " " << m_orbits[m_index_b].l2/2 << " " << j2_b << " " << m_orbits[m_index_b].tz2 << endl;

  m_index_a += (j2_a - m_orbits[m_index_a].mj2)/2;
  m_index_b += (j2_b - m_orbits[m_index_b].mj2)/2;

  
  // find m-scheme orbits so that m_a = m_b, which will work for mu=0
  double obd = 0;
  int ma_min = max(-j2_a,mu-j2_b);
  int ma_max = min(j2_a,mu+j2_b);

//  vector<vector<mvec_type>> keys_i;
  vector<key_type> keys_i;
  vector<double> amp_vec_i;
//  cout << "size of amplitudes_i = " << amplitudes_i.size() << endl;
  for (auto& it_amp : amplitudes_i )
  {
     double amp_i = it_amp.second[J_index_i][eigvec_i];
     if (abs(amp_i)<1e-7) continue;
     auto& key = it_amp.first;
     keys_i.push_back( key );
     amp_vec_i.push_back(amp_i);
  }

  double clebsch_fi = CG(J2i,Mi,Lambda2,mu,J2f,Mf);
  if (abs(clebsch_fi)<1e-9)
  {
     clebsch_fi = CG(J2i,Mi+2,Lambda2,mu-2,J2f,Mf);
     if (abs(clebsch_fi)<1e-9)
     {
        cout << " ERROR:    Still got zero Clebsch" << endl;
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

//    uint64_t mask_a = (0x1L<<ia);
//    uint64_t mask_b = (0x1L<<ib);

    // convention: tilded destruction operator b~(m) = (-1)**(jb + mb) b(-m)
    //                                         b(m)  = (-1)**(jb -mb) b~(-m)
    int phase_b = (1-(j2_b-mb)%4);
    double clebsch = CG(j2_a,ma,j2_b,-mb,Lambda2,mu) ;
//    cout << "clebsch: " << j2_a << " " << ma << " " << j2_b << " " << -mb << " " << Lambda2 << " " << mu << endl;

    for (size_t iamp=0; iamp<amp_vec_i.size(); ++iamp)
    {
      auto& key = keys_i[iamp];
      auto& amp_i = amp_vec_i[iamp];

//      if ( not( key[0] & mask_b )) continue;
//      if ( ia != ib and   ( key[0] & mask_a) ) continue;
      if ( not key[ib] ) continue;
      auto new_key = key;
      new_key.set(ib,0);
      if ( new_key[ia] ) continue;
      new_key.set(ia,1);
//      new_key[0] &= ~mask_b;  // remove orbit b
//      new_key[0] |=  mask_a;  // add to orbit a
      if (amplitudes_f.find(new_key) == amplitudes_f.end() ) continue;


      int phase_ladder = 1;
//      for (int iphase=min(ia,ib)+1;iphase<max(ia,ib);++iphase) if( (key[0] >>iphase )&0x1L) phase_ladder *=-1;
      for (int iphase=min(ia,ib)+1;iphase<max(ia,ib);++iphase) if( key[iphase]) phase_ladder *=-1;
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
  int mu = Lambda2%2;
  int Mi = J2i%2;
  int Mf = J2f%2;
  
  int j2_a = m_orbits[m_index_a].j2;
  int j2_b = m_orbits[m_index_b].j2;
  int j2_c = m_orbits[m_index_c].j2;
  int j2_d = m_orbits[m_index_d].j2;

  // start out with the maximally-projected m state
  m_index_a += (j2_a - m_orbits[m_index_a].mj2)/2;
  m_index_b += (j2_b - m_orbits[m_index_b].mj2)/2;
  m_index_c += (j2_c - m_orbits[m_index_c].mj2)/2;
  m_index_d += (j2_d - m_orbits[m_index_d].mj2)/2;

//  bool printout = false;
//  if (j2_a==7 and j2_b==7 and j2_c==7 and j2_d==7 and J2ab==4 and J2cd==4) printout=true;


  // check some triangle conditions
  if ((J2f+Lambda2 < J2i) or (abs(J2f-Lambda2)>J2i)) return 0;
  if ((j2_a+j2_b<J2ab) or (abs(j2_a-j2_b)>J2ab)) return 0;
  if ((j2_c+j2_d<J2cd) or (abs(j2_c-j2_d)>J2cd)) return 0;
  if ((J2ab+J2cd < Lambda2) or (abs(J2ab-J2cd)>Lambda2)) return 0;

  
//  vector<vector<mvec_type>> keys_i;
  vector<key_type> keys_i;
  vector<double> amp_vec_i;
  for (auto& it_amp : amplitudes_i )
  {
     double amp_i = it_amp.second[J_index_i][eigvec_i];
     if (abs(amp_i)<1e-7) continue;
     auto& key = it_amp.first;
     keys_i.push_back( key );
     amp_vec_i.push_back(amp_i);
  }

  double clebsch_fi = CG(J2i,Mi,Lambda2,mu,J2f,Mf);
  if (abs(clebsch_fi)<1e-9)
  {
     clebsch_fi = CG(J2i,Mi+2,Lambda2,mu-2,J2f,Mf);
     if (abs(clebsch_fi)<1e-9)
     {
        cout << " ERROR:    Still got zero Clebsch" << endl;
        cout << "           J2i=" << J2i << " M2i=" << Mi+2 << " Lambda2=" << Lambda2 << " mu=" << mu-2 << " J2f=" << J2f << " M2f=" << Mf << endl;
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
  int Mab_min = max(-J2ab,mu-J2cd);
  int Mab_max = min( J2ab,mu+J2cd);
  for (int Mab=Mab_min;Mab<=Mab_max;Mab+=2)
  {
    int Mcd = mu - Mab;
    int phasecd = (1- abs(J2cd - Mcd)%4); // phase from getting rid of the tildes
    double clebsch_abcd = CG(J2ab,Mab,J2cd,Mcd,Lambda2,mu);
    if ( abs(clebsch_abcd)<1e-7) continue;
    int ma_min = max(-j2_a, Mab-j2_b);
    int ma_max = min(j2_a, Mab+j2_b);
    if (m_index_a==m_index_b) ma_max = min(ma_max, Mab/2);

    for (int ma=ma_min;ma<=ma_max;ma+=2)
    {
      int mb=Mab-ma;
      double clebsch_ab = CG(j2_a,ma, j2_b,mb, J2ab,Mab);
      if ( abs(clebsch_ab)<1e-7) continue;
      int ia = m_index_a - ( j2_a -ma )/2;
      int ib = m_index_b - ( j2_b -mb )/2;
      if (ib==ia) continue;
      int mc_min = max(-j2_c, Mab-j2_d);
      int mc_max = min(j2_c, Mab+j2_d);
      if (m_index_c == m_index_d) mc_max = min(j2_c, Mab/2);

      for (int mc=mc_min;mc<=mc_max;mc+=2)
      {
        int md = -Mcd-mc;
        double clebsch_cd = CG(j2_c,mc, j2_d, md, J2cd,-Mcd); // Mcd = -Mab, and another (-) comes from getting rid of the tildes
        if ( abs(clebsch_cd)<1e-7) continue;
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
          for (int iphase = min(ic,id)+1;iphase<max(ic,id);++iphase)  if( key[iphase]) phase_ladder *=-1;
          for (int iphase = min(ia,ib)+1;iphase<max(ia,ib);++iphase)  if( new_key[iphase]) phase_ladder *=-1;

          tbd += amp_i * amp_f * clebsch_abcd * clebsch_ab * clebsch_cd * phasecd * phase_ladder;
        }
      }
    }
  }

  tbd *= sqrt((J2f+1.)/(Lambda2+1.)) / clebsch_fi * norm;
  return tbd;

}



arma::mat TransitionDensity::CalcOBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2)
{

  double t_start = omp_get_wtime();

  int Ji = Jlist_i[J_index_i];
  int Jf = Jlist_f[J_index_f];
  size_t njorb = jorbits.size();
  arma::mat obtd(njorb,njorb,arma::fill::zeros);
  ofstream densout( densfile_name, ios::app );
  if ( densfile_name != "none")
  {
     densout << endl;
     densout << "Jf nJf  Ji nJi  Lambda = " << setw(3) << setprecision(1) << Jf*0.5 << " " << eigvec_f+1
             << "    " << setw(3) << setprecision(1) << Ji*0.5 << " " << eigvec_i+1
             << "    " << setw(3) << setprecision(1) << Lambda2*0.5  << endl;
     densout << "-------------- OBTD ---------------------" << endl;
  }

  for (size_t i=0; i<njorb; ++i)
  {
    int j2i = m_orbits[jorbits[i]].j2;
    int jmin = 0;
    if (basename_i==basename_f and Ji==Jf and eigvec_i==eigvec_f) jmin = i;
    for (size_t j=jmin; j<njorb; ++j)
    {
      obtd(i,j) = OBTD( J_index_i, eigvec_i, J_index_f, eigvec_f, jorbits[i], jorbits[j], Lambda2);
      if (Lambda2 == 0)  obtd(i,j) /= sqrt( Ji+1 );

      if (densfile_name != "none" and abs(obtd(i,j))>1e-7)
      {
         densout << setw(3) << i << " " << setw(3) << j << " "  << setw(12) << fixed << setprecision(8) << obtd(i,j)  << endl;
      }
      if (basename_i==basename_f and Ji==Jf and eigvec_i==eigvec_f)
      {
        int j2j = m_orbits[jorbits[j]].j2;
        obtd(j,i) = (1-abs(j2j-j2i)%4) * obtd(i,j);
      }
    }
  }

  densout.close();
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
  if (abs(Ji-Jf)>Lambda2 or Ji+Jf<Lambda2) return tbtd;


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
        tbtd(iket,ibra) = tbtd(ibra,iket) * (1-abs(J2ab-J2cd)%4); 
      }
    }
  }

  if (Lambda2==0) tbtd /= sqrt( Ji+1.);
  if ( densfile_name != "none" )
  {
    ofstream densout(densfile_name, ios::app);
    densout << endl;
    densout << "-------------- TBTD ---------------------" << endl;
    for (size_t i=0;i<tbtd.n_rows;++i)
    {
      for (size_t j=0;j<tbtd.n_cols;++j)
      {
         if (abs(tbtd(i,j))>1e-7)
         densout << setw(3) << i << " " << setw(3) << j << " "  << setw(12) << fixed << setprecision(8) << tbtd(i,j) << endl;
      }
    }
    densout.close();
  }

  profiler.timer["CalcTBTD"] += omp_get_wtime() - t_start;
  return tbtd;
}



arma::mat TransitionDensity::ReadOBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, string fname)
{
  bool found_it = false;
  int Ji = Jlist_i[J_index_i];
  int Jf = Jlist_f[J_index_f];
  size_t njorb = jorbits.size();
  ifstream densfile(fname);
  string line;
  ostringstream line_to_find;
  line_to_find <<  "Jf nJf  Ji nJi  Lambda = "   << setw(3) << setprecision(1) << Jf*0.5 << " " << eigvec_f+1
          << "    " << setw(3) << setprecision(1) << Ji*0.5 << " " << eigvec_i+1
          << "    " << setw(3) << setprecision(1) << Lambda2*0.5;


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
    istringstream(line) >> i >> j >> obd;
    obtd(i,j) = obd;
    if (basename_i==basename_f and Ji==Jf and eigvec_i==eigvec_f)
    {
      int j2i = m_orbits[jorbits[i]].j2;
      int j2j = m_orbits[jorbits[j]].j2;
      obtd(j,i) = (1-abs(j2j-j2i)%4) * obtd(i,j);
    }

  }
  densfile.close();
  if (not found_it)
  {
    cout << "WARNING!! Didn't find " << line_to_find.str() << "  in nutbar_densities.dat " << endl;
  }
  return obtd;
}


arma::mat TransitionDensity::ReadTBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, string fname)
{
  bool found_it = false;
  int Ji = Jlist_i[J_index_i];
  int Jf = Jlist_f[J_index_f];
  ifstream densfile(fname);
  string line;
  ostringstream line_to_find;
  line_to_find << "Jf nJf  Ji nJi  Lambda = "  << setw(3) << setprecision(1) << Jf*0.5 << " " << eigvec_f+1
          << "    " << setw(3) << setprecision(1) << Ji*0.5 << " " << eigvec_i+1
          << "    " << setw(3) << setprecision(1) << Lambda2*0.5;


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
    istringstream(line) >> ibra >> iket >> tbd;

    int J2ab = ket_J[ibra];
    int J2cd = ket_J[iket];
    tbtd(ibra,iket) = tbd;

    if  ((basename_i==basename_f) and (Ji == Jf) and (eigvec_i==eigvec_f))
    {
      tbtd(iket,ibra) = tbtd(ibra,iket) * (1-abs(J2ab-J2cd)%4); 
    }
  }
  densfile.close();
  if (not found_it)
  {
    cout << "WARNING!! Didn't find " << line_to_find.str() << "  in nutbar_densities.dat " << endl;
  }
  return tbtd;
}


arma::mat TransitionDensity::GetOneBodyTransitionOperator( string filename, int& Rank_J, int& Rank_T, int& parity )
{

  ifstream opfile(filename);

  if (not opfile.good())
  {
    cout << "Trouble reading file " << filename << endl;
    Rank_J = -1;
    return arma::mat();
  }

  string line;

  while ( line.find("Rank_J") == string::npos)   getline(opfile, line);
  istringstream( line.substr( line.rfind(":")+1 ) ) >> Rank_J;
  while ( line.find("Rank_T") == string::npos)   getline(opfile, line);
  istringstream( line.substr( line.rfind(":")+1 ) ) >> Rank_T;
  while ( line.find("Parity") == string::npos)   getline(opfile, line);
  istringstream( line.substr( line.rfind(":")+1 ) ) >> parity;


  while ( line.find("index") == string::npos)   getline(opfile, line);
  vector<int> orbits_in;
  getline(opfile, line);
  while ( line.size() > 5 and line.find("a")==string::npos )  // no specific reason why 5. Just looking for the empty comment line
  {
    int index,n,l,j2,tz2;
    istringstream(line.substr(1)) >> index >> n >> l >> j2 >> tz2;
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
    if ( (abs(j2a-j2b) > 2*Rank_J) or (j2a+j2b < 2*Rank_J) and abs(Op_ab)>1e-6 )
    {
      std::cout << "WARNING: mat. el. violates triangle contition for rank-" << Rank_J << " operator. -- " << ain << " " << bin << " " << Op_ab << std::endl;
      continue;
    }
    if ( (abs(tz2a-tz2b) != 2*Rank_T) and abs(Op_ab)>1e-6   )
    {
      std::cout << "WARNING: mat. el. has wrong isospin projection. dTz = " << Rank_T << " operator. -- " << ain << " " << bin << " " << Op_ab << std::endl;
      continue;
    }
    if ( (la+lb)%2 != (1-parity)/2 and abs(Op_ab)>1e-6  )
    {
      std::cout << "WARNING: mat. el. has wrong parity (" << parity << ")-- " << ain << " " << bin << " " << Op_ab << std::endl;
      continue;
    }
    Op1b(a,b) = Op_ab;
    Op1b(b,a) = (1 - abs(j2a-j2b)%4) * Op_ab; // phase factor (-1)^(ja-jb)
  }

  return Op1b;
}



arma::mat TransitionDensity::GetTwoBodyTransitionOperator( string filename , int& Rank_J, int& Rank_T, int& parity)
{

  ifstream opfile(filename);
  if (not opfile.good())
  {
    cout << "Trouble reading file " << filename << endl;
    Rank_J = -2;
    return arma::mat();
  }
  string line;

  while ( line.find("Rank_J") == string::npos)   getline(opfile, line);
  istringstream( line.substr( line.rfind(":")+1 ) ) >> Rank_J;
  while ( line.find("Rank_T") == string::npos)   getline(opfile, line);
  istringstream( line.substr( line.rfind(":")+1 ) ) >> Rank_T;
  while ( line.find("Parity") == string::npos)   getline(opfile, line);
  istringstream( line.substr( line.rfind(":")+1 ) ) >> parity;


  // Read in the single particle basis used in the file
  while ( line.find("index") == string::npos)   getline(opfile, line);
  vector<int> orbits_in;
  getline(opfile, line);
  while ( line.size() > 5 and line.find("a")==string::npos )  // no specific reason why 5. Just looking for the empty comment line
  {
    int index,n,l,j2,tz2;
    istringstream(line.substr(1)) >> index >> n >> l >> j2 >> tz2;
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

  if ( line.find("Jab")==string::npos )
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
    if ( (abs(Jab-Jcd) > Rank_J) or (Jab+Jcd < Rank_J) and abs(Op_abcd)>1e-6 )
    {
      std::cout << "WARNING: mat. el. violates triangle contition for rank-" << Rank_J << " operator. -- "
               << ain << " " << bin << " " << cin << " " << din << " " << Jab << " " << Jcd << " " << Op_abcd << std::endl;
      continue;
    }
    if ( (abs(tz2a+tz2b-tz2c-tz2d) != 2*Rank_T)  and abs(Op_abcd)>1e-6)
    {
      std::cout << "WARNING: mat. el. has wrong isospin projection. dTz = " << Rank_T << " operator. -- " 
               << ain << " " << bin << " " << cin << " " << din << " " << Jab << " " << Jcd << " " << Op_abcd << std::endl;
      continue;
    }
    if ( (la+lb+lc+ld)%2 != (1-parity)/2 and abs(Op_abcd)>1e-6)
    {
      std::cout << "WARNING: mat. el. has wrong parity -- " 
               << ain << " " << bin << " " << cin << " " << din << " " << Jab << " " << Jcd << " " << Op_abcd << std::endl;
      continue;
    }
    Jab *=2;
    Jcd *=2;
    if (a>b)
    {
      swap(a,b);
      Op_abcd *= -(1-abs(m_orbits[jorbits[a]].j2 + m_orbits[jorbits[b]].j2 - Jab)%4);
    }
    if (c>d)
    {
      swap(c,d);
      Op_abcd *= -(1-abs(m_orbits[jorbits[c]].j2 + m_orbits[jorbits[d]].j2 - Jcd)%4);
    }

    size_t ibra=0,iket=0;
    while( ibra<ket_a.size() and not( (ket_a[ibra]==a) and (ket_b[ibra]==b) and ket_J[ibra]==Jab) ) ibra++;
    while( iket<ket_a.size() and not( (ket_a[iket]==c) and (ket_b[iket]==d) and ket_J[iket]==Jcd) ) iket++;

    if (ibra >= ket_a.size() or iket >= ket_a.size())
    {
      cout << "trouble:  " << a << " " << b << " " << c << " " << d << " " <<Jab << " " << Jcd << " " << Op_abcd << endl;
      cout << "   ibra = " << ibra << "  iket = " << iket << "  size = " << ket_a.size() << endl;
    }

    Op2b(ibra,iket) = Op_abcd;
    if (ibra!=iket)
      Op2b(iket,ibra) = (1 - abs(Jab+Jcd )%4) * Op_abcd; // phase factor (-1)^(Jab-Jcd)
  }

  return Op2b;
}



//void TransitionDensity::GetScalarTransitionOperator( string filename, arma::mat& Op1b, arma::mat& Op2b)
void TransitionDensity::GetScalarTransitionOperator( string filename, double& Op0b, arma::mat& Op1b, arma::mat& Op2b)
{


  // generate all the two body states that are needed
  SetupKets();

  Op1b.set_size( jorbits.size(), jorbits.size() );
  Op2b.set_size( ket_a.size(), ket_a.size() );
  Op1b.zeros();
  Op2b.zeros();
  
  ifstream opfile(filename);
  if (! opfile.good() )
  {
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "!!! ERROR: Trouble reading " << filename << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    return;
  }
  string line = "!";

  int dummy;

  while (line.find("!") != string::npos)
  {
    getline(opfile,line);
    auto zb_ptr = line.find("Zero body term:");
    if ( zb_ptr != string::npos)
       istringstream( line.substr(zb_ptr + 16) ) >> Op0b;
  }

  istringstream iss( line );
  iss >> dummy;
  for (size_t i=0; i<jorbits.size(); ++i)
  {
     iss >> Op1b(i,i);
     Op1b(i,i) *= sqrt( m_orbits[jorbits[i]].j2 + 1); // convert to a reduced matrix element
  }

//  cout << Op1b << endl;

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
      ME *= -(1-abs(oa.j2 + ob.j2 -J*2 )%4);
      swap(a,b);
    }
    if (c>d)
    {
      ME *= -(1-abs(oc.j2 + od.j2 -J*2 )%4);
      swap(c,d);
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


void TransitionDensity::SetupKets()
{
  if (ket_a.size()>0) return;
  for (size_t a=0;a<jorbits.size();++a)
  {
    int ja = m_orbits[jorbits[a]].j2;
    for (size_t b=a; b<jorbits.size();++b)
    {      
      int jb = m_orbits[jorbits[b]].j2;
      int Jmin = abs(ja-jb);
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





void TransitionDensity::Jplus(vector<key_type>& mvecs_in, vector<double>& amp_in, int J2, int M2)
{
  if (M2+2 > J2)
  {
     mvecs_in.clear();
     amp_in.clear();
     return;
  }
  unordered_map<key_type,double> amps_out;

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
    if ( abs(it_amp.second)>1e-6)
    {
      mvecs_in.push_back(it_amp.first);
      amp_in.push_back(it_amp.second);
    }
  }
  // don't forget to normalize...
  double norm = sqrt( inner_product(begin(amp_in),end(amp_in),begin(amp_in),0.0) );
  for (size_t i=0;i<amp_in.size();++i) amp_in[i] /= norm;
}





// Write out the eigenvectors in the Darmstadt MBPT/NCSM format
// so it can be read in by Petr Navratil's TRDENS code
void TransitionDensity::WriteEGV(string fname)
{
  vector<MschemeOrbit> mscheme_orbits;

  // find all the core orbits
  for (int N=0;N<Nshell;++N) 
  {
    int sumN = (N+1)*(N+2)*(N+3)/3;
    for (int n=0;2*n<=N;++n)
    {
      int l2 = 2*(N-2*n);
      for (int j2=l2+1;j2>=max(0,l2-1);j2-=2)
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
  for (auto& morbit : jbasis_list_i[0].m_orbits)
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
  
  ofstream output(fname);
  
  string interaction_id = "imsrg";
  int Nmax = 0;
  float hw = 20; // this should be irrelevant

  output << left << setw(10) << Z_i                            << " !  number of protons" << endl;
  output << left << setw(10) << A_i-Z_i                        << " !  number of neutrons" << endl;
  output << left << setw(max(10,(int)interaction_id.size()+3)) << interaction_id << " ! interaction id" << endl;
  output << left << setw(10) << hw                             << " ! hbar omega" << endl;
  output << left << setw(10) << Nshell                         << " ! N shells " << endl;
  output << left << setw(10) << mscheme_orbits.size()          << " ! m-scheme orbits " << endl;
  output << left << setw(10) << Nmax                           << " ! Nmax " << endl;
  output << left << setw(10) << amplitudes_i.size()            << " ! number of basis Slater determinants" << endl;
  output << left << setw(10) << showpos << parity << noshowpos << " ! parity " << endl;
  output << left << setw(10) << MJtot_i                        << " ! total MJ*2" << endl;
  output << left << setw(10) << total_number_levels_i          << " ! Number of eigenstates" << endl;
  
  // write out energies, J, T of eigenstates
  for (auto& nuvec : nuvec_list_i)
  {
   int imax = nuvec.no_level;
   if ( max_states_per_J_i.find(nuvec.J2) != max_states_per_J_i.end() ) imax = min(imax,max_states_per_J_i[nuvec.J2]);
   for (int ilevel=0;ilevel<imax;++ilevel)
   {
     output << right << fixed << setw(12) << setprecision(4) << nuvec.alpha[ilevel] << " " << nuvec.J2/2.0 << " " << 0 << " " << 0.0 << endl;
   }
  }
  
  
  output << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  output << "!!!  now list mscheme single-particle basis !!!" << endl;
  output << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  // write single-particle basis in mscheme
  
  
  
  int s = 1;
  for (auto& morbit : mscheme_orbits)
  {
    output << setw(4) << s << " " << setw(3) << morbit.n << " " << setw(3) << morbit.l2 << " " << setw(3) << morbit.j2 << " " << setw(3) << morbit.mj2 << " " << setw(3) << morbit.tz2 << endl;;
    s++;
  }
  
  output << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  output << "!!!  now list mscheme basis states and amplitudes !!!" << endl;
  output << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  
  
  
  // Write m-scheme basis occupations and eigenvector coefficients
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
        output << scientific << setw(13) << amp << " ";
      }
    }
    output << endl;
  
  }

}




void TransitionDensity::WriteTRDENS_input(string fname)
{
  ofstream outfile(fname);

  outfile << "T        ! specify which states to take?" << endl;
  outfile << "1   " << total_number_levels_i << "   ! ki,nki" << endl;
  for (int i=1;i<=total_number_levels_i;++i) outfile << i << " ";
  outfile << endl;
  outfile << "1   " << total_number_levels_i << "   ! kf,nkf" << endl;
  for (int i=1;i<=total_number_levels_i;++i) outfile << i << " ";
  outfile << endl;
  outfile << *max_element(begin(Jlist_i),end(Jlist_i)) << "      ! jtotal2max" << endl;
  outfile << "0                   ! irestart" << endl;
  outfile << "2                   ! majortot" << endl;
  outfile << "F                   ! irem_cal" << endl;
  outfile << "F                   ! formfcal" << endl;
  outfile << "F                   ! radial  " << endl;
  outfile << "F                   ! momdist " << endl;
  outfile << "T                   ! twobdcal" << endl;
  outfile << "1                   ! ipn     " << endl;
  outfile << "IMSRG.int_iso                 " << endl;
  outfile << "IMSRG_E2_1b.op_iso            " << endl;
  outfile << "IMSRG_E2_2b.op_iso            " << endl;
  outfile << "F                   ! cluster " << endl;
  outfile << "F                   ! antoine " << endl;
  outfile << "24.d0               ! hbar*Omega for antoine file" << endl;
  outfile << "1                   ! #init states in antoine file" << endl;
  outfile << "0  0 " << endl;
  outfile << "1                   ! #fin states in antoine file" << endl;
  outfile << "0  0 " << endl;
  outfile << "T                   ! mbpt_ncsm " << endl;
  outfile << "F                   ! redstick" << endl;
  outfile << "F                   ! mfd_james" << endl;
  outfile << "F                   ! threebdcal" << endl;
  outfile << "F                   ! NCSMC_kernels" << endl;
  outfile << "3                   ! num_of_interaction_files" << endl;
  outfile << "../vrelnp_rgm.int_H2srg-n3lo2.2_2014 " << endl;
  outfile << "../vrelpp_rgm.int_H2srg-n3lo2.2_2014 " << endl;
  outfile << "../vrelnn_rgm.int_H2srg-n3lo2.2_2014 " << endl;
  outfile << "F                   ! V3Nint " << endl;
  outfile << "14 14 14            ! N1_max,N12_max,N123_max " << endl;
  outfile << "chi2b3b400cD-02cE0098_srg0800ho40C_eMax14_EMax14_hwHO016.me3j_bin" << endl;


  outfile.close();


}




void TransitionDensity::GetAZFromFileName(  )
{
  string trimmed_basename_i = (basename_i.find("/")==string::npos) ? basename_i : basename_i.substr( basename_i.find_last_of("/")+1 );
  string trimmed_basename_f = (basename_f.find("/")==string::npos) ? basename_f : basename_f.substr( basename_f.find_last_of("/")+1 );
  string element_i = trimmed_basename_i.substr( 0,2);
  string element_f = trimmed_basename_f.substr( 0,2);

  auto el_position_i = find( periodic_table.begin(),periodic_table.end(), element_i);
  auto el_position_f = find( periodic_table.begin(),periodic_table.end(), element_f);
  if (el_position_i == periodic_table.end())
  {
   cout << "ERROR! : could not find " << element_i << " in periodic table" << endl;
   return;
  }
  if (el_position_f == periodic_table.end())
  {
   cout << "ERROR! : could not find " << element_i << " in periodic table" << endl;
   return;
  }
  Z_i = el_position_i - periodic_table.begin();
  Z_f = el_position_f - periodic_table.begin();
  istringstream( trimmed_basename_i.substr(2,2) ) >> A_i;
  istringstream( trimmed_basename_f.substr(2,2) ) >> A_f;


  // now guess at what the core should be, assuming a full major oscillator shell valence space
  for (int N=0;(N+1)*(N+2)*(N+3)/3<=Z_i;++N) Zcore = (N+1)*(N+2)*(N+3)/3; 
  for (int N=0;(N+1)*(N+2)*(N+3)/3<=(A_i-Z_i);++N) Acore = Zcore + (N+1)*(N+2)*(N+3)/3; 

  if ( A_i < Acore ) A_i +=100;
  if ( A_f < Acore ) A_f +=100;
  
}



void TransitionDensity::ReadSPfile()
{
  string sp_file_name = sps_file_name.substr(0,sps_file_name.find_last_of(".")) + ".sp";
  ifstream infile(sp_file_name);
  if (not infile.good()) return;
  infile.ignore(256,'\n');
  infile.ignore(256,'\n');
  infile >> Acore >> Zcore;

}


void TransitionDensity::SetDensFile( string fname )
{
  densfile_name = fname;
  ofstream densout(densfile_name);
  densout << "# One and two body transition densities" << endl << "#" << endl;
  densout << "# One-body basis:" << endl;
  densout << "# i   n   l   2j  2tz (proton=+1)" << endl;
  for (size_t i=0;i<jorbits.size();++i)
  {
    auto morb = m_orbits[jorbits[i]];
    densout << setw(3) << i << " " << setw(3) << morb.n << " " << setw(3) << morb.l2/2 << " "
            << setw(3) << morb.j2 << " " << setw(3) << morb.tz2 << endl;
  }

  densout << "# Two-body basis: " << endl;
  densout << "# i   a   b   J" << endl;
  for (size_t i=0;i<ket_a.size();++i)
  {
    densout << setw(3) << i << " " << setw(3) << ket_a[i] << " " << setw(3) << ket_b[i] << " " << setw(3) << ket_J[i]/2 << endl;
  }

  densout.close(); 

}


vector<float> operator*(const float lhs, const vector<float>& rhs)
{
  vector<float> vout = rhs;
  for (size_t i=0;i<vout.size();++i) vout[i] *= lhs;
  return vout;
}

vector<float> operator*(const vector<float>& lhs, const vector<float>& rhs)
{
  vector<float> vout = rhs;
  for (size_t i=0;i<vout.size();++i) vout[i] *= lhs[i];
  return vout;
}

vector<float>& operator+=(vector<float>& lhs, const vector<float>& rhs)
{
  for (size_t i=0;i<lhs.size();++i) lhs[i] += rhs[i];
  return lhs;
}







