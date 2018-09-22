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




TransitionDensity::TransitionDensity()
{}

TransitionDensity::TransitionDensity(std::vector<int> jlist)
:
 Jlist_i(jlist), Jlist_f(jlist)
{}

TransitionDensity::TransitionDensity(std::vector<int> jlist_i, std::vector<int> jlist_f)
:
 Jlist_i(jlist_i), Jlist_f(jlist_f)
{}








/// Calculate the Mscheme wave functions for each eigenstate.
/// If the final and initial states are the same, then don't calculate it twice.
//void TransitionDensity::CalculateMschemeAmplitudes()
void TransitionDensity::CalculateMschemeAmplitudes(Settings& settings)
{
  CalculateMschemeAmplitudes_fi( nuvec_list_i, settings.jbasis_i,  blank_vector_i, amplitudes_i);

  if (settings.same_states_fi)
  {
    amplitudes_f = amplitudes_i;
  }
  else
  {
    CalculateMschemeAmplitudes_fi( nuvec_list_f, settings.jbasis_f,  blank_vector_f, amplitudes_f);
  }
}


// vec   : labels a sub-block of good J,pi
// state : basis state with good J
// level : eigenstate of the Hamiltonian
//
// if we call a level |J> and a state |j> and an m-scheme determinant |m>,
// this function evaluates   |J> = sum_{j,m} |m><m|j><j|J>.
//
void TransitionDensity::CalculateMschemeAmplitudes_fi(std::vector<NuVec>& nuvec_list, JBasis& jbasis, std::vector<std::vector<float>>& blank_vector, std::unordered_map< key_type, std::vector<std::vector<float>> >& amplitudes)
{
  double t_start = omp_get_wtime();
  int nthreads = omp_get_max_threads();

  // Choose the MJ2 value equal to the minimum J in the basis. I do this because m-scheme dimension grows with decreasing M.
  int MJ2 = (*std::min_element(begin(jbasis.basis_states),end(jbasis.basis_states),jbasis.basis_states.value_comp())).first % 2 ;
  for (size_t ivec=0; ivec<nuvec_list.size(); ++ivec)
  {
   const auto& nuvec = nuvec_list[ivec];
   int J2 = nuvec.J2;
   int imax = nuvec.no_level;
  #ifdef VERBOSE
   std::cout << "Naive m-scheme dimension = " << jbasis.GetNaiveMschemeDimension(J2) << "  ";
   std::cout << "imax = " << imax << std::endl;
  #endif


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
     const JMState jmst = jbasis.GetBasisState(istate,J2,MJ2);  // can we be smarter about picking an initial projection?
     for (auto& it_mstate : jmst.m_coefs)
     {
       auto& key = it_mstate.first;
       const float& m_coef = it_mstate.second;

       if( local_amplitudes[thread_num].find(key) == local_amplitudes[thread_num].end() ) local_amplitudes[thread_num][key] = std::vector<float>(imax,0.);
       for (size_t ilevel=0;ilevel<imax;++ilevel) local_amplitudes[thread_num][key][ilevel] += m_coef * nuvec.coefT[ilevel][istate];
     }
   }
   Profiler::timer["amplitudes_calc"] += omp_get_wtime() - t_start_inner;
   t_start_inner = omp_get_wtime();

   // Now accumulate amplitudes from all threads
   
   for (int ith=0;ith<nthreads;++ith)
   {
     for (auto& it_amp : local_amplitudes[ith] )
     {
       auto& key = it_amp.first;
       if( amplitudes.find(key) == amplitudes.end() ) amplitudes[key] = blank_vector;
       for (size_t ilevel=0;ilevel<imax;++ilevel) amplitudes[key][ivec][ilevel] +=  it_amp.second[ilevel];
     }
   }
   
   Profiler::timer["amplitudes_accumulate"] += omp_get_wtime() - t_start_inner;
  }// for ivec

  Profiler::timer["CalculateMschemeAmplitudes"] += omp_get_wtime() - t_start;
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
  int ma = Mf-Mi;

  double td = 0;

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


  int ia = m_index_a - ( j2_a -ma )/2;

  // loop through initial state amplitudes
  for (size_t iamp=0; iamp<amp_vec_i.size(); ++iamp)
  {
    auto& key = keys_i[iamp];
    auto& amp_i = amp_vec_i[iamp];

    if (key[ia]) continue; // if orbit a is occupied in the initial configuration, we get zero


    auto new_key = key;
    new_key.set(ia,1);
    if (amplitudes_f.find(new_key) == amplitudes_f.end() ) continue;


    // pick up a phase from commuting the a+ over to its place
    int phase_ladder = 1;
    for (int iphase=0;iphase<ia;++iphase) if( key[iphase]) phase_ladder *=-1;


    double amp_f = amplitudes_f[new_key][J_index_f][eigvec_f];

    td += amp_i * amp_f * phase_ladder;
  }
  
  td *= - sqrt((J2f+1.)/(j2_a+1.)) / clebsch_fi;  // minus sign from Wigner-Eckart convention of a phase (-1)^{2*lambda}, and lambda=ja is half-integer in this case.

  return td;
}











double TransitionDensity::TD_axaxa(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a, int m_index_b, int m_index_c, int J2ab, int Lambda2 )
{
  double td = 0;

  int J2i = Jlist_i[J_index_i];
  int J2f = Jlist_f[J_index_f];
  int Mi = J2i%2;
  int Mf = J2f%2;
  int mu = Mf-Mi;
  
  int j2_a = m_orbits[m_index_a].j2;
  int j2_b = m_orbits[m_index_b].j2;
  int j2_c = m_orbits[m_index_c].j2;
//  int j2_d = m_orbits[m_index_d].j2;

  // start out with the maximally-projected m state
  m_index_a += (j2_a - m_orbits[m_index_a].mj2)/2;
  m_index_b += (j2_b - m_orbits[m_index_b].mj2)/2;
  m_index_c += (j2_c - m_orbits[m_index_c].mj2)/2;
//  m_index_d += (j2_d - m_orbits[m_index_d].mj2)/2;


  // check some triangle conditions
  if ((J2f+Lambda2 < J2i) or (std::abs(J2f-Lambda2)>J2i)) return 0;
  if ((j2_a+j2_b<J2ab) or (std::abs(j2_a-j2_b)>J2ab)) return 0;
//  if ((j2_c+j2_d<J2cd) or (std::abs(j2_c-j2_d)>J2cd)) return 0;
  if ((J2ab+j2_c < Lambda2) or (std::abs(J2ab-j2_c)>Lambda2)) return 0;

  
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
//  if (m_index_c==m_index_d) norm *= SQRT2;  

  int Mab_min = std::max(-J2ab,mu-j2_c);
  int Mab_max = std::min( J2ab,mu+j2_c);
  for (int Mab=Mab_min;Mab<=Mab_max;Mab+=2)
  {
//    int Mcd = mu - Mab;
/*
    mc = Mab - mu;
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
*/
  }

  td *= sqrt((J2f+1.)/(Lambda2+1.)) / clebsch_fi * norm;




  return td;
}




arma::mat TransitionDensity::CalcOBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, Settings& settings)
{

  double t_start = omp_get_wtime();
  auto& jorbits = settings.jorbits;

  int Ji = Jlist_i[J_index_i];
  int Jf = Jlist_f[J_index_f];
  size_t njorb = jorbits.size();
  arma::mat obtd(njorb,njorb,arma::fill::zeros);


  for (size_t i=0; i<njorb; ++i)
  {
    int j2i = m_orbits[jorbits[i]].j2;
    int jmin = 0;
    if (settings.same_basename_fi and Ji==Jf and eigvec_i==eigvec_f) jmin = i;
    for (size_t j=jmin; j<njorb; ++j)
    {
      obtd(i,j) = OBTD( J_index_i, eigvec_i, J_index_f, eigvec_f, jorbits[i], jorbits[j], Lambda2);
      if (Lambda2 == 0)  obtd(i,j) /= sqrt( Ji+1 );

      if (settings.same_basename_fi and Ji==Jf and eigvec_i==eigvec_f)
      {
        int j2j = m_orbits[jorbits[j]].j2;
        obtd(j,i) = (1-std::abs(j2j-j2i)%4) * obtd(i,j);
      }
    }
  }

  Profiler::timer["CalcOBTD"] += omp_get_wtime() - t_start;
  return obtd;

}






arma::mat TransitionDensity::CalcTBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2, Settings& settings)
{

  double t_start = omp_get_wtime();
  // generate all the two body states that are needed
//   SetupKets();
  auto& ket_a = settings.ket_a;
  auto& ket_b = settings.ket_b;
  auto& ket_J = settings.ket_J;
  auto& jorbits = settings.jorbits;

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
    size_t iket_min = (settings.same_basename_fi and (Ji == Jf) and (eigvec_i==eigvec_f)) ? ibra : 0;
    for (size_t iket=iket_min;iket<ket_J.size();++iket)
    {
      int c = ket_a[iket];
      int d = ket_b[iket];
      int J2cd = ket_J[iket];
      tbtd(ibra,iket) = TBTD(  J_index_i,  eigvec_i,  J_index_f,  eigvec_f,
                               jorbits[a], jorbits[b],  jorbits[c], jorbits[d],
                                                          J2ab,  J2cd,  Lambda2 );

      // if final == initial, only calculate half of the matrix
      if  (settings.same_basename_fi and (Ji == Jf) and (eigvec_i==eigvec_f))
      {
        tbtd(iket,ibra) = tbtd(ibra,iket) * (1-std::abs(J2ab-J2cd)%4); 
      }
    }
  }

  if (Lambda2==0) tbtd /= sqrt( Ji+1.);

  Profiler::timer["CalcTBTD"] += omp_get_wtime() - t_start;
  return tbtd;
}



arma::vec TransitionDensity::CalcTransitionDensity_ax( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, Settings& settings)
{
  double t_start = omp_get_wtime();

  auto& jorbits = settings.jorbits;

  int Ji = Jlist_i[J_index_i];
  int Jf = Jlist_f[J_index_f];
  size_t njorb = jorbits.size();
  arma::vec td(njorb,arma::fill::zeros);
  std::cout << __func__ << "   Ji Jf eig_i, eig_f " << Ji << " " << Jf << " " << eigvec_i << " " << eigvec_f << std::endl;

//  std::ofstream densout( densfile_name, std::ios::app );
//  if ( densfile_name != "none")
//  {
//     densout << std::endl;
//     densout << "Jf nJf  Ji nJi = " << std::setw(3) << std::setprecision(1) << Jf*0.5 << " " << eigvec_f+1
//             << "    " << std::setw(3) << std::setprecision(1) << Ji*0.5 << " " << eigvec_i+1
//             << std::endl;
//     densout << "-------------- <f||a+||i> --------------------" << std::endl;
//  }

  for (size_t i=0; i<njorb; ++i)
  {
    int j2i = m_orbits[jorbits[i]].j2;
//    int jmin = 0;
//    if (basename_i==basename_f and Ji==Jf and eigvec_i==eigvec_f) jmin = i;
//    for (size_t j=jmin; j<njorb; ++j)
//    {
      td(i) = TD_ax( J_index_i, eigvec_i, J_index_f, eigvec_f, jorbits[i] );
//      if (Lambda2 == 0)  obtd(i,j) /= sqrt( Ji+1 );

//      if (densfile_name != "none" and std::abs(td(i))>1e-7)
//      {
//         densout << std::setw(3) << i << " "  << std::setw(12) << std::fixed << std::setprecision(8) << td(i)  << std::endl;
//      }

  }

//  densout.close();
  Profiler::timer["CalcTD_ax"] += omp_get_wtime() - t_start;
  return td;

}




arma::mat TransitionDensity::CalcTransitionDensity_axaxa( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, Settings& settings)
{
  arma::mat td;

  return td;

}









/*

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
*/




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







