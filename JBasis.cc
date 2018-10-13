
#include <string>
#include <vector>
#include <set>
#include <map>

#include "JBasis.hh"

//#define VERBOSE true

//using namespace std;

JBasis::JBasis()
{}



JBasis::JBasis( std::string sps_file,  std::vector<std::string> proton_files, std::vector<std::string> neutron_files, std::vector<int>& J2list)
// : J2(j2), M2(m2)
{
  SetupBasis( sps_file, proton_files, neutron_files, J2list );
}



void JBasis::SetupBasis( std::string sps_file,  std::vector<std::string> A_files, std::vector<std::string> B_files, std::vector<int>& J2list  )
{
  #ifdef VERBOSE
    cout << "JBasis::SetupBasis -- begin. sps_file = " << sps_file << endl;
  #endif
  NuProj nuproj;
  NuBasis nubasis;
  nubasis.ReadSPS( sps_file);
  m_orbits = nubasis.m_orbits;
  #ifdef VERBOSE
    cout << "JBasis::SetupBasis -- done reading sps_file " << endl;
  #endif

  std::vector<int> offsets_a;
  std::vector<int> offsets_b;

//  std::cout << "   Calling SetUpJMState_ab" << std::endl << std::endl;

  SetUpJMState_ab( nubasis, nuproj, A_files, jmstates_a, offsets_a);
  SetUpJMState_ab( nubasis, nuproj, B_files, jmstates_b, offsets_b);
  
  for (int J2 : J2list) AddBasisStates_J( A_files.size(), B_files.size(), offsets_a, offsets_b, J2);
}


void JBasis::AddBasisStates_J( size_t nA, size_t nB, std::vector<int>& offsets_a, std::vector<int>& offsets_b, int J2  )
{
  

//  std::cout << "     Adding a basis state with nA,nB = " << nA << " " << nB << std::endl;
//  std::cout << "    The sizes of the jmstates are  " << jmstates_a.size() << "   " << jmstates_b.size() << std::endl;
//  std::cout << "   the offsets are " << std::endl;
//  for (auto aoff : offsets_a ) std::cout << aoff << " ";
//  std::cout << std::endl;
//  for (auto aoff : offsets_b ) std::cout << aoff << " ";
//  std::cout << std::endl;
  basis_states[J2] = {};
  // Loop over A files, which contain, e.g. proton configurations for each proton J.
  for (int ind_Afile=0; ind_Afile<nA; ++ind_Afile)
  {
   // i_a loops over A-type basis states
   int ia_min = offsets_a[ind_Afile];
//   int ia_max = (ind_Afile<nA-1) ? offsets_a[ind_Afile+1] : nA;
//   int ia_max = (ind_Afile<nA-1) ? offsets_a[ind_Afile+1] : nA-1;
   int ia_max = (ind_Afile<nA-1) ? offsets_a[ind_Afile+1] : std::min(nA,jmstates_a.size()-1);
   // this weird ordering is dictated by the NuShellX data format. Sorry...
   for (int ind_Bfile=0; ind_Bfile<nB; ++ind_Bfile)
   {
     int ib_min = offsets_b[ind_Bfile];
//     int ib_max = (ind_Bfile<nB-1) ? offsets_b[ind_Bfile+1] : nB;
//     int ib_max = (ind_Bfile<nB-1) ? offsets_b[ind_Bfile+1] : nB-1;
     int ib_max = (ind_Bfile<nB-1) ? offsets_b[ind_Bfile+1] : std::min(nB,jmstates_b.size()-1);
//     std::cout << " looping ib_min to ib_max " << ib_min << " " << ib_max << std::endl;
     for (int i_b=ib_min; i_b<ib_max; ++i_b)
     {
//      int JB = jmstates_b[i_b].J2;
      int JB = jmstates_b.at(i_b).J2;
      for (int i_a=ia_min; i_a<ia_max; ++i_a)
      {
        int JA = jmstates_a[i_a].J2;
        if (JA+JB<J2 or abs(JA-JB)>J2) continue;
//        std::cout << " About to push back a basis state with ia,ib = " << i_a << " " << i_b << std::endl;
        basis_states[J2].emplace_back( ab_pair(i_a,i_b));
      }
     }
   }
  }
}



void JBasis::SetUpJMState_ab( NuBasis& nubasis, NuProj& nuproj,  std::vector<std::string> filenames, std::vector<JMState>& jmstates,  std::vector<int>& offsets)
{
//  std::cout << "In " << __func__ << "  looping over filenames. There are " << filenames.size() << "  of them " << std::endl;
  for (std::string afile : filenames)
  {
    #ifdef VERBOSE
      std::cout << "JBasis::SetupBasis -- afile =  " << afile << std::endl;
     #endif
//    std::cout << "Reading from " << afile  << " .npa, .prj" << std::endl;
    nubasis.ReadFile(afile + ".nba");
    nuproj.ReadFile(afile + ".prj");
    int ngood = nuproj.ngood;

//    std::cout << "ngood = " << ngood << std::endl;

    offsets.push_back( jmstates.size() );
    for (int i=0;i<ngood;++i)
    {
      if (nubasis.ibf[nuproj.pindx[i]-1]<1) continue;
//      std::cout << "  adding jmstate i = " << i << " .  size of jmstates = " << jmstates.size() << std::endl;
      jmstates.emplace_back( nubasis, nuproj, i);
//      std::cout << "    we have "  << jmstates.back().m_coefs.size() <<  " mstates " << std::endl;
//      std::cout << "  this state has M2 = " << jmstates.back().M2 << std::endl;
    }
  }
}

JMState JBasis::GetBasisState(size_t index, int J2, int M2) const
{
  const int i_a = basis_states.at(J2).at(index).ia;
  const int i_b = basis_states.at(J2).at(index).ib;
//  std::cout << "    In " << __func__ << "  calling tensor product with  ia,ib = " <<  i_a << " " << i_b << " j,m =  " << jmstates_a[i_a].J2 << " " << jmstates_a[i_a].M2 << ", " << jmstates_b[i_b].J2 << " " << jmstates_b[i_b].M2 << std::endl;
  return TensorProduct( jmstates_a[i_a], jmstates_b[i_b], J2, M2 );
}



int JBasis::GetNaiveMschemeDimension(int J2) const
{
  std::map<int,std::set<uint64_t>> mstates_a, mstates_b;
  for (auto& jst : jmstates_a)
  {
    if (mstates_a.count(jst.M2)<1) mstates_a[jst.M2] = std::set<uint64_t>();
    for (auto& itm: jst.m_coefs)
    {
      uint64_t k = itm.first.to_ulong();
//      if (mstates_a.count(itm.first)<1) mstates_a.insert(itm.first);
      if (mstates_a[jst.M2].count(k)<1 and jst.M2==3)
      {
        mstates_a[jst.M2].insert(k);
      }
    }
  }
  for (auto& jst : jmstates_b)
  {
    if (mstates_b.count(jst.M2)<1) mstates_b[jst.M2] = std::set<uint64_t>();
    for (auto& itm: jst.m_coefs)
    {
      uint64_t k = itm.first.to_ulong();
      if (mstates_b[jst.M2].count(k)<1) mstates_b[jst.M2].insert(k);
    }
  }
  int nm = 0;
  for (auto& it_a : mstates_a)
  {
    int ma = it_a.first;
    int mb = abs(J2 - ma);
    if (mstates_b.count(mb))
    {
      nm += (std::min(ma,mb)+1)*it_a.second.size() * mstates_b[mb].size();
    }
  }
  return mstates_a.size() * mstates_b.size();

}

