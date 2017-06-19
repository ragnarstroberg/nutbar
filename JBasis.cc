
#include <string>
#include <set>

#include "JBasis.hh"

//#define VERBOSE true

using namespace std;

JBasis::JBasis()
{}

JBasis::JBasis(int j2, int m2)
 : J2(j2), M2(m2)
{}

JBasis::JBasis( string sps_file,  vector<string> proton_files, vector<string> neutron_files, int j2, int m2 )
 : J2(j2), M2(m2)
{
  SetupBasis( sps_file, proton_files, neutron_files );
}


void JBasis::SetupBasis( string sps_file,  vector<string> A_files, vector<string> B_files  )
{
  #ifdef VERBOSE
    cout << "JBasis::SetupBasis -- begin. sps_file = " << sps_file << endl;
  #endif
  nubasis_a.ReadSPS( sps_file);
  nubasis_b.ReadSPS( sps_file);
  #ifdef VERBOSE
    cout << "JBasis::SetupBasis -- done reading sps_file " << endl;
  #endif

  for (string afile : A_files)
  {
    #ifdef VERBOSE
      cout << "JBasis::SetupBasis -- afile =  " << afile << endl;
     #endif
    nubasis_a.Clear();
    nuproj_a.Clear();
    nubasis_a.ReadFile(afile + ".nba");
    nuproj_a.ReadFile(afile + ".prj");
    int ngood_a = nuproj_a.ngood;
    for (string bfile : B_files)
    {
    #ifdef VERBOSE
      cout << "  JBasis::SetupBasis -- bfile =  " << bfile << endl;
     #endif
      nubasis_b.Clear();
      nuproj_b.Clear();
      nubasis_b.ReadFile(bfile + ".nba");
      nuproj_b.ReadFile(bfile + ".prj");
      int ngood_b = nuproj_b.ngood;


     #ifdef VERBOSE
      cout << "  JBasis::SetupBasis -- ngood_a, ngood_b =  " << ngood_a << " , " << ngood_b << endl;
     #endif
//      cout << "at begin, size of jmstates_a = " << jmstates_a.size() << "  and b = " << jmstates_b.size() << endl;
      int offset_a = jmstates_a.size();
      int offset_b = jmstates_b.size();
//      size_t offset_a = jmstates_a.size();
//      size_t offset_b = jmstates_b.size();
      for (int i_a=0;i_a<ngood_a;++i_a)
        jmstates_a.emplace_back( nubasis_a, nuproj_a, i_a);
      for (int i_b=0;i_b<ngood_b;++i_b)
        jmstates_b.emplace_back( nubasis_b, nuproj_b, i_b);

      for (int i_b=0;i_b<ngood_b;++i_b)
      {
//        JMState jmstate_b( nubasis_b, nuproj_b, i_b);
//        int JB = jmstate_b.J2;
//        cout << "i_b: " << i_b << "  JB= " << JB << "  or " << jmstates_b[i_b+offset_b].J2 << endl;
        int JB = jmstates_b[i_b+offset_b].J2;
        if (nubasis_b.ibf[nuproj_b.pindx[i_b]-1]<1) continue;
        for (int i_a=0;i_a<ngood_a;++i_a)
        {
//          JMState jmstate_a( nubasis_a, nuproj_a, i_a);
//          int JA = jmstate_a.J2;
          int JA = jmstates_a[i_a+offset_a].J2;
          if (nubasis_a.ibf[nuproj_a.pindx[i_a]-1]<1) continue;
          if (JA+JB<J2 or abs(JA-JB)>J2) continue;

         #ifdef VERBOSE
//          cout << "  JBasis::SetupBasis -- about to add TensorProduct( " << JA << ", " << JB << ", " << J2 << " " << M2 << " )  " << i_a << ", " << i_b << endl;
         #endif
//          basis_states.emplace_back( TensorProduct( jmstate_a, jmstate_b, J2,M2) );
//          basis_states.emplace_back( TensorProduct( jmstates_a[i_a+offset_a], jmstates_b[i_b+offset_b], J2,M2) );
          basis_states.emplace_back( array<int,4>({ i_a+offset_a, i_b+offset_b, J2, M2 }) );
         #ifdef VERBOSE
//          cout << "  JBasis::SetupBasis -- done. size of basis state = " << basis_states.back().m_coefs.size() << endl;
         #endif

        }
      }

    }
  }
}


JMState JBasis::GetBasisState(size_t index) const
{
  int i_a = basis_states[index][0];
  int i_b = basis_states[index][1];
  int J2 = basis_states[index][2];
  int M2 = basis_states[index][3];
  return TensorProduct( jmstates_a[i_a], jmstates_b[i_b], J2, M2 );
//  return basis_states[index];
}



int JBasis::GetNaiveMschemeDimension() const
{
//  set<key_type> mstates_a, mstates_b;
  map<int,set<uint64_t>> mstates_a, mstates_b;
  for (auto& jst : jmstates_a)
  {
    if (mstates_a.count(jst.M2)<1) mstates_a[jst.M2] = set<uint64_t>();
    for (auto& itm: jst.m_coefs)
    {
      uint64_t k = itm.first.to_ulong();
//      if (mstates_a.count(itm.first)<1) mstates_a.insert(itm.first);
      if (mstates_a[jst.M2].count(k)<1 and jst.M2==3)
      {
        mstates_a[jst.M2].insert(k);
        cout << "Inserted " << itm.first <<  " J,M = (" << jst.J2 <<  " " << jst.M2 << ")   k = " << k << endl;
      }
    }
  }
  for (auto& jst : jmstates_b)
  {
    if (mstates_b.count(jst.M2)<1) mstates_b[jst.M2] = set<uint64_t>();
    for (auto& itm: jst.m_coefs)
    {
      uint64_t k = itm.first.to_ulong();
      if (mstates_b[jst.M2].count(k)<1) mstates_b[jst.M2].insert(k);
//      if (mstates_b.count(itm.first)<1) mstates_b.insert(itm.first);
    }
  }
  int nm = 0;
  cout << "J2,M2 = " << J2 << " " << M2 << endl;
  for (auto& it_a : mstates_a)
  {
    int ma = it_a.first;
    int mb = abs(J2 - ma);
    cout << "ma = " << ma << endl;
    if (mstates_b.count(mb))
    {
      cout << "ma = " << ma << "  mb = " << mb << "  sizes = " << mstates_a[ma].size() << " " << mstates_b[mb].size() << endl;
      nm += (min(ma,mb)+1)*it_a.second.size() * mstates_b[mb].size();
      for ( auto vecb : mstates_b[mb]) cout << bitset<64>(vecb) << endl;
    }
  }
//  cout << "dim a: " << mstates_a.size() << "   dim b: " << mstates_b.size() << "  product = " << mstates_a.size() * mstates_b.size() << endl;
  cout << "dim = " << nm << endl;
  return mstates_a.size() * mstates_b.size();

}

