
#include <string>

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

// This is the old version. It worked for most cases, but I think it breaks for 76Ge 0vBB decay...
/*
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
      {
        jmstates_a.emplace_back( nubasis_a, nuproj_a, i_a);
//        cout << "Added jmstate_a  " << jmstates_a.size() << "  i_a = " << i_a << endl;
      }
     #ifdef VERBOSE
      cout << "  JBasis::SetupBasis -- done looping of a and emplacing jmstates"  << endl;
     #endif
      for (int i_b=0;i_b<ngood_b;++i_b)
      {
        jmstates_b.emplace_back( nubasis_b, nuproj_b, i_b);
//        cout << "Added jmstate_b  " << jmstates_b.size() << "  i_b = " << i_b << endl;
      }

      #ifdef VERBOSE
       cout << "  JBasis::SetupBasis -- begin looping over ngood_b,ngood_a" << endl;
      #endif
      for (int i_b=0;i_b<ngood_b;++i_b)
      {
//        JMState jmstate_b( nubasis_b, nuproj_b, i_b);
//        int JB = jmstate_b.J2;
//        cout << "i_b: " << i_b << "  JB= " << JB << "  or " << jmstates_b[i_b+offset_b].J2 << endl;
        int JB = jmstates_b[i_b+offset_b].J2;
        if (nubasis_b.ibf[nuproj_b.pindx[i_b]-1]<1) continue;
        cout << "inner i_a loop: ngood_a = " << ngood_a << endl;
        for (int i_a=0;i_a<ngood_a;++i_a)
        {
//          JMState jmstate_a( nubasis_a, nuproj_a, i_a);
//          int JA = jmstate_a.J2;
        cout << "ia,ib: " << i_a << "," << i_b << endl;
          int JA = jmstates_a[i_a+offset_a].J2;
          if (nubasis_a.ibf[nuproj_a.pindx[i_a]-1]<1) continue;
          if (JA+JB<J2 or abs(JA-JB)>J2) continue;
        cout << "   ia,ib: " << i_a << "," << i_b << endl;

//          basis_states.emplace_back( TensorProduct( jmstate_a, jmstate_b, J2,M2) );
//          basis_states.emplace_back( TensorProduct( jmstates_a[i_a+offset_a], jmstates_b[i_b+offset_b], J2,M2) );
          basis_states.emplace_back( array<int,4>({ i_a+offset_a, i_b+offset_b, J2, M2 }) );

        cout << " basis state: " << basis_states.size()-1 << "   " << i_a+offset_a << " " << jmstates_a[i_a+offset_a].J2 << "  " << jmstates_a[i_a+offset_a].pindx[0] << "  x  "
             << i_b+offset_b << " " << jmstates_b[i_b+offset_b].J2 << "  " << jmstates_b[i_b+offset_b].pindx[0] << "   -> " << J2 << " " << M2 << endl;
        }
      }
         #ifdef VERBOSE
          cout << "  JBasis::SetupBasis -- done looping over ngood_b,ngood_a" << endl;
         #endif

    }
  }
  cout << "Size of jmstates_a = " << jmstates_a.size() << endl;
  cout << "Size of jmstates_b = " << jmstates_b.size() << endl;
}
*/


void JBasis::SetupBasis( string sps_file,  vector<string> A_files, vector<string> B_files  )
{
  #ifdef VERBOSE
    cout << "JBasis::SetupBasis -- begin. sps_file = " << sps_file << endl;
  #endif
  NuBasis nubasis;
  nubasis.ReadSPS( sps_file);
  #ifdef VERBOSE
    cout << "JBasis::SetupBasis -- done reading sps_file " << endl;
  #endif

  vector<int> offsets_a;
  vector<int> offsets_b;
  SetUpJMState_ab( sps_file, A_files, jmstates_a, offsets_a);
  SetUpJMState_ab( sps_file, B_files, jmstates_b, offsets_b);
  

  for (int ind_Afile=0; ind_Afile<A_files.size(); ++ind_Afile)
  {
   int ia_min = offsets_a[ind_Afile];
   int ia_max = (ind_Afile<A_files.size()-1) ? offsets_a[ind_Afile+1] : A_files.size();
   for (int ind_Bfile=0; ind_Bfile<B_files.size(); ++ind_Bfile)
   {

     int ib_min = offsets_b[ind_Bfile];
     int ib_max = (ind_Bfile<B_files.size()-1) ? offsets_b[ind_Bfile+1] : B_files.size();
     for (int i_b=ib_min; i_b<ib_max; ++i_b)
     {
      int JB = jmstates_b[i_b].J2;
      for (int i_a=ia_min; i_a<ia_max; ++i_a)
      {
        int JA = jmstates_a[i_a].J2;
        if (JA+JB<J2 or abs(JA-JB)>J2) continue;
        basis_states.emplace_back( array<int,4>({ i_a, i_b, J2, M2 }) );
      }
     }
   }
  }

}

void JBasis::SetUpJMState_ab(  string sps_file,  vector<string> filenames, vector<JMState>& jmstates,  vector<int>& offsets)
{
  NuBasis nubasis;
  nubasis.ReadSPS( sps_file);
  m_orbits = nubasis.m_orbits;
  for (string afile : filenames)
  {
    #ifdef VERBOSE
      cout << "JBasis::SetupBasis -- afile =  " << afile << endl;
     #endif
    NuProj nuproj;
    nubasis.ReadFile(afile + ".nba");
    nuproj.ReadFile(afile + ".prj");
    m_orbits = nubasis.m_orbits;
    int ngood = nuproj.ngood;

    offsets.push_back( jmstates.size() );
    for (int i=0;i<ngood;++i)
    {
      if (nubasis.ibf[nuproj.pindx[i]-1]<1) continue;
      jmstates.emplace_back( nubasis, nuproj, i);
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

