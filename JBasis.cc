
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
      for (int i_b=0;i_b<ngood_b;++i_b)
      {
        JMState jmstate_b( nubasis_b, nuproj_b, i_b);
        int JB = jmstate_b.J2;
        if (nubasis_b.ibf[nuproj_b.pindx[i_b]-1]<1) continue;
        for (int i_a=0;i_a<ngood_a;++i_a)
        {
          JMState jmstate_a( nubasis_a, nuproj_a, i_a);
          int JA = jmstate_a.J2;
          if (nubasis_a.ibf[nuproj_a.pindx[i_a]-1]<1) continue;
          if (JA+JB<J2 or abs(JA-JB)>J2) continue;
//          cout << endl << "basis state " << basis_states.size()
//               << "  JA,JB = " << JA << "," << JB
//               << "  ngood = " << ngood_a << " " << ngood_b
//               << " ibf = " << nubasis_a.ibf[nuproj_a.pindx[i_a]-1] << " " << nubasis_b.ibf[nuproj_b.pindx[i_b]-1] 
//               << "   norm = " << jmstate_a.Norm() << " " << jmstate_b.Norm() << "  = > "
//               << endl;
         #ifdef VERBOSE
          cout << "  JBasis::SetupBasis -- about to add TensorProduct( " << JA << ", " << JB << ", " << J2 << " " << M2 << " )  " << i_a << ", " << i_b << endl;
         #endif
          basis_states.emplace_back( TensorProduct( jmstate_a, jmstate_b, J2,M2) );
         #ifdef VERBOSE
          cout << "  JBasis::SetupBasis -- done. size of basis state = " << basis_states.back().m_coefs.size() << endl;
         #endif
//          cout << "basis state " << basis_states.size()-1
//               << "  JA,JB = " << JA << "," << JB
//               << "  ngood = " << ngood_a << " " << ngood_b
//               << " ibf = " << nubasis_a.ibf[i_a] << " " << nubasis_b.ibf[i_b] 
//               << "   norm = " << jmstate_a.Norm() << " " << jmstate_b.Norm() << "  = > "
//               << basis_states.back().Norm()
//               << endl;
        }
      }

    }
  }
//  cout << "Set up basis with Jtot = " << J2 << ". Number of basis states = " << basis_states.size() << endl;
}
