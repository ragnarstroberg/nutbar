
#include <string>

#include "JBasis.hh"

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
  nubasis_a.ReadSPS( sps_file);
  nubasis_b.ReadSPS( sps_file);

  int JAmax = 2*(A_files.size()-1);
  int JBmax = 2*(B_files.size()-1);

  int JAmin = max(0,J2-JBmax);

  for (int JA = JAmin; JA<=JAmax;JA+=2)
  {
//    cout << endl;
//    cout << "================== basis size = " << basis_states.size() << " =================" << endl;
//    cout << "JA = " << JA << "  file = " << A_files[JA/2] << endl;
    nubasis_a.Clear();
    nuproj_a.Clear();
    nubasis_a.ReadFile(A_files[JA/2] + ".nba");
    nuproj_a.ReadFile(A_files[JA/2] + ".prj");
//    cout << "Read the A files" << endl;
    int ngood_a = nuproj_a.ngood;
//    cout << "ngood_a = " << ngood_a << endl;
    int JBmin = abs(J2-JA);
    JBmax = min(J2+JA,2*((int)B_files.size()-1));
//    cout << " JB min/max = " << JBmin << " " << JBmax << endl;
    for (int JB = JBmin; JB<=JBmax;JB+=2)
    {
      nubasis_b.Clear();
      nuproj_b.Clear();
//      cout << " JB = " << JB << "  file = " << B_files[JB/2] << endl;
      nubasis_b.ReadFile(B_files[JB/2] + ".nba");
      nuproj_b.ReadFile(B_files[JB/2] + ".prj");
//      cout << " Read the A files" << endl;
      int ngood_b = nuproj_b.ngood;
//      cout << " ngood_b = " << ngood_b << endl;
   
      for (int i_b=0;i_b<ngood_b;++i_b)
      {
        JMState jmstate_b( nubasis_b, nuproj_b, i_b);
        if (jmstate_b.J2!= JB)
          cout << "AAAH, trouble matching JB: " << jmstate_b.J2 << " != " << JB << endl;
        for (int i_a=0;i_a<ngood_a;++i_a)
        {
          JMState jmstate_a( nubasis_a, nuproj_a, i_a);
          if (jmstate_a.J2!= JA)
            cout << "AAAH, trouble matching JA: " << jmstate_a.J2 << " != " << JA << endl;
          basis_states.emplace_back( TensorProduct( jmstate_a, jmstate_b, J2,M2) );
        }
      }
    }
  }
//  cout << "Set up basis with Jtot = " << J2 << ". Number of basis states = " << basis_states.size() << endl;
}
