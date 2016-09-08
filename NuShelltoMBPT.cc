#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <string>

#include "NuBasis.hh"
#include "NuProj.hh"
#include "NuVec.hh"
#include "JMState.hh"
#include "JBasis.hh"
#include "TransitionDensity.hh"


using namespace std;






int main(int argc, char** argv)
{

TransitionDensity trans;

if (argc>1)
{
  trans.ReadInputFromFile( argv[1]);
}
else
{
  trans.ReadInputInteractive();
}


trans.ReadFiles();

trans.CalculateMschemeAmplitudes();

trans.WriteEGV("mbpt.egv");
trans.WriteTRDENS_input("trdens.in");

cout << "obtd:" << endl << trans.CalcOBTD(0,0,0,0,0) << endl;
cout << "obtd:" << endl << trans.CalcOBTD(0,0,0,1,0) << endl;
cout << "obtd:" << endl << trans.CalcOBTD(0,0,1,0,4) << endl;
cout << "obtd:" << endl << trans.CalcOBTD(0,0,1,1,4) << endl;

trans.CalcTBTD(0,0,1,0,4);

int J_index_i = 1;
for (int J_index_f : {0, 1})
{
  for (int eigvec_i=0;eigvec_i<2;++eigvec_i)
  {
    for (int eigvec_f=0;eigvec_f<2;++eigvec_f)
    {
      int m_index_a = 0;
      int m_index_b = 5;
      int m_index_c = 0;
      int m_index_d = 10;
      int J2ab = 4;
      int J2cd = 4;
      int Lambda2 = 4;
      double tbd = trans.TBTD( J_index_i,  eigvec_i,  J_index_f,  eigvec_f,  m_index_a,  m_index_b,  m_index_c,  m_index_d,  J2ab,  J2cd,  Lambda2 );
      cout << J_index_i << " " << eigvec_i << " " << J_index_f << " " << eigvec_f << " " << m_index_a << " " << m_index_b << " " << m_index_c << " " << m_index_d << " " << J2ab << " " << J2cd << " " << Lambda2 << "     " << tbd << endl;
    }
  }
}

  return 0;
}

