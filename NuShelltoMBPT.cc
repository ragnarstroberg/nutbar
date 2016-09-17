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

cout << "about to read transition operator" << endl;
arma::mat Op1b = trans.GetOneBodyTransitionOperator("testfiles/IMSRG_E2_1b.op");
arma::mat Op2b = trans.GetTwoBodyTransitionOperator("testfiles/IMSRG_E2_2b.op");
cout << "Op1b " << endl << Op1b << endl;

//cout << "obtd:" << endl << trans.CalcOBTD(0,0,0,0,0) << endl;
cout << "obtd:" << endl << trans.CalcOBTD(0,0,0,0,0) << endl;
arma::mat obtd = trans.CalcOBTD(1,0,1,0,4);
cout << "obtd:" << endl << obtd << endl;



arma::mat tbtd = trans.CalcTBTD(1,0,1,0,4);
//cout << "Op2b " << endl << Op2b << endl;
double me2b = 0;
for (size_t ibra=0;ibra<Op2b.n_cols;++ibra)
{
  for (size_t iket=0;iket<Op2b.n_rows;++iket)
  {
    me2b += Op2b(ibra,iket) * tbtd(ibra,iket);
    if (abs(Op2b(ibra,iket))>1e-8 and abs(tbtd(ibra,iket))>1e-8)
     cout << setw(3) << ibra << " " << setw(3) << iket << " " << setw(12) << setprecision(8) << fixed << Op2b(ibra,iket) << " " << setw(12) << setprecision(8) << tbtd(ibra,iket) 
    << setw(12) << setprecision(8) << me2b << endl;
  }
}
cout << "<Op1b>" << arma::accu( Op1b % obtd ) << endl;
cout << "<Op2b>" << arma::accu( Op2b % tbtd) << endl;
//cout << "obtd:" << endl << trans.CalcOBTD(0,0,1,1,4) << endl;
//
//trans.CalcTBTD(0,0,1,0,4);
//
//int J_index_i = 1;
//for (int J_index_f : {0, 1})
//{
//  for (int eigvec_i=0;eigvec_i<2;++eigvec_i)
//  {
//    for (int eigvec_f=0;eigvec_f<2;++eigvec_f)
//    {
//      int m_index_a = 0;
//      int m_index_b = 5;
//      int m_index_c = 0;
//      int m_index_d = 10;
//      int J2ab = 4;
//      int J2cd = 4;
//      int Lambda2 = 4;
//      double tbd = trans.TBTD( J_index_i,  eigvec_i,  J_index_f,  eigvec_f,  m_index_a,  m_index_b,  m_index_c,  m_index_d,  J2ab,  J2cd,  Lambda2 );
//      cout << J_index_i << " " << eigvec_i << " " << J_index_f << " " << eigvec_f << " " << m_index_a << " " << m_index_b << " " << m_index_c << " " << m_index_d << " " << J2ab << " " << J2cd << " " << Lambda2 << "     " << tbd << endl;
//    }
//  }
//}

  return 0;
}

