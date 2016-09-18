//
//  Program nutbar  (NuShellX Transtitions from Binary Arrays)
//
//

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
//cout << "obtd:" << endl << trans.CalcOBTD(0,0,0,0,0) << endl;
//cout << "obtd:" << endl << obtd << endl;

int Lambda2 = 4;
arma::mat OpScalar1b;
arma::mat OpScalar2b;
//vector<double> Janky_2J = { 3, 5, 1, 3, 5 , 1 };
vector<double> Janky_2J = { 5, 3, 1, 5, 3 , 1 };

trans.GetScalarTransitionOperator("testfiles/IMSRG.int",OpScalar1b,OpScalar2b);
//trans.GetScalarTransitionOperator("testfiles/usdbpn.int",OpScalar1b,OpScalar2b);

for (int Ji : {0,1} )
{
 for (int Jf : {0,1} )
// for (int Jf : {0} )
 {
  for (int ivec = 0; ivec<=1; ++ivec)
  {
   for (int fvec = 0; fvec<=1; ++fvec)
   {
    if (Ji==Jf and fvec>ivec) continue;
    cout << Ji << " " << ivec << " ->  " << Jf << " " << fvec << endl;
    arma::mat obtd,tbtd;
    obtd = trans.CalcOBTD(Ji,ivec,Jf,fvec,Lambda2);
    tbtd = trans.CalcTBTD(Ji,ivec,Jf,fvec,Lambda2);
    cout << "obtd" << endl << obtd << endl;
    cout << "<Op1b>" << arma::accu( Op1b % obtd ) << endl;
    cout << "<Op2b>" << arma::accu( Op2b % tbtd) << endl;
//    double sum_op = 0;
//    for (int i=obtd.n_rows-1;i>=0;--i)
//    {
//     for (int j=obtd.n_rows-1;j>=0;--j)
//     {
//      sum_op += Op1b(i,j) * obtd(i,j);
//      if (abs(Op1b(i,j))>1e-6)
//      cout << i << " " << j << " : " << fixed << setprecision(8) << setw(14) << Op1b(i,j) << " " << setprecision(8) << setw(14) <<  obtd(i,j) << "  " << setprecision(8) << setw(14) << Op1b(i,j) * obtd(i,j) << "   " << setprecision(8) << setw(10) << sum_op << endl;
//     }
//    }
    if (Ji==Jf )
    {
      obtd = trans.CalcOBTD(Ji,ivec,Jf,fvec,0);
      tbtd = trans.CalcTBTD(Ji,ivec,Jf,fvec,0);
//      cout << "occupations:" << obtd.diag().t() << endl;
      double sum_occ=0;
      cout << "occupations:";
      for (int i=0; i<obtd.n_rows; ++i)
      {
        sum_occ += obtd(i,i) * sqrt( Janky_2J[i]+1 ) / sqrt(4*Jf+1.);
        cout   <<  obtd(i,i) * sqrt( Janky_2J[i]+1 ) / sqrt(4.*Jf+1.)<< "  ";
      }
      cout << "   sum = " << sum_occ  << endl;
      double scop1 = arma::accu( OpScalar1b % obtd ) / sqrt(4*Jf+1.) ;
      double scop2 = arma::accu( OpScalar2b % tbtd ) / sqrt(4*Jf+1.);
//      double summed_me = 0;
//      for (size_t i=0; i< tbtd.n_rows; ++i)
//      {
//        for (size_t j=0; j<tbtd.n_rows; ++j)
//        {
//          summed_me += tbtd(i,j) * OpScalar2b(i,j) ;
//          if ( abs(tbtd(i,j)*OpScalar2b(i,j))>1e-6)
//          cout << i << " " << j << ":  " << setw(10) << setprecision(6) << tbtd(i,j) << "  " << setw(10) << setprecision(6) << OpScalar2b(i,j)
//               << "  " << setw(10) << setprecision(6) << summed_me
//               << endl;
//        }
//      }
      cout << "<ScalarOp1b>" << scop1  << endl;
      cout << "<ScalarOp2b>" << scop2 << endl;
      cout << "<ScalarOp1b+2b>" << scop1+scop2 << endl;
    }

   }
  }
 }
}

//cout << "Op2b " << endl << Op2b << endl;
//double me2b = 0;
//for (size_t ibra=0;ibra<Op2b.n_cols;++ibra)
//{
//  for (size_t iket=0;iket<Op2b.n_rows;++iket)
//  {
//    me2b += Op2b(ibra,iket) * tbtd(ibra,iket);
//    if (abs(Op2b(ibra,iket))>1e-8 and abs(tbtd(ibra,iket))>1e-8)
//     cout << setw(3) << ibra << " " << setw(3) << iket << " " << setw(12) << setprecision(8) << fixed << Op2b(ibra,iket) << " " << setw(12) << setprecision(8) << tbtd(ibra,iket) 
//    << setw(12) << setprecision(8) << me2b << endl;
//  }
//}
//cout << trans.TBTD(0,0,0,0,11,11,11,11,0,0,0) << endl;
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

