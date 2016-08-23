#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>
#include <iomanip>
#include <vector>
#include <cstdint>
//#include <bitset>
#include <string>
#include <sstream>

#include <iostream>

#include "NuBasis.hh"
#include "NuProj.hh"
#include "NuVec.hh"
#include "JMState.hh"
#include "JBasis.hh"


using namespace std;


// some handy operator overloads for vectors
//
vector<float> operator*(const float lhs, const vector<float>& rhs)
{
  vector<float> vout = rhs;
  for (size_t i=0;i<vout.size();++i) vout[i] *= lhs;
  return vout;
}
vector<float>& operator+=(vector<float>& lhs, const vector<float>& rhs)
{
  for (size_t i=0;i<lhs.size();++i) lhs[i] += rhs[i];
  return lhs;
}





int main(int argc, char** argv)
{

ostringstream ostr;
vector<string> Afiles,Bfiles;
for (int i=0;i<5;++i)
{
 ostr << "ne200a" << i << "0B";
 Afiles.push_back( ostr.str() );
// cout << ostr.str() << endl;
 ostr.str("");
 ostr.str().clear();
}

for (int i=0;i<5;++i)
{
 ostr << "ne200b" << i << "0z";
 Bfiles.push_back( ostr.str() );
// cout << ostr.str() << endl;
 ostr.str("");
 ostr.str().clear();
}


int nprotons = 2;
int nneutrons = 2;
string interaction_id = "imsrg";
float hw = 20;
int Nshell = 6;
int Nmax = 0;
int Jtot = 4;
int MJtot = 0;
int parity = +1;


JBasis jbasis_0( "sdpn.sps", Afiles, Bfiles, Jtot, MJtot);


NuVec nuvec;
string vecfile = "ne200222.xvc";
nuvec.ReadFile(vecfile);
//cout << "Read file " << vecfile << ". Number of states = " << nuvec.no_state << "  Number of levels = " << nuvec.no_level << endl;




// create a hash table of the mstate basis
unordered_map<vector<mvec_type>,vector<float>,KeyHash> amplitudes;

// loop over J-coupled basis states
for (int istate=0;istate<nuvec.no_state;++istate)
{
   vector<float> level_coefs(nuvec.no_level,0.0);
   for (int ilevel=0;ilevel<nuvec.no_level;++ilevel)
   {
     level_coefs[ilevel] = nuvec.coefT[ilevel][istate];
   }

   JMState& jmst = jbasis_0.basis_states[istate];
   for (auto& it_mstate : jmst.m_coefs)
   {
     auto& key = it_mstate.first;
     float& m_coef = it_mstate.second;
     if( amplitudes.find(key) == amplitudes.end() ) amplitudes[key] = vector<float>(nuvec.no_level,0.0);
     amplitudes[key] += m_coef * level_coefs;
   }
}


cout << left << setw(10) << nprotons  << " !  number of protons" << endl;
cout << left << setw(10) << nneutrons << " !  number of neutrons" << endl;
cout << left << setw(max(10,(int)interaction_id.size()+3)) << interaction_id << " ! interaction id" << endl;
cout << left << setw(10) << hw << " ! hbar omega" << endl;
cout << left << setw(10) << Nshell << " ! N shells " << endl;
cout << left << setw(10) << jbasis_0.nubasis_a.m_orbits.size() +16 << " ! m-scheme orbits " << endl;
cout << left << setw(10) << Nmax << " ! Nmax " << endl;
cout << left << setw(10) << amplitudes.size() << " ! number of basis Slater determinants" << endl;
cout << left << setw(10) << showpos << parity << noshowpos << " ! parity " << endl;
cout << left << setw(10) << MJtot << " ! total MJ*2" << endl;
cout << left << setw(10) << nuvec.no_level << " ! Number of eigenstates" << endl;

//cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
//cout << "!!!  now list eigenvalues, J, T, deltaE !!!!!!" << endl;
//cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

// write out energies, J, T of eigenstates
for (int ilevel=0;ilevel<nuvec.no_level;++ilevel)
{
  cout << right << fixed << setw(12) << setprecision(4) << nuvec.alpha[ilevel] << " " << Jtot/2.0 << " " << 0 << " " << 0.0 << endl;
}

cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
cout << "!!!  now list mscheme single-particle basis !!!" << endl;
cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
// write single-particle basis in mscheme

int s=1;
for (int tz2=-1;tz2<=1;tz2+=2)
{
 for (int N=0;N<=1;++N)
 {
  for (int n=0;2*n<=N;++n)
  {
    int l2 = 2*(N-2*n);
    for (int j2=l2+1;j2>=max(0,l2-1);j2-=2)
    {
      for (int mj2=-j2;mj2<=j2;mj2+=2)
      {
      cout << setw(4) << s << " " << setw(3) << n << " " << setw(3) << l2 << " " << setw(3) << j2 << " " << setw(3) << mj2 << " " << setw(3) << tz2 << endl;;
      s++;
      }
    }
  }
 }
}
int s_core = s;
for (auto& morbit : jbasis_0.nubasis_a.m_orbits)
{
  cout << setw(4) << s << " " << setw(3) << morbit.n << " " << setw(3) << morbit.l2 << " " << setw(3) << morbit.j2 << " " << setw(3) << morbit.mj2 << " " << setw(3) << morbit.tz2 << endl;;
  s++;
}

cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
cout << "!!!  now list mscheme basis states and amplitudes !!!" << endl;
cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;


// Write m-scheme basis occupations and eigenvector coefficients
size_t bits_per_word = 8*sizeof(mvec_type);
for (auto it_amp : amplitudes)
{
  auto& mvec = it_amp.first;
  for (int i=0;i<s_core;++i) cout << 1 << " ";
  cout << " " ;
  for (size_t i=0;i<jbasis_0.nubasis_a.m_orbits.size();++i)
  {
    cout << ((mvec[i/bits_per_word] >> (i%bits_per_word))&0x1) << " ";
  }
//  cout << endl;
  cout << "   ";
  for (auto amp : it_amp.second)
  {
    cout << scientific << setw(13) << amp << " ";
  }
  cout << endl;

}


/*
  NuBasis nubasis_p, nubasis_n;
  NuProj  nuproj_p, nuproj_n;

  nubasis_p.ReadSPS("sdpn.sps");
  nubasis_n.ReadSPS("sdpn.sps");
  cout << "done reading sps" << endl;

  nubasis_p.ReadFile( "ne200a10B.nba" );
  nubasis_n.ReadFile( "ne200b10z.nba" );
  cout << "done reading nba" << endl;


  nuproj_p.ReadFile( "ne200a10B.prj" );
  nuproj_n.ReadFile( "ne200b10z.prj" );
  cout << "done reading prj" << endl;
  
  JMState jm_p( nubasis_p, nuproj_p, 0);
  JMState jm_n( nubasis_n, nuproj_n, 0);

  cout << "=========== jm_p ===============" << endl;
  jm_p.Print();

  cout << endl;
  cout << "=========== jm_n ===============" << endl;
  jm_n.Print();

  JMState jm_coupled = TensorProduct(jm_p,jm_n,0,0);
  jm_coupled.Print();

*/

/*
  nubasis.ReadSPS( argv[4] );
  nubasis.PrintSPS( );

  nubasis.ReadFile( argv[1] );
  nubasis.PrintBasis();

  cout << endl;
  cout << "================================================" << endl;
  cout << "========= Done with nubasis, start nuproj ======" << endl;
  cout << "================================================" << endl;
  cout << endl;

  NuProj nuproj;
  nuproj.ReadFile( argv[2] );
  nuproj.PrintProj();

  cout << endl;
  cout << "================================================" << endl;
  cout << "========= Done with nuproj, start nuvec ========" << endl;
  cout << "================================================" << endl;
  cout << endl;

  NuVec nuvec;
  nuvec.ReadFile( argv[3] );
  nuvec.PrintVectors();
  nuvec.CheckOrthoNormal();

  cout << endl;
  cout << "================================================" << endl;
  cout << "====== Done with nuvec, construct JMState ======" << endl;
  cout << "================================================" << endl;
  cout << endl;

  JMState jmstate( nubasis, nuproj, 2);
  jmstate.Print();
  JMState jmlower = jmstate.Jminus();
  jmlower.Print();
  JMState jmlower_upper = jmlower.Jplus();
  jmlower_upper.Print();
  JMState jmreversed = jmstate.TimeReverse();
  jmreversed.Print();

  JMState jmsum = jmstate + jmreversed;
  jmsum.Print();
  jmstate *= 2.0;
  jmstate.Print();

*/

  return 0;
}

