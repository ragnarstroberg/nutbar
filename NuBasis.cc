
#include <cmath>
#include <numeric>
#include "NuBasis.hh"

#define DELIMITER_LENGTH 4

const int NuBasis::gwords = 1; // this is good if the modelspace has <=64 m-states. Otherwise, we may need to do some tinkering.


NuBasis::NuBasis()
{}


void NuBasis::ReadFile(string fname)
{
  Clear();
  ifstream infile(fname, ios::binary );
  if (not infile.good() )
  {
   cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
   cout << "!! ERROR: Trouble in NuBasis::ReadFile --" << fname << "  !!" << endl;
   cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    exit(EXIT_FAILURE) ;
  }

 int nwords;
 sum_nJT=0;

 infile.read((char*)&nwords, sizeof(nwords));
// cout << "nwords = " << nwords << endl;
 infile.read((char*)&no_spart, sizeof(no_spart));
 infile.read((char*)&nwords, sizeof(nwords));
// cout << "nwords = " << nwords << endl;
// cout <<"no_spart: " << no_spart << endl;
 pindx.resize(no_spart);
 nJT.resize(no_spart);
 j.resize(no_spart);
 t.resize(no_spart);
 ibf.resize(no_spart);
 max_cM2.resize(no_spart);
 max_cTz2.resize(no_spart);
 partition.resize(no_spart);
 vec.resize(no_spart);

 for (int isp=0;isp<no_spart;isp++)
 {

    infile.read((char*)&nwords, sizeof(nwords));
    infile.read((char*)&pindx[isp], sizeof(pindx[isp]));
    infile.read((char*)&nJT[isp], sizeof(nJT[isp]));
    infile.read((char*)&j[isp], sizeof(j[isp]));
    infile.read((char*)&t[isp], sizeof(t[isp]));
    infile.read((char*)&ibf[isp], sizeof(ibf[isp]));
    infile.read((char*)&max_cM2[isp], sizeof(max_cM2[isp]));
    infile.read((char*)&max_cTz2[isp], sizeof(max_cTz2[isp]));
    infile.read((char*)&nwords, sizeof(nwords));

//    cout << "pindx: " << pindx[isp] << endl;
//    cout << "nJT: " << nJT[isp] << endl;
//    cout << "j: " << j[isp] << endl;
//    cout << "t: " << t[isp] << endl;
//    cout << "ibf: " << ibf[isp] << endl;
//    cout << "max_cM2: " << max_cM2[isp] << endl;
//    cout << "max_cTz2: " << max_cTz2[isp] << endl;
    sum_nJT += nJT[isp];
   
    infile.read((char*)&nwords, sizeof(nwords));
//    cout << "nwords = " << nwords << endl;
    partition[isp].resize(nwords/sizeof(part_type));
    infile.read((char*)&partition[isp][0], nwords);
    infile.read((char*)&nwords, sizeof(nwords));
   // infile.seekg(DELIMITER_LENGTH,ios::cur);
   
//    cout << "parition: ";
//    for (auto c : partition[isp]) cout << int(c) << " ";
//    cout << endl;
    
    if (ibf[isp]>0)
    {    
    infile.read((char*)&nwords, sizeof(nwords));
//    cout << "for vec, nwords = " << nwords  << "  compare with " << gwords*sizeof(mvec_type)*ibf[isp] << endl;
      vec[isp].resize(ibf[isp]);
//      cout << "vec:" << endl;
      for (int i=0;i<ibf[isp];i++) 
      {
//        vec[isp][i].resize(gwords);
//        infile.read((char*)&vec[isp][i][0], gwords*sizeof(mvec_type));
        infile.read((char*)&vec[isp][i], gwords*sizeof(mvec_type));
//        cout << "  i: " << i << endl;
//        for (int g=0;g<gwords;++g)
//        {
//           cout << "  " << setw(10) << vec[isp][i][g] << "   " << bitset<8*gwords*sizeof(mvec_type)>(vec[isp][i][g]) << endl;;
//        }
//        cout << endl;
      }
      infile.read((char*)&nwords, sizeof(nwords));
//      cout << "done with vec. nwords = " << nwords << endl;
    }
//    cout << endl;
 }
// cout << "sum nJT = " << sum_nJT << endl;

}






void NuBasis::ReadSPS(string fname)
{

  ifstream infile( fname );
  if ( not infile.good() )
  {
   cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
   cout << "!! ERROR: Trouble reading SPS file " << fname << "  !!" << endl;
   cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
   exit(EXIT_FAILURE);
  }
  int nshells;
  infile >> nshells;
  char pn[nshells];
  infile >> pn;
  for (int i=0;i<nshells;++i)
  {
    int dummy1,dummy2,dummy3,n,l,j2;
    infile >> dummy1 >> n >> l >> j2 >> dummy2 >> dummy3;
    int tz2 = pn[i]=='p' ? 1 : -1;
    for (int mj2=-j2;mj2<=j2;mj2+=2)
    {
       m_orbits.push_back( MschemeOrbit(n,2*l,j2,mj2,tz2) );
    }
  }


}






void NuBasis::PrintBasis()
{
  cout << "Number of Partitions: " << no_spart << endl;
  cout << "Total number with MJ=J: " << sum_nJT << endl;
  for (int isp=0;isp<no_spart;++isp)
  {
    cout << endl;
    cout << "isp = " << isp << endl;
    cout << "partition index: " << pindx[isp] << endl;
    cout << "Number with MJ=J: " << nJT[isp] << endl;
    cout << "J*2: " << j[isp] << endl;
    cout << "T*2: " << t[isp] << endl;
    cout << "Number of M-scheme configurations: " << ibf[isp] << endl;
    cout << "Maximum M*2 for this partition: " << max_cM2[isp] << endl;
    cout << "Maximum Tz*2 for this parition: " << max_cTz2[isp] << endl;
    cout << "Parition: ";
    for (auto c : partition[isp]) cout << int(c) << " ";
    cout << endl;
    for (int i=0;i<ibf[isp];i++) 
    {
        cout << "  i: " << i << endl;
        for (int g=0;g<gwords;++g)
        {
           cout << "  " << setw(10) << vec[isp][i][g] << "   " << bitset<8*gwords*sizeof(mvec_type)>(vec[isp][i][g]) << endl;
//           cout << "  " << setw(10) << vec[isp][i][g] << "   " << bitset<8*gwords*sizeof(mvec_type)>(TimeReverse(vec[isp][i])[g]) << endl;
//           vector<vector<mvec_type>> mvec_out;
           vector<key_type> mvec_out;
           vector<float> coef_st;
//           LoweringOperator( vec[isp][i], mvec_out, coef_st);
//           cout << "------------- Lowering operator --------------------" << endl;
//           for (size_t ilow=0;ilow<mvec_out.size();++ilow)
//           {
//             cout << "  " << setw(10) << coef_st[ilow] << " x " << bitset<8*gwords*sizeof(mvec_type)>(mvec_out[ilow][0]) << endl;
//           }
        }
        vector<int> occ;
//        for (int iword=0;iword<gwords;++iword)
//        {
         for (int ibit=0;ibit<sizeof(mvec_type)*8;++ibit)
         {
//           if ( (vec[isp][i][iword] >> ibit)&1 )  occ.push_back(ibit + iword*sizeof(mvec_type)*8);
//           if ( vec[isp][i] >> ibit).to_ulong() & 0x1L )  occ.push_back(ibit);
           if ( vec[isp][i][ibit] )  occ.push_back(ibit);
         }
//        }
        for (auto o : occ ) cout << o << " ";
        cout << endl;
        cout << endl;
    }
  }


}


void NuBasis::PrintSPS()
{

  for (size_t i=0;i<m_orbits.size();++i)
  {
    auto& mo = m_orbits[i];
    cout << i << ": " << mo.n << " " << mo.l2 << " " << mo.j2 << " " << mo.mj2 << " " << mo.tz2 << endl;
  }
}








/*
// Apply a lowering operator to an m-scheme state
// This is just a sum of lowering operators to each orbit, with the Pauli principle enforced.
//void NuBasis::LoweringOperator( vector<mvec_type>& mvec_in, vector<vector<mvec_type>>& mvec_out, vector<float>& coefs)
void NuBasis::LoweringOperator( key_type& mvec_in, vector<key_type>& mvec_out, vector<float>& coefs)
{
   for (size_t i_m=1;i_m<m_orbits.size();++i_m)  // don't bother starting with 0, since it's already in the lowest m_j state
   {
//      int iword = i_m/(sizeof(mvec_type)*8);
//      int ibit = i_m%(sizeof(mvec_type)*8);
//      if (ibit<1) continue;
//      if ( (not((mvec_in[iword]>>ibit)&0x1L)) or (mvec_in[iword]>>(ibit-1))&0x1L  ) continue; // Pauli principle
      if ( (not((mvec_in>>i_m)&0x1L)) or (mvec_in>>(i_m-1))&0x1L  ) continue; // Pauli principle
      
      int j2  = m_orbits[i_m].j2;
      int mj2 = m_orbits[i_m].mj2;
      if (mj2==-j2) continue;
//      vector<mvec_type> temp_mvec_out = mvec_in;
      key_type temp_mvec_out = mvec_in;
//      temp_mvec_out[iword] &= ~(0x1L << (ibit));
      temp_mvec_out &= ~(0x1L << (i_m));
//      temp_mvec_out[iword] |=  (0x1L << (ibit-1));
      temp_mvec_out |=  (0x1L << (i_m-1));
      mvec_out.push_back( temp_mvec_out );
      coefs.push_back( sqrt( j2*(j2+1)-mj2*(mj2-1) )*0.5 );
   }
   float norm = sqrt( inner_product(begin(coefs),end(coefs),begin(coefs),0.));
   if (abs(norm)>1e-9)
     for (size_t i=0;i<coefs.size();i++) coefs[i] /= norm;

}
*/


/*
// Apply a raising operator to an m-scheme state
// This is just a sum of raising operators to each orbit, with the Pauli principle enforced.
//void NuBasis::RaisingOperator( vector<mvec_type>& mvec_in, vector<vector<mvec_type>>& mvec_out, vector<float>& coefs)
void NuBasis::RaisingOperator( key_type& mvec_in, vector<key_type>& mvec_out, vector<float>& coefs)
{
   for (size_t i_m=0;i_m<m_orbits.size()-1;++i_m)  // don't bother with the last one, since it's already in the highest m_j state
   {
      int iword = i_m/(sizeof(mvec_type)*8);
      int ibit = i_m%(sizeof(mvec_type)*8);
      if (ibit==sizeof(mvec_type)*8) continue;
      if ( not(mvec_in[iword]>>ibit)&1 or (mvec_in[iword]>>(ibit-1))&1  ) continue;
      
      int j2 = m_orbits[i_m].j2;
      int mj2 = m_orbits[i_m].mj2;
      if (mj2==j2) continue;
//      vector<mvec_type> temp_mvec_out = mvec_in;
      key_type temp_mvec_out = mvec_in;
      temp_mvec_out[iword] &= ~(0x1L << (ibit));
      temp_mvec_out[iword] |=  (0x1L << (ibit+11));
      mvec_out.push_back( temp_mvec_out );
      coefs.push_back( sqrt( j2*(j2+1)-mj2*(mj2-1) )*0.5 );
   }

}
*/

/*
//vector<mvec_type> NuBasis::TimeReverse( vector<mvec_type>& mvec_in)
key_type NuBasis::TimeReverse( key_type& mvec_in)
{
   
   size_t i_m=0;
   vector<mvec_type> mvec_out = mvec_in;
   while (i_m<gwords*sizeof(mvec_type)*8) 
   {
      int iword = i_m /(sizeof(mvec_type)*8);
      size_t i_m_local = i_m % (sizeof(mvec_type)*8);
      int mj2 = m_orbits[i_m_local].mj2;
      int j2 = m_orbits[i_m_local].j2;
      if (mj2>0)
      {
        i_m += (j2-mj2)/2+1;
        continue;
      }
      size_t i_r = i_m - mj2; // time reversed m-state << TODO: This isn't right
      int iword_r = i_r /(sizeof(mvec_type)*8);
      i_r = i_r %(sizeof(mvec_type)*8);
      if ( (mvec_out[iword]>>i_m_local)&0x1L ^ (mvec_out[iword_r]>>i_r)&0x1L)
      {
        mvec_out[iword] ^= (0x1L << i_m_local);
        mvec_out[iword_r] ^= (0x1L << i_r);
      }
      i_m++;
   }
   return mvec_out;
}
*/



void NuBasis::Clear()
{
  no_spart = 0;
  sum_nJT = 0;
  pindx.clear();
  nJT.clear();
  j.clear();
  t.clear();
  ibf.clear();
  max_cM2.clear();
  max_cTz2.clear();
  partition.clear();
  vec.clear();
}

