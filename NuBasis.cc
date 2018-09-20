
#include <cmath>
#include <numeric>
#include <vector>
#include <string>
#include <bitset>
#include <iostream>
#include <iomanip>
#include "NuBasis.hh"


const int NuBasis::gwords = 1; // this is good if the modelspace has <=64 m-states. Otherwise, we may need to do some tinkering.
//const int NuBasis::delimiter_length = 4;


NuBasis::NuBasis()
{}


void NuBasis::ReadFile(std::string fname)
{
  Clear();
  std::ifstream infile(fname, std::ios::binary );
  if (not infile.good() )
  {
   std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
   std::cout << "!! ERROR: Trouble in NuBasis::ReadFile --" << fname << "  !!" << std::endl;
   std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    exit(EXIT_FAILURE) ;
  }

 int nwords;
 sum_nJT=0;

 infile.read((char*)&nwords, sizeof(nwords));
 infile.read((char*)&no_spart, sizeof(no_spart));
 infile.read((char*)&nwords, sizeof(nwords));
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

    sum_nJT += nJT[isp];
   
    infile.read((char*)&nwords, sizeof(nwords));
    partition[isp].resize(nwords/sizeof(part_type));
    infile.read((char*)&partition[isp][0], nwords);
    infile.read((char*)&nwords, sizeof(nwords));
    
    if (ibf[isp]>0)
    {    
    infile.read((char*)&nwords, sizeof(nwords));
      vec[isp].resize(ibf[isp]);
      for (int i=0;i<ibf[isp];i++) 
      {
        infile.read((char*)&vec[isp][i], gwords*sizeof(mvec_type));
      }
      infile.read((char*)&nwords, sizeof(nwords));
    }
 }

}






void NuBasis::ReadSPS(std::string fname)
{

  std::ifstream infile( fname );
  if ( not infile.good() )
  {
   std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
   std::cout << "!! ERROR: Trouble reading SPS file " << fname << "  !!" << std::endl;
   std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
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
  std::cout << "Number of Partitions: " << no_spart << std::endl;
  std::cout << "Total number with MJ=J: " << sum_nJT << std::endl;
  for (int isp=0;isp<no_spart;++isp)
  {
    std::cout << std::endl;
    std::cout << "isp = " << isp << std::endl;
    std::cout << "partition index: " << pindx[isp] << std::endl;
    std::cout << "Number with MJ=J: " << nJT[isp] << std::endl;
    std::cout << "J*2: " << j[isp] << std::endl;
    std::cout << "T*2: " << t[isp] << std::endl;
    std::cout << "Number of M-scheme configurations: " << ibf[isp] << std::endl;
    std::cout << "Maximum M*2 for this partition: " << max_cM2[isp] << std::endl;
    std::cout << "Maximum Tz*2 for this parition: " << max_cTz2[isp] << std::endl;
    std::cout << "Parition: ";
    for (auto c : partition[isp]) std::cout << int(c) << " ";
    std::cout << std::endl;
    for (int i=0;i<ibf[isp];i++) 
    {
        std::cout << "  i: " << i << std::endl;
        for (int g=0;g<gwords;++g)
        {
           std::cout << "  " << std::setw(10) << vec[isp][i][g] << "   " << std::bitset<8*gwords*sizeof(mvec_type)>(vec[isp][i][g]) << std::endl;
//           std::cout << "  " << std::setw(10) << vec[isp][i][g] << "   " << std::bitset<8*gwords*sizeof(mvec_type)>(TimeReverse(vec[isp][i])[g]) << std::endl;
//           std::vector<std::vector<mvec_type>> mvec_out;
           std::vector<key_type> mvec_out;
           std::vector<float> coef_st;
//           LoweringOperator( vec[isp][i], mvec_out, coef_st);
//           std::cout << "------------- Lowering operator --------------------" << std::endl;
//           for (size_t ilow=0;ilow<mvec_out.size();++ilow)
//           {
//             std::cout << "  " << std::setw(10) << coef_st[ilow] << " x " << std::bitset<8*gwords*sizeof(mvec_type)>(mvec_out[ilow][0]) << std::endl;
//           }
        }
        std::vector<int> occ;
//        for (int iword=0;iword<gwords;++iword)
//        {
         for (int ibit=0;ibit<sizeof(mvec_type)*8;++ibit)
         {
//           if ( (vec[isp][i][iword] >> ibit)&1 )  occ.push_back(ibit + iword*sizeof(mvec_type)*8);
//           if ( vec[isp][i] >> ibit).to_ulong() & 0x1L )  occ.push_back(ibit);
           if ( vec[isp][i][ibit] )  occ.push_back(ibit);
         }
//        }
        for (auto o : occ ) std::cout << o << " ";
        std::cout << std::endl;
        std::cout << std::endl;
    }
  }


}


void NuBasis::PrintSPS()
{

  for (size_t i=0;i<m_orbits.size();++i)
  {
    auto& mo = m_orbits[i];
    std::cout << i << ": " << mo.n << " " << mo.l2 << " " << mo.j2 << " " << mo.mj2 << " " << mo.tz2 << std::endl;
  }
}




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

