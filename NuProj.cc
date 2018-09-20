
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdint>

#include "NuProj.hh"


const int NuProj::delimiter_length = 4;


NuProj::NuProj()
 : nptt(0),ngood(0)
{}


void NuProj::ReadFile(std::string fname)
{

 int nwords;
 Clear();

 std::ifstream infile(fname, std::ios::binary );

 infile.seekg(delimiter_length,std::ios::cur);
 infile.read((char*)&nptt, sizeof(nptt));
 infile.seekg(delimiter_length,std::ios::cur);

 for (int isp=0;isp<nptt;isp++)
 {

    int32_t pindx_in=-1, ngood_in=-1, j_in=-1, t_in=-1, dim_in=-1;
    float x_in=-1;
    std::vector<part_type> partition;

    infile.read((char*)&nwords, sizeof(nwords));
    infile.read((char*)&pindx_in, sizeof(pindx_in));
    infile.read((char*)&ngood_in, sizeof(ngood_in));
    infile.read((char*)&j_in, sizeof(j_in));
    infile.read((char*)&t_in, sizeof(t_in));
    infile.read((char*)&dim_in, sizeof(dim_in));
    infile.read((char*)&x_in, sizeof(x_in));
    infile.read((char*)&nwords, sizeof(nwords));
   
    infile.read((char*)&nwords, sizeof(nwords));
    partition.resize(nwords/sizeof(part_type));
    infile.read((char*)&partition[0], nwords);
    infile.read((char*)&nwords, sizeof(nwords));
   
    if (dim_in>0 and ngood_in>0)
    {
      coef_st.resize(ngood+ngood_in);
      for (int ng=0;ng<ngood_in;++ng)
      {
        pindx.push_back(pindx_in);
        j.push_back(j_in);
        t.push_back(t_in);
        dim.push_back(dim_in);
        x.push_back(x_in);
        coef_st[ngood+ng].resize(dim_in);
        infile.read((char*)&nwords, sizeof(nwords));
        infile.read((char*)&coef_st[ngood+ng][0], coef_st[ngood+ng].size()*sizeof(coef_st[ngood+ng][0]));
        infile.read((char*)&nwords, sizeof(nwords));

      }
    }
    ngood += ngood_in;
 }

}



void NuProj::PrintProj()
{
  std::cout << "ngood = " << ngood << std::endl;
  std::cout << std::endl;

  for (int igood=0;igood<ngood;++igood)
  {
    std::cout << "basis state: " << igood << std::endl;
    std::cout << "pindx: " << pindx[igood] << std::endl;
    std::cout << "j: " << j[igood] << std::endl;
    std::cout << "t: " << t[igood] << std::endl;
    std::cout << "dim: " << dim[igood] << std::endl;
    std::cout << "x: " << x[igood] << std::endl;
    std::cout << "coef: ";
    for (auto c : coef_st[igood]) std::cout << c << " ";
    std::cout << std::endl;
    std::cout << std::endl;
  }
  std::cout << std::endl;


}



void NuProj::Clear()
{
  nptt = 0;
  ngood = 0;
  no_spart.clear();
  pindx.clear();
  j.clear();
  t.clear();
  dim.clear();
  max_cM2.clear();
  max_cTz2.clear();

}
