
#include <cmath> // for std::abs(float)
#include <numeric> // for inner_product
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "NuVec.hh"


const int NuVec::delimiter_length = 4;

NuVec::NuVec()
{}

NuVec::NuVec(int J2)
: J2(J2)
{}

void NuVec::ReadFile(std::string fname)
{

 int nwords=0;

 std::ifstream infile(fname, std::ios::binary );
 
  if (not infile.good() )
  {
    std::cout << "NuVec::ReadFile -- Trouble opening file " << fname << std::endl;
    return ;
  }


// infile.read((char*)&nwords, delimiter_length);
 infile.read((char*)&nwords, delimiter_length);
 if (nwords != sizeof(no_state)+sizeof(no_level))
 {
   std::cout << "NuVec::ReadFile -- Trouble in file " << fname
        << " :  Want to read " << sizeof(no_state)+sizeof(no_level)
        << " words, but record length is " << nwords << std::endl;
   return;
 }
 infile.read((char*)&no_state, sizeof(no_state));
 infile.read((char*)&no_level, sizeof(no_level));
 infile.read((char*)&nwords, delimiter_length);

 alpha.resize(no_level);
 coefT.resize(no_level,std::vector<float>(no_state));

// std::cout << "no_state: " << no_state << std::endl;
// std::cout << "no_level: " << no_level << std::endl;
 size_t blocksize = 2 * delimiter_length + no_state*sizeof(alpha[0]) + no_state*sizeof(coefT[0][0]);

 for (int ilevel=0;ilevel<no_level;++ilevel)
 {
   infile.read((char*)&nwords, delimiter_length);
   if (nwords != ( sizeof(alpha[ilevel]) + no_state*sizeof(coefT[ilevel][0]) ) )
   {
      std::cout << "NuVec::ReadFile -- Trouble in file " << fname
           << " :  Want to read " <<  sizeof(alpha[ilevel]) + no_state*sizeof(coefT[ilevel][0])
           << " words, but record length is " << nwords << std::endl;
   }
   infile.read((char*)&(alpha[ilevel]), sizeof(alpha[ilevel]));
   infile.read((char*)&(coefT[ilevel][0]), no_state*sizeof(coefT[ilevel][0]));
  
   infile.read((char*)&nwords, delimiter_length);
 }

}


void NuVec::PrintVectors()
{
 for (int ilevel=0;ilevel<no_level;++ilevel)
 {
     std::cout << ilevel << ": ( ";
     for (auto c : coefT[ilevel]) std::cout << std::fixed << std::setprecision(6) << std::setw(9) << c << ", ";
     std::cout << " ) " << std::endl;
 }

}
void NuVec::PrintDetailedVectors()
{

 for (int ilevel=0;ilevel<no_level;++ilevel)
 {
     std::cout << "------------ " << ilevel << " ---------------" << std::endl;
     std::cout << "alpha: " << std::fixed << std::setw(12) << std::setprecision(8) << alpha[ilevel] << std::endl;
     float sumcoef2 = 0;
     for (int istate=0;istate<no_state;++istate)
     {
      sumcoef2 += coefT[ilevel][istate]*coefT[ilevel][istate];
      std::cout << istate << " :  " << std::setw(12) << std::setprecision(8) << coefT[ilevel][istate] 
           << "    " << std::setw(12) << std::setprecision(8) << coefT[ilevel][istate]*coefT[ilevel][istate]
           << "    " << std::setw(12) << std::setprecision(8) << sumcoef2
           << std::endl;
     }
 }

}


bool NuVec::CheckOrthoNormal()
{
  for (int ilevel=0;ilevel<no_level;++ilevel)
  {
   for (int jlevel=ilevel;jlevel<no_level;++jlevel)
   {
     float sumcoef2 = inner_product(begin(coefT[ilevel]),end(coefT[ilevel]),begin(coefT[jlevel]),0.0);
      
     if ((ilevel!=jlevel and (std::abs(sumcoef2)>1.e-4)) or (ilevel==jlevel and (std::abs(sumcoef2-1.0)>1.e-4)) )
     {
       std::cout << "Trouble in CheckOrthoNormal. < state " << ilevel << " | state " << jlevel << " >  =  " << sumcoef2 << std::endl;
       return false;
     }
   }
  }
  return true; // passed
}
