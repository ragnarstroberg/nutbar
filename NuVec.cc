#include "NuVec.hh"
#include <cmath> // for abs(float)
#include <numeric> // for inner_product


NuVec::NuVec()
{}

NuVec::NuVec(int J2)
: J2(J2)
{}


void NuVec::ReadFile(string fname)
{

 int nwords=0;

 ifstream infile(fname, ios::binary );
  if (not infile.good() )
  {
    cout << "NuVec::ReadFile -- Trouble opening file " << fname << endl;
    return ;
  }


 infile.read((char*)&nwords, DELIMITER_LENGTH);
 if (nwords != sizeof(no_state)+sizeof(no_level))
 {
   cout << "NuVec::ReadFile -- Trouble in file " << fname
        << " :  Want to read " << sizeof(no_state)+sizeof(no_level)
        << " words, but record length is " << nwords << endl;
   return;
 }
 infile.read((char*)&no_state, sizeof(no_state));
 infile.read((char*)&no_level, sizeof(no_level));
 infile.read((char*)&nwords, DELIMITER_LENGTH);

 alpha.resize(no_level);
 coefT.resize(no_level,vector<float>(no_state));

// cout << "no_state: " << no_state << endl;
// cout << "no_level: " << no_level << endl;
 for (int ilevel=0;ilevel<no_level;++ilevel)
 {
   infile.read((char*)&nwords, DELIMITER_LENGTH);
   if (nwords != ( sizeof(alpha[ilevel]) + no_state*sizeof(coefT[ilevel][0]) ) )
   {
      cout << "NuVec::ReadFile -- Trouble in file " << fname
           << " :  Want to read " <<  sizeof(alpha[ilevel]) + no_state*sizeof(coefT[ilevel][0])
           << " words, but record length is " << nwords << endl;
   }
   infile.read((char*)&(alpha[ilevel]), sizeof(alpha[ilevel]));
   infile.read((char*)&(coefT[ilevel][0]), no_state*sizeof(coefT[ilevel][0]));
  
   infile.read((char*)&nwords, DELIMITER_LENGTH);
 }

}


void NuVec::PrintVectors()
{
 for (int ilevel=0;ilevel<no_level;++ilevel)
 {
     cout << ilevel << ": ( ";
     for (auto c : coefT[ilevel]) cout << fixed << setprecision(6) << setw(9) << c << ", ";
     cout << " ) " << endl;
 }

}
void NuVec::PrintDetailedVectors()
{

 for (int ilevel=0;ilevel<no_level;++ilevel)
 {
     cout << "------------ " << ilevel << " ---------------" << endl;
     cout << "alpha: " << fixed << setw(12) << setprecision(8) << alpha[ilevel] << endl;
     float sumcoef2 = 0;
     for (int istate=0;istate<no_state;++istate)
     {
      sumcoef2 += coefT[ilevel][istate]*coefT[ilevel][istate];
      cout << istate << " :  " << setw(12) << setprecision(8) << coefT[ilevel][istate] 
           << "    " << setw(12) << setprecision(8) << coefT[ilevel][istate]*coefT[ilevel][istate]
           << "    " << setw(12) << setprecision(8) << sumcoef2
           << endl;
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
      
     if ((ilevel!=jlevel and (abs(sumcoef2)>1.e-4)) or (ilevel==jlevel and (abs(sumcoef2-1.0)>1.e-4)) )
     {
       cout << "Trouble in CheckOrthoNormal. < state " << ilevel << " | state " << jlevel << " >  =  " << sumcoef2 << endl;
       return false;
     }
   }
  }
  return true; // passed
}
