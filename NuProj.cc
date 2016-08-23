#include "NuProj.hh"


NuProj::NuProj()
 : nptt(0),ngood(0)
{}


void NuProj::ReadFile(string fname)
{

 int nwords;

 ifstream infile(fname, ios::binary );

 infile.seekg(DELIMITER_LENGTH,ios::cur);
 infile.read((char*)&nptt, sizeof(nptt));
 infile.seekg(DELIMITER_LENGTH,ios::cur);

 for (int isp=0;isp<nptt;isp++)
 {

    int32_t pindx_in=-1, ngood_in=-1, j_in=-1, t_in=-1, dim_in=-1;
    float x_in=-1;
    vector<part_type> partition;

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
  cout << "ngood = " << ngood << endl;
  cout << endl;

  for (int igood=0;igood<ngood;++igood)
  {
    cout << "basis state: " << igood << endl;
    cout << "pindx: " << pindx[igood] << endl;
    cout << "j: " << j[igood] << endl;
    cout << "t: " << t[igood] << endl;
    cout << "dim: " << dim[igood] << endl;
    cout << "x: " << x[igood] << endl;
    cout << "coef: ";
    for (auto c : coef_st[igood]) cout << c << " ";
    cout << endl;
    cout << endl;
  }
  cout << endl;


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
