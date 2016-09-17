#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <armadillo>

#include "TransitionDensity.hh"

#define SQRT2 1.4142135623730950488

using namespace std;


vector<char> TransitionDensity::an_code = {  '0','1','2','3','4','5','6','7','8','9',
                                             '_','-','c','d','e','f','g','h','i','j',  
                                             'k','l','m','n','o','p','q','r','s','t',  
                                             'u','v','w','x','y','z','A','B','C','D',  
                                             'E','F','G','H','I','J','K','L','M','N',  
                                             'O','P','Q','R','S','T','U','V','W','X',  
                                             '0','1','2','3','4','5','6','7','8','9',  
                                             '_','-','c'};

vector<string> TransitionDensity::periodic_table = {
   "N","h_","he",
       "li","be","b_","c_","n_","o_","f_","ne",
       "na","mg","al","si","p_","s_","cl","ar",
       "k_","ca","sc","ti","v_","cr","mn","fe","co","ni","cu","zn","ga","ge","as","se","br","kr",
       "rb","sr","y_","zr","nb","mo","tc","ru","rh","pd","ag","cd","in","sn","sb","te","i_","xe",
       "cs","ba",
         "la","ce","pr","nd","pm","sm","eu","gd","tb","dy","ho","er","tm","yb",
                 "lu","hf","ta","w_","re","os","ir","pr","au","hg","tl","pb","bi","po","at","rn",
       "fr","ra",
         "ac","th","u_","np","pu","am","cm","bk","cf","es","fm","md","no","lr",
                 "rf","db","sg","bh","hs","mt","ds","rg","cn"};


TransitionDensity::TransitionDensity()
: total_number_levels(0)
{}

TransitionDensity::TransitionDensity(vector<int> jlist)
:total_number_levels(0), Jlist(jlist)
{}




void TransitionDensity::ReadInputInteractive()
{

  cout << endl;
  cout << "Name of of NuShellFiles, (e.g. ne200): ";
  cout.flush();
  cin >> basename;

  cout << "Name of sps file, (e.g. sdpn.sps): ";
  cout.flush();
  cin >> sps_file_name;

  cout << "sps_file_name = " << sps_file_name << endl;



  istringstream iss;
  string line;

  cout << "J*2 values of eigenvectors (separated by space): ";
  cout.flush();
  getline(cin,line);
  getline(cin,line);
  iss.clear();
  iss.str(line);
  int j;
  while( iss >> j ) Jlist.push_back(j);

//  cout << "J values: ";
//  for (auto j : Jlist) cout << j << " ";
//  cout << endl;


  cout << "Number of eigenstates for each J (separated by space, or single value if same for all J): ";
  cout.flush();
  getline(cin,line);
  iss.clear();
  iss.str(line);
  vector<int> spj;
  while( iss >> j ) spj.push_back(j);
  for (size_t i=0;i<Jlist.size();++i)
  {
    if (spj.size() < Jlist.size()) max_states_per_J[ Jlist[i] ] = spj[0];
    else                           max_states_per_J[ Jlist[i] ] = spj[i];
  }


  ofstream outfile("NuShelltoMBPT.input");
  outfile << basename << endl;
  outfile << sps_file_name << endl;
  for (auto j : Jlist) outfile << j << " ";
  outfile << endl;
  for (auto j : spj) outfile << j << " ";
  outfile << endl;


}



void TransitionDensity::ReadInputFromFile(string filename)
{
  ifstream infile(filename);
  istringstream iss;

  const int BUFFERSIZE = 500;
  char line[BUFFERSIZE];


  // basename for *.nba *.prj and *.xvc files
  infile.getline(line, BUFFERSIZE);
  iss.str(line);
  iss >> basename;

  cout << "basename = " << basename << endl;

  // basename for *.nba *.prj and *.xvc files
  infile.getline(line, BUFFERSIZE);
  iss.clear();
  iss.str(line);
  iss >> sps_file_name;

  cout << "sps_file_name = " << sps_file_name << endl;


  // list of 2*J values
  infile.getline(line, BUFFERSIZE);
  iss.clear();
  iss.str(line);
  int j;
  while( iss >> j ) Jlist.push_back(j);

  cout << "J values: ";
  for (auto j : Jlist) cout << j << " ";
  cout << endl;

  // number of states to treat for each J
  infile.getline(line, BUFFERSIZE);
  iss.clear();
  iss.str(line);
  vector<int> spj;
  while( iss >> j ) spj.push_back(j);
  for (size_t i=0;i<Jlist.size();++i)
  {
    if (spj.size() < Jlist.size()) max_states_per_J[ Jlist[i] ] = spj[0];
    else                           max_states_per_J[ Jlist[i] ] = spj[i];
  }




}





//void TransitionDensity::ReadFiles( string spsfile, string abfile_base )
void TransitionDensity::ReadFiles( )
{

  ifstream testread; // for checking if files exist
  ostringstream ostr;
  vector<string> Afiles,Bfiles;
  GetAZFromFileName(); // Make a default guess of the core
  ReadSPfile(); // If the sp file is there, use that to get the core
  cout << "Acore,Zcore = " << Acore << " " << Zcore << endl;

  int nvalence_protons = Z - Zcore;
  int nvalence_neutrons = A-Z - (Acore-Zcore);
  for ( auto j : Jlist ) cout << j << " ";
  cout << endl;
  MJtot = (*min_element(begin(Jlist),end(Jlist)))%2;

  for (auto j : Jlist )
  {
    if (j%2 != A%2)
    {
      cout << "Warning A=" << A << " and J*2 = " << j << ".  This isn't good" << endl;
    }
  }
  
  
  // Find all the prj and nba files
  int iJ;
  for (int icode : { nvalence_protons, -nvalence_protons, nvalence_neutrons, -nvalence_neutrons  } )
  {
    iJ=0;
    while (true)
    {
      ostr.str("");
      ostr.str().clear();
      ostr << basename << "a" << iJ << "0" << an_code[36 + icode/2];
      testread.open(ostr.str()+".nba");
      if ( testread.fail() ) break;
      if ( find(Afiles.begin(),Afiles.end(), ostr.str() ) == Afiles.end())
        Afiles.push_back( ostr.str() );
      testread.close();
      iJ++;
    }
    iJ=0;
    while (true)
    {
      ostr.str("");
      ostr.str().clear();
      ostr << basename << "b" << iJ << "0" << an_code[36 + icode/2];
      testread.open(ostr.str()+".nba");
      if ( testread.fail() ) break;
      if ( find(Bfiles.begin(),Bfiles.end(), ostr.str() ) == Bfiles.end())
      Bfiles.push_back( ostr.str() );
      testread.close();
      iJ++;
    }
  }

  for (int Jtot : Jlist )
  {
    jbasis_list.emplace_back( JBasis( sps_file_name, Afiles, Bfiles, Jtot, MJtot));
  
  
  // Guess the name of the xvc file
    ostr.str("");
    ostr.clear();
    ostr << basename << an_code[Jtot/2] << an_code[nvalence_protons] << an_code[nvalence_neutrons] << ".xvc";
    string vecfile = ostr.str();
    testread.open(vecfile);
    if ( testread.fail() )
    {
      ostr.str("");
      ostr.clear();
      ostr << basename << an_code[Jtot/2] << an_code[nvalence_neutrons] << an_code[nvalence_protons] << ".xvc";
      vecfile = ostr.str();
      testread.close();
      testread.open(vecfile);
      if ( testread.fail() )
      {
        cout << "ERROR! I cant figure out what the *.xvc file should be. Exiting." << endl;
        cout << "( " << ostr.str() << " ) didnt work. abfile_base = " << basename << endl;
        return ;
      }
    }
    testread.close();
  
    nuvec_list.emplace_back( NuVec(Jtot) );
    nuvec_list.back().ReadFile(vecfile);
  
  }

  m_orbits = jbasis_list[0].nubasis_a.m_orbits;
  Nshell=0;
  for (auto morbit : m_orbits)
  {
    Nshell = max(Nshell, 2*morbit.n + morbit.l2/2);
  }
  Nshell++;

  for (size_t ivec=0; ivec<nuvec_list.size(); ++ivec)
  {
   int imax = nuvec_list[ivec].no_level;
   if ( max_states_per_J.find(nuvec_list[ivec].J2) != max_states_per_J.end() ) imax = min(imax,max_states_per_J[nuvec_list[ivec].J2]);
    blank_vector.push_back(vector<float>(imax, 0.));
    total_number_levels += blank_vector.back().size();
  }

}


void TransitionDensity::CalculateMschemeAmplitudes()
{
  for (size_t ivec=0; ivec<nuvec_list.size(); ++ivec)
  {
   cout << "ivec = " << ivec << endl;
   auto& nuvec = nuvec_list[ivec];
   auto& jbasis = jbasis_list[ivec];
   cout << "no_state = " << nuvec.no_state << endl;
   cout << "no_level = " << nuvec.no_level << endl;
  
   // loop over J-coupled basis states
   for (int istate=0;istate<nuvec.no_state;++istate)
   {

     int imax = nuvec.no_level;
     if ( max_states_per_J.find(nuvec.J2) != max_states_per_J.end() ) imax = min(imax,max_states_per_J[nuvec.J2]);
     vector<float> level_coefs(imax,0.0);
     for (int ilevel=0;ilevel<imax;++ilevel)
     {
       level_coefs[ilevel] = nuvec.coefT[ilevel][istate];
     }
  
     if (istate>=jbasis.basis_states.size())
     {
       cout << "ERROR istate = " << istate << ",  basis_states.size() = " << jbasis.basis_states.size() << endl;
       return;
     }
     JMState& jmst = jbasis.basis_states[istate];
     for (auto& it_mstate : jmst.m_coefs)
     {
       auto& key = it_mstate.first;
       float& m_coef = it_mstate.second;
       if( amplitudes.find(key) == amplitudes.end() ) amplitudes[key] = blank_vector;
       amplitudes[key][ivec] += m_coef * level_coefs;
     }
   }
  }

}


double TransitionDensity::OBTD(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a, int m_index_b, int Lambda2 )
{

  int J2i = Jlist[J_index_i];
  int J2f = Jlist[J_index_f];
  
  int j2_a = m_orbits[m_index_a].j2;
  int j2_b = m_orbits[m_index_b].j2;

  if ((J2i+Lambda2 < J2f) or (abs(J2i-Lambda2)>J2f)) return 0;
  
  double clebsch_fi = CG(J2i,0,Lambda2,0,J2f,0);
  if (abs(clebsch_fi)<1e-9)
  {
     cout << "Warning got zero Clebsch while inverting the Wigner-Eckart theorem" << endl;
     return 0;
  }
  
  // find m-scheme orbits so that m_a = m_b, which will work for mu=0
  double obd = 0;
  while ( m_orbits[m_index_a].mj2 < min(j2_a, m_orbits[m_index_b].mj2) ) m_index_a++;
  while ( m_orbits[m_index_b].mj2 < min(j2_b, m_orbits[m_index_a].mj2) ) m_index_b++;
  if (m_orbits[m_index_a].j2 != j2_a) return 0;
  if (m_orbits[m_index_b].j2 != j2_b) return 0;
  m_index_a--;
  m_index_b--;

  for ( int im=-min(j2_a,j2_b);im<=min(j2_a,j2_b);im+=2)
  {
    m_index_a++;
    m_index_b++;

    // convention: tilded destruction operator b~(m) = (-1)**(jb + mb) b(-m)
    //                                         b(m)  = (-1)**(jb -mb) b~(-m)
    int phase_b = (1-(j2_b-im)%4);
    double clebsch = CG(j2_a,im,j2_b,-im,Lambda2,0) ;

    for ( auto& it_amp : amplitudes )
    {
      auto& key = it_amp.first;
      if ( not( (key[0] >> m_index_a)&0x1) ) continue;
      if ( m_index_a != m_index_b and   ( (key[0] >> m_index_b)&0x1) ) continue;
      double amp_i = it_amp.second[J_index_i][eigvec_i];
      if (abs(amp_i)<1e-7) continue;
      auto new_key = key;
      new_key[0] &= ~( 0x1 << (m_index_a));
      new_key[0] |=  ( 0x1 << (m_index_b));
      if (amplitudes.find(new_key) == amplitudes.end() ) continue;
      int phase_ladder = 0;
      for (int iphase=min(m_index_a,m_index_b)+1;iphase<max(m_index_a,m_index_b);++iphase) phase_ladder += ( key[0] >>iphase )&0x1;
      phase_ladder = 1-2*(phase_ladder%2);
      double amp_f = amplitudes[new_key][J_index_f][eigvec_f];
      obd += clebsch * amp_i * amp_f * phase_ladder * phase_b;
    }
  }
  
  obd *= sqrt((J2f+1.)/(Lambda2+1)) / clebsch_fi;

  return obd;

}



double TransitionDensity::TBTD(int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int m_index_a, int m_index_b, int m_index_c, int m_index_d, int J2ab, int J2cd, int Lambda2 )
{

  int J2i = Jlist[J_index_i];
  int J2f = Jlist[J_index_f];
  int mu = 0;
  int Mi = 0;
  int Mf = 0;
  
  int j2_a = m_orbits[m_index_a].j2;
  int j2_b = m_orbits[m_index_b].j2;
  int j2_c = m_orbits[m_index_c].j2;
  int j2_d = m_orbits[m_index_d].j2;

  // start out with the maximally-projected m state
  m_index_a += (j2_a - m_orbits[m_index_a].mj2)/2;
  m_index_b += (j2_b - m_orbits[m_index_b].mj2)/2;
  m_index_c += (j2_c - m_orbits[m_index_c].mj2)/2;
  m_index_d += (j2_d - m_orbits[m_index_d].mj2)/2;


  // check some triangle conditions
  if ((J2f+Lambda2 < J2i) or (abs(J2f-Lambda2)>J2i)) return 0;
  if ((j2_a+j2_b<J2ab) or (abs(j2_a-j2_b)>J2ab)) return 0;
  if ((j2_c+j2_d<J2cd) or (abs(j2_c-j2_d)>J2cd)) return 0;
  if ((J2ab+J2cd < Lambda2) or (abs(J2ab-J2cd)>Lambda2)) return 0;
  

  double clebsch_fi = CG(J2i,Mi,Lambda2,mu,J2f,Mf);
  if (abs(clebsch_fi)<1e-9)
  {
     cout << "Warning got zero Clebsch while inverting the Wigner-Eckart theorem" << endl;
     return 0;
  }

  double tbd = 0;
  int Mab_max = min(J2ab,J2cd);

  // Restricted sum a<=b c<=d gives 1/[(1+delta_ab)(1+delta_cd)], while
  // using normalized TBMEs gives sqrt[ (1+delta_ab)(1+delta_cd) ].
  // The additional factor of 2 comes from being able to limit ma < mb or mc < md
  double norm = 1;
  if (m_index_a==m_index_b) norm *= SQRT2;  
  if (m_index_c==m_index_d) norm *= SQRT2;  

  for (int Mab=-Mab_max;Mab<=Mab_max;Mab+=2)
  {
    int Mcd = mu - Mab;
    int phasecd = (1- abs(J2cd - Mcd)%4); // phase from getting rid of the tildes
    double clebsch_abcd = CG(J2ab,Mab,J2cd,Mcd,Lambda2,mu);
    if ( abs(clebsch_abcd)<1e-7) continue;
    int ma_min = max(-j2_a, Mab-j2_b);
    int ma_max = min(j2_a, Mab+j2_b);
    if (m_index_a==m_index_b) ma_max = min(ma_max, Mab/2);

    for (int ma=ma_min;ma<=ma_max;ma+=2)
    {
      int mb=Mab-ma;
      double clebsch_ab = CG(j2_a,ma, j2_b,mb, J2ab,Mab);
      if ( abs(clebsch_ab)<1e-7) continue;
      int ia = m_index_a - ( j2_a -ma )/2;
      int ib = m_index_b - ( j2_b -mb )/2;
      if (ib==ia) continue;
      int mc_min = max(-j2_c, Mab-j2_d);
      int mc_max = min(j2_c, Mab+j2_d);
      if (m_index_c == m_index_d) mc_max = min(j2_c, Mab/2);

      for (int mc=mc_min;mc<=mc_max;mc+=2)
      {
        int md = -Mcd-mc;
        double clebsch_cd = CG(j2_c,mc, j2_d, md, J2cd,-Mcd); // Mcd = -Mab, and another (-) comes from getting rid of the tildes
        if ( abs(clebsch_cd)<1e-7) continue;
        int ic = m_index_c - ( j2_c -mc )/2;
        int id = m_index_d - ( j2_d -md )/2;
        if (ic==id) continue;

        uint64_t mask_ab = ((0x1 << ia ) + (0x1 << ib) );
        uint64_t mask_cd = ((0x1 << ic ) + (0x1 << id) );

        for (auto& it_amp : amplitudes )
        {
          double amp_i = it_amp.second[J_index_i][eigvec_i];
          if (abs(amp_i)<1e-7) continue;
          auto& key = it_amp.first;

          if (not(( (key[0] >> ic) & (key[0] >> id) )&0x1) ) continue;
          auto new_key = key;
          new_key[0] &= (~ mask_cd); // remove particles from c and then from d  (d-c-)
          if ( new_key[0] & mask_ab) continue;
          new_key[0] |= mask_ab; // add particles to b and then to a  (a+b+)

          auto amp_newkey = amplitudes.find(new_key);
          if (amp_newkey == amplitudes.end() 
            or amp_newkey->second.size() < J_index_f
            or amp_newkey->second[J_index_f].size() < eigvec_f) continue;

          double amp_f = amp_newkey->second[J_index_f][eigvec_f];

          // pick up a phase from commuting the ladder operators
          int phase_ladder = (ia>ib xor id<ic) ? -1 : 1;
          for (int iphase = min(ic,id)+1;iphase<max(ic,id);++iphase)  if(( key[0] >>iphase )&0x1) phase_ladder *=-1;
          for (int iphase = min(ia,ib)+1;iphase<max(ia,ib);++iphase)  if(( new_key[0] >>iphase )&0x1) phase_ladder *=-1;

          tbd += amp_i * amp_f * clebsch_abcd * clebsch_ab * clebsch_cd * phasecd * phase_ladder;

        }
      }
    }
  }

  tbd *= sqrt((J2f+1.)/(Lambda2+1.)) / clebsch_fi * norm;
  return tbd;

}



arma::mat TransitionDensity::CalcOBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2)
{
  vector<int> jorbits;
  for (size_t i=0;i<m_orbits.size();++i )
  {
    if (m_orbits[i].mj2 == -m_orbits[i].j2) jorbits.push_back(i);
  }

  cout << "Initial state J2 = " << nuvec_list[J_index_i].J2 << "   E = " << nuvec_list[J_index_i].alpha[eigvec_i] << endl;
  cout << "  Final state J2 = " << nuvec_list[J_index_f].J2 << "   E = " << nuvec_list[J_index_f].alpha[eigvec_f] << endl;

  size_t njorb = jorbits.size();
  arma::mat obtd(njorb,njorb,arma::fill::zeros);
  for (size_t i=0; i<njorb; ++i)
  {
    int j2i = m_orbits[jorbits[i]].j2;
    for (size_t j=0; j<njorb; ++j)
    {
      obtd(i,j) = OBTD( J_index_i, eigvec_i, J_index_f, eigvec_f, jorbits[i], jorbits[j], Lambda2);
    }
  }

  return obtd;

}





arma::mat TransitionDensity::CalcTBTD( int J_index_i, int eigvec_i, int J_index_f, int eigvec_f, int Lambda2)
{
  // make a list of m-scheme indices for the beginning of each j-shell
  vector<int> jorbits;
  for (size_t i=0;i<m_orbits.size();++i )
  {
    if (m_orbits[i].mj2 == -m_orbits[i].j2) jorbits.push_back(i);
  }

  // generate all the two body states that are needed
  vector<int> ket_a, ket_b, ket_J;
  for (size_t a=0;a<jorbits.size();++a)
  {
    int ja = m_orbits[jorbits[a]].j2;
    for (size_t b=a; b<jorbits.size();++b)
    {      
      int jb = m_orbits[jorbits[b]].j2;
      int Jmin = abs(ja-jb);
      int Jmax = ja+jb;
      for (int J2=Jmin;J2<=Jmax;J2+=2)
      {
        if (a==b and (J2%4)>0) continue;
        ket_a.push_back(a);
        ket_b.push_back(b);
        ket_J.push_back(J2);
      }
    }
  }

  arma::mat tbtd(ket_J.size(), ket_J.size(), arma::fill::zeros);
  #pragma omp parallel for
  for (size_t ibra=0;ibra<ket_J.size();++ibra)
  {
    int a = ket_a[ibra];
    int b = ket_b[ibra];
    int J2ab = ket_J[ibra];
    size_t iket_min = ((J_index_i == J_index_f) and (eigvec_i==eigvec_f)) ? ibra : 0;
    for (size_t iket=iket_min;iket<ket_J.size();++iket)
    {
      int c = ket_a[iket];
      int d = ket_b[iket];
      int J2cd = ket_J[iket];
      tbtd(ibra,iket) = TBTD(  J_index_i,  eigvec_i,  J_index_f,  eigvec_f,
                               jorbits[a], jorbits[b],  jorbits[c], jorbits[d],
                                                          J2ab,  J2cd,  Lambda2 );
      if  ((J_index_i == J_index_f) and (eigvec_i==eigvec_f))
      {
        tbtd(iket,ibra) = tbtd(ibra,iket) * (1-abs(J2ab-J2cd)%4);
      }
    }
  }
  return tbtd;
}


arma::mat TransitionDensity::GetOneBodyTransitionOperator( string filename)
{

  ifstream opfile(filename);
  string line,bufstring;
  int Rank_J, Rank_T, parity;

  while ( line.find("Rank_J") == string::npos)   getline(opfile, line);
  istringstream( line.substr( line.rfind(":")+1 ) ) >> Rank_J;
  while ( line.find("Rank_T") == string::npos)   getline(opfile, line);
  istringstream( line.substr( line.rfind(":")+1 ) ) >> Rank_T;
  while ( line.find("Parity") == string::npos)   getline(opfile, line);
  istringstream( line.substr( line.rfind(":")+1 ) ) >> parity;

  vector<int> jorbits;
  for (size_t i=0;i<m_orbits.size();i+= m_orbits[i].j2+1)
  {
    if (m_orbits[i].mj2 == -m_orbits[i].j2) jorbits.push_back(i);
  }

  while ( line.find("index") == string::npos)   getline(opfile, line);
  vector<int> orbits_in;
  getline(opfile, line);
  while ( line.size() > 5 )  // no specific reason why 5. Just looking for the empty comment line
  {
    int index,n,l,j2,tz2;
    istringstream(line.substr(1)) >> index >> n >> l >> j2 >> tz2;
    // find the corresponding j-orbit
    for (size_t i=0;i<jorbits.size(); ++i)
    { // we use -tz2 because we switch isospin conventions
      if (    m_orbits[jorbits[i]].n==n  and m_orbits[jorbits[i]].l2==2*l and m_orbits[jorbits[i]].j2==j2 and m_orbits[jorbits[i]].tz2==-tz2)
      {
         orbits_in.push_back(i);
         cout << orbits_in.size() << " :  " << i << "  " << jorbits[i] << "  " << n << " " << l << " " << j2 << " " << tz2 << endl;
      }
    }

    getline(opfile, line);
  }


  for (size_t i=0;i<orbits_in.size();++i) cout << i << " :  " << orbits_in[i] << endl;

  arma::mat Op1b(jorbits.size(),jorbits.size(),arma::fill::zeros);

  getline(opfile, line); // skip final header
  int a,b;
  double Op_ab;
  while ( opfile.good() )
  {
    opfile >> a >> b >> Op_ab;
    a = orbits_in[a-1]; // fortran indexing...
    b = orbits_in[b-1];
    int j2a = m_orbits[jorbits[a]].j2;
    int j2b = m_orbits[jorbits[b]].j2;
    Op1b(a,b) = Op_ab;
//    Op1b(b,a) = (1 - abs( m_orbits[jorbits[a]].j2 - m_orbits[jorbits[b]].j2 )%4) * Op_ab; // phase factor (-1)^(ja-jb)
    Op1b(b,a) = (1 - abs(j2a-j2b)%4) * Op_ab; // phase factor (-1)^(ja-jb)
  }

  return Op1b;
}



arma::mat TransitionDensity::GetTwoBodyTransitionOperator( string filename)
{

  ifstream opfile(filename);
  string line,bufstring;
  int Rank_J, Rank_T, parity;

  while ( line.find("Rank_J") == string::npos)   getline(opfile, line);
  istringstream( line.substr( line.rfind(":")+1 ) ) >> Rank_J;
  while ( line.find("Rank_T") == string::npos)   getline(opfile, line);
  istringstream( line.substr( line.rfind(":")+1 ) ) >> Rank_T;
  while ( line.find("Parity") == string::npos)   getline(opfile, line);
  istringstream( line.substr( line.rfind(":")+1 ) ) >> parity;

  // set up the j scheme basis from the m scheme basis
  vector<int> jorbits;
  for (size_t i=0;i<m_orbits.size();i+= m_orbits[i].j2+1)
  {
    if (m_orbits[i].mj2 == -m_orbits[i].j2) jorbits.push_back(i);
  }

  // Read in the single particle basis used in the file
  while ( line.find("index") == string::npos)   getline(opfile, line);
  vector<int> orbits_in;
  getline(opfile, line);
  while ( line.size() > 5 and line.find("a")==string::npos )  // no specific reason why 5. Just looking for the empty comment line
  {
    int index,n,l,j2,tz2;
    istringstream(line.substr(1)) >> index >> n >> l >> j2 >> tz2;
    // find the corresponding j-orbit
    for (size_t i=0;i<jorbits.size(); ++i)
    {// we use -tz2 because we switch isospin conventions
      if (    m_orbits[jorbits[i]].n==n  and m_orbits[jorbits[i]].l2==2*l and m_orbits[jorbits[i]].j2==j2 and m_orbits[jorbits[i]].tz2==-tz2)
      {
         orbits_in.push_back(i);
//         cout << orbits_in.size() << " :  " << i << "  " << jorbits[i] << "  " << n << " " << l << " " << j2 << " " << tz2 << endl;
      }
    }
    getline(opfile, line);
  }

  cout << "------------------------------------" << endl;

  for (size_t i=0;i<orbits_in.size();++i) cout << i << " :  " << orbits_in[i] << endl;


  vector<int> ket_a, ket_b, ket_J;

  for (size_t a=0;a<jorbits.size();++a)
  {
    int ja = m_orbits[jorbits[a]].j2;
    for (size_t b=a; b<jorbits.size();++b)
    {      
      int jb = m_orbits[jorbits[b]].j2;
      int Jmin = abs(ja-jb);
      int Jmax = ja+jb;
      for (int J2=Jmin;J2<=Jmax;J2+=2)
      {
        if (a==b and (J2%4)>0) continue;
//        cout << setw(3) << ket_a.size() << ": " << a << " " << b << " " << J2 << endl;
        ket_a.push_back(a);
        ket_b.push_back(b);
        ket_J.push_back(J2);
      }
    }
  }

  arma::mat Op2b(ket_J.size(), ket_J.size(), arma::fill::zeros);

  if ( line.find("Jab")==string::npos )
    getline(opfile, line); // skip final header
  int a,b,c,d,Jab,Jcd;
  double Op_abcd;
  cout << "line before reading matrix elements: " << line << endl;
  while ( opfile >> a >> b >> c >> d >> Jab >> Jcd >> Op_abcd )
  {
    a = orbits_in[a-1]; // fortran indexing...
    b = orbits_in[b-1];
    c = orbits_in[c-1]; // fortran indexing...
    d = orbits_in[d-1];
    Jab *=2;
    Jcd *=2;
    if (a>b)
    {
      swap(a,b);
      Op_abcd *= -(1-abs(m_orbits[jorbits[a]].j2 + m_orbits[jorbits[b]].j2 - Jab)%4);
    }
    if (c>d)
    {
      swap(c,d);
      Op_abcd *= -(1-abs(m_orbits[jorbits[c]].j2 + m_orbits[jorbits[d]].j2 - Jcd)%4);
    }
//    if (a==b) Op_abcd *= 0.5; // since we only store a<=b, c<=d we get a factor 4/(1+delta_ab)(1_delta_cd)
//    if (c==d) Op_abcd *= 0.5; // which combines with the 1/4 out in front

    size_t ibra=0,iket=0;
    while( ibra<ket_a.size() and not( (ket_a[ibra]==a) and (ket_b[ibra]==b) and ket_J[ibra]==Jab) ) ibra++;
    while( iket<ket_a.size() and not( (ket_a[iket]==c) and (ket_b[iket]==d) and ket_J[iket]==Jcd) ) iket++;

    if (ibra==0 and iket==1)
    {
     cout << "aaaaa   " << a << " " << b << " " << c << " " << d << "  " << Jab << " " << Jcd << " "  << Op_abcd << endl;
    }

    if (a==3 and b==3 and c==3 and d==3 and Jab==0 and Jcd==4)
    {
      cout << "bbbb   " << a << " " << b << " " << c << " " << d << " " << Jab << " " << Jcd << " " << Op_abcd << "   " << ibra << " " << iket << endl;
    }

    Op2b(ibra,iket) = Op_abcd;
    if (ibra!=iket)
      Op2b(iket,ibra) = (1 - abs(Jab+Jcd )%4) * Op_abcd; // phase factor (-1)^(Jab-Jcd)
  }

  return Op2b;
}





void TransitionDensity::WriteEGV(string fname)
{

  vector<MschemeOrbit> mscheme_orbits;

  // find all the core orbits
  for (int N=0;N<Nshell;++N) 
  {
    int sumN = (N+1)*(N+2)*(N+3)/3;
    for (int n=0;2*n<=N;++n)
    {
      int l2 = 2*(N-2*n);
      for (int j2=l2+1;j2>=max(0,l2-1);j2-=2)
      {
       for (int mj2=-j2;mj2<=j2;mj2+=2)
       {
        if ( sumN<=Zcore )        mscheme_orbits.emplace_back( MschemeOrbit(n,l2,j2,mj2,+1 ) );
        if ( sumN<=Acore-Zcore )  mscheme_orbits.emplace_back( MschemeOrbit(n,l2,j2,mj2,-1 ) );
       }
      }
    }
  }


  // loop over the valence orbits
  int n_core_orbits = mscheme_orbits.size();
  for (auto& morbit : jbasis_list[0].nubasis_a.m_orbits)
  {
     mscheme_orbits.push_back ( morbit );
  }
  // calculate parity of 0hw configurations
  int N=0;
  for (N=0;(N+1)*(N+2)*(N+3)/3 < Z;++N){};
  int parity_protons = (N%2);
  for (N=0;(N+1)*(N+2)*(N+3)/3 < A-Z;++N){};
  int parity_neutrons = (N%2);
  int parity = 1-2*(( parity_protons*(Z-Zcore) + parity_neutrons*(A-Acore-(Z-Zcore)))%2);
  
  ofstream output(fname);
  
  string interaction_id = "imsrg";
  int Nmax = 0;
  float hw = 20; // this should be irrelevant

  output << left << setw(10) << Z                              << " !  number of protons" << endl;
  output << left << setw(10) << A-Z                            << " !  number of neutrons" << endl;
  output << left << setw(max(10,(int)interaction_id.size()+3)) << interaction_id << " ! interaction id" << endl;
  output << left << setw(10) << hw                             << " ! hbar omega" << endl;
  output << left << setw(10) << Nshell                         << " ! N shells " << endl;
  output << left << setw(10) << mscheme_orbits.size()          << " ! m-scheme orbits " << endl;
  output << left << setw(10) << Nmax                           << " ! Nmax " << endl;
  output << left << setw(10) << amplitudes.size()              << " ! number of basis Slater determinants" << endl;
  output << left << setw(10) << showpos << parity << noshowpos << " ! parity " << endl;
  output << left << setw(10) << MJtot                          << " ! total MJ*2" << endl;
  output << left << setw(10) << total_number_levels            << " ! Number of eigenstates" << endl;
  
  // write out energies, J, T of eigenstates
  for (auto& nuvec : nuvec_list)
  {
   int imax = nuvec.no_level;
   if ( max_states_per_J.find(nuvec.J2) != max_states_per_J.end() ) imax = min(imax,max_states_per_J[nuvec.J2]);
   for (int ilevel=0;ilevel<imax;++ilevel)
   {
     output << right << fixed << setw(12) << setprecision(4) << nuvec.alpha[ilevel] << " " << nuvec.J2/2.0 << " " << 0 << " " << 0.0 << endl;
   }
  }
  
  
  output << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  output << "!!!  now list mscheme single-particle basis !!!" << endl;
  output << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  // write single-particle basis in mscheme
  
  
  
  int s = 1;
  for (auto& morbit : mscheme_orbits)
  {
    output << setw(4) << s << " " << setw(3) << morbit.n << " " << setw(3) << morbit.l2 << " " << setw(3) << morbit.j2 << " " << setw(3) << morbit.mj2 << " " << setw(3) << morbit.tz2 << endl;;
    s++;
  }
  
  output << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  output << "!!!  now list mscheme basis states and amplitudes !!!" << endl;
  output << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  
  
//  cout << "n_core_orbits = " << n_core_orbits << "  mscheme_orbits.size = " << mscheme_orbits.size() << endl;
  
  // Write m-scheme basis occupations and eigenvector coefficients
  size_t bits_per_word = 8*sizeof(mvec_type);
  for (auto it_amp : amplitudes)
  {
    auto& mvec = it_amp.first;
    for (int i=0;i<n_core_orbits;++i) output << i+1 << " ";
    output << " " ;
    for (size_t i=0;i<mscheme_orbits.size()-n_core_orbits;++i)
    {
      if ( (mvec[i/bits_per_word] >> (i%bits_per_word))&0x1)
          output << n_core_orbits + i+1 << " ";
    }
    output << "   ";
    for (auto& ampvec : it_amp.second)
    {
      for (auto amp : ampvec )
      {
        output << scientific << setw(13) << amp << " ";
      }
    }
    output << endl;
  
  }

}




void TransitionDensity::WriteTRDENS_input(string fname)
{
  ofstream outfile(fname);

  outfile << "T        ! specify which states to take?" << endl;
  outfile << "1   " << total_number_levels << "   ! ki,nki" << endl;
  for (int i=1;i<=total_number_levels;++i) outfile << i << " ";
  outfile << endl;
  outfile << "1   " << total_number_levels << "   ! kf,nkf" << endl;
  for (int i=1;i<=total_number_levels;++i) outfile << i << " ";
  outfile << endl;
  outfile << *max_element(begin(Jlist),end(Jlist)) << "      ! jtotal2max" << endl;
  outfile << "0                   ! irestart" << endl;
  outfile << "2                   ! majortot" << endl;
  outfile << "F                   ! irem_cal" << endl;
  outfile << "F                   ! formfcal" << endl;
  outfile << "F                   ! radial  " << endl;
  outfile << "F                   ! momdist " << endl;
  outfile << "T                   ! twobdcal" << endl;
  outfile << "1                   ! ipn     " << endl;
  outfile << "IMSRG.int_iso                 " << endl;
  outfile << "IMSRG_E2_1b.op_iso            " << endl;
  outfile << "IMSRG_E2_2b.op_iso            " << endl;
  outfile << "F                   ! cluster " << endl;
  outfile << "F                   ! antoine " << endl;
  outfile << "24.d0               ! hbar*Omega for antoine file" << endl;
  outfile << "1                   ! #init states in antoine file" << endl;
  outfile << "0  0 " << endl;
  outfile << "1                   ! #fin states in antoine file" << endl;
  outfile << "0  0 " << endl;
  outfile << "T                   ! mbpt_ncsm " << endl;
  outfile << "F                   ! redstick" << endl;
  outfile << "F                   ! mfd_james" << endl;
  outfile << "F                   ! threebdcal" << endl;
  outfile << "F                   ! NCSMC_kernels" << endl;
  outfile << "3                   ! num_of_interaction_files" << endl;
  outfile << "../vrelnp_rgm.int_H2srg-n3lo2.2_2014 " << endl;
  outfile << "../vrelpp_rgm.int_H2srg-n3lo2.2_2014 " << endl;
  outfile << "../vrelnn_rgm.int_H2srg-n3lo2.2_2014 " << endl;
  outfile << "F                   ! V3Nint " << endl;
  outfile << "14 14 14            ! N1_max,N12_max,N123_max " << endl;
  outfile << "chi2b3b400cD-02cE0098_srg0800ho40C_eMax14_EMax14_hwHO016.me3j_bin" << endl;


  outfile.close();


}




void TransitionDensity::GetAZFromFileName(  )
{
  string trimmed_basename = basename.substr( basename.find_last_of("/")+1 );
  string element = trimmed_basename.substr( 0,2);
  cout << trimmed_basename << " " << element << endl;
  auto el_position = find( periodic_table.begin(),periodic_table.end(), element);
  if (el_position == periodic_table.end())
  {
   cout << "ERROR! : could not find " << element << " in periodic table" << endl;
   return;
  }
  Z = el_position - periodic_table.begin();
  istringstream( trimmed_basename.substr(2,2) ) >> A;

  // now guess at what the core should be, assuming a full major oscillator shell valence space
  for (int N=0;(N+1)*(N+2)*(N+3)/3<=Z;++N) Zcore = (N+1)*(N+2)*(N+3)/3; 
  for (int N=0;(N+1)*(N+2)*(N+3)/3<=(A-Z);++N) Acore = Zcore + (N+1)*(N+2)*(N+3)/3; 
  
}



void TransitionDensity::ReadSPfile()
{
  string sp_file_name = sps_file_name.substr(0,sps_file_name.find_last_of(".")) + ".sp";
  ifstream infile(sp_file_name);
  if (not infile.good()) return;
  infile.ignore(256,'\n');
  infile.ignore(256,'\n');
  infile >> Acore >> Zcore;

}





vector<float> operator*(const float lhs, const vector<float>& rhs)
{
  vector<float> vout = rhs;
  for (size_t i=0;i<vout.size();++i) vout[i] *= lhs;
  return vout;
}

vector<float> operator*(const vector<float>& lhs, const vector<float>& rhs)
{
  vector<float> vout = rhs;
  for (size_t i=0;i<vout.size();++i) vout[i] *= lhs[i];
  return vout;
}

vector<float>& operator+=(vector<float>& lhs, const vector<float>& rhs)
{
  for (size_t i=0;i<lhs.size();++i) lhs[i] += rhs[i];
  return lhs;
}







