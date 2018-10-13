
#include <numeric>
#include <cmath>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <bitset>
#include <gsl/gsl_sf_coupling.h>
#include <omp.h>

#include "JMState.hh"


const int JMState::gwords = 1;


JMState::JMState()
: J2(0),T2(0),M2(0),pindx(0),m_coefs({})
{
}


JMState::JMState(const NuBasis& nubasis, const NuProj& nuproj, int istate)
 : J2(nuproj.j[istate]), T2(nuproj.t[istate]), M2(nuproj.j[istate]),pindx(1,nuproj.pindx[istate]),
  partition(nubasis.partition),m_orbits(nubasis.m_orbits)
{
//  std::cout << "making a JMState.   looping iibf from zero to  " << nubasis.ibf[pindx[0]-1] << std::endl;
//  std::cout << "   int the JMState constructor, J2,M2 = " << J2 << " " << M2 << std::endl;
  for (int iibf=0;iibf<nubasis.ibf[pindx[0]-1];++iibf)
  {
    m_coefs[nubasis.vec[pindx[0]-1][iibf]] += nuproj.coef_st[istate][iibf];
//    std::cout << "just added an element iibf = " << iibf << " to m_coefs.  size of m_coefs = " << m_coefs.size() << std::endl;
  }
#ifdef VERBOSE
  std::cout << "JMState::JMstate  about to call EliminateZeros. istate = " << istate << std::endl;
#endif
  EliminateZeros();
#ifdef VERBOSE
  std::cout << "JMState::JMstate done" << std::endl;
#endif
}




void JMState::Print() const
{
  std::cout << "2J = " << J2 << std::endl;
  std::cout << "2MJ = " << M2 << std::endl;
  std::cout << "Partition index = ";
  for (size_t i=0;i<pindx.size()-1;++i) std::cout << pindx[i] << " x ";
  std::cout << pindx.back() << std::endl;
  float sum_coef = 0;
  for ( auto it_state : m_coefs )
  {
    std::cout << std::fixed << std::setw(11) << std::setprecision(8) << it_state.first << " :  "
         << std::fixed << std::setw(11) << std::setprecision(8) << it_state.second;
    sum_coef += it_state.second*it_state.second;
    for (int g=0;g<gwords;++g)
    {
       std::cout << "   " << PrintMstate(it_state.first) << std::endl;
    }
           
  }
  std::cout << "normalization: " << sqrt(sum_coef) << std::endl;
  std::cout << std::endl;

}

std::string JMState::PrintMstate(key_type vec_in) const
{
  std::ostringstream mstr;
  std::vector<int> spaces;
  int iorb = 0;
  for (auto orb : m_orbits )
  {
    if (orb.mj2 == -orb.j2)  spaces.push_back(iorb);
    iorb++;
  }
  spaces.push_back(iorb);
  
  for (int ispace=0;ispace<spaces.size()-1;ispace++)
  {
    for (int ibit=spaces[ispace];ibit<spaces[ispace+1];++ibit)
    {
      mstr << vec_in[ibit];
    }
    mstr << " ";
  }
  return mstr.str();
}


JMState JMState::Jplus()
{
  JMState jmout = *this;
  jmout.m_coefs.clear();
  jmout.M2 +=2;
  if (jmout.M2 > jmout.J2) return jmout;
  for (auto& it_mstate : m_coefs)
  {
   auto& mvec_in = it_mstate.first;
   float m_in_coef = it_mstate.second;
   std::vector<key_type> mvec_out;
   std::vector<float> coefs;
   for (size_t i_m=0;i_m<m_orbits.size();++i_m)  
   {
      if ( (not mvec_in[i_m]) or (mvec_in[i_m+1])  ) continue; // Pauli principle
      
      int j2  = m_orbits[i_m].j2;
      int mj2 = m_orbits[i_m].mj2;
      if (mj2==j2) continue;
      key_type temp_mvec_out = mvec_in;
      temp_mvec_out.set(i_m,0).set(i_m+1,1);
      jmout.m_coefs[ temp_mvec_out ] += m_in_coef * sqrt( j2*(j2+2)-mj2*(mj2+2) )*0.5;
   }

  }
  jmout.EliminateZeros();
  jmout.Normalize();

  return jmout;
}



JMState JMState::Jminus()
{
  JMState jmout = *this;
  jmout.m_coefs.clear();
  jmout.M2 -=2;
  if (jmout.M2 < -jmout.J2) return jmout;
  for (auto& it_mstate : m_coefs)
  {
   auto& mvec_in = it_mstate.first;
   float m_in_coef = it_mstate.second;
   std::vector<key_type> mvec_out;
   std::vector<float> coefs;
   for (size_t i_m=1;i_m<m_orbits.size();++i_m)  // don't bother starting with 0, since it's already in the lowest m_j state
   {
      if ( (not mvec_in[i_m]) or mvec_in[i_m-1]  ) continue; // Pauli principle
      
      int j2  = m_orbits[i_m].j2;
      int mj2 = m_orbits[i_m].mj2;
      if (mj2==-j2) continue;
      key_type temp_mvec_out = mvec_in;
      temp_mvec_out.set(i_m,0).set(i_m-1,1);
      jmout.m_coefs[ temp_mvec_out ] += m_in_coef * sqrt( j2*(j2+2)-mj2*(mj2-2) )*0.5;
   }

  }

  jmout.EliminateZeros();
  jmout.Normalize();
  
  return jmout;

}


// TODO: This must be missing a phase factor
// or else it's not exactly doing what I want...
JMState JMState::TimeReverse()
{
  JMState jmout(*this);
  jmout.m_coefs.clear();
  jmout.M2 = -M2;

  size_t i_m;
  for (auto& it_mstate : m_coefs)
  {
   auto& mvec_in = it_mstate.first;
   float m_in_coef = it_mstate.second;
   i_m=0;
   key_type mvec_out = mvec_in;
   while (i_m<m_orbits.size()) 
   {
//      int iword = i_m /(sizeof(mvec_type)*8);
//      size_t i_m_local = i_m % (sizeof(mvec_type)*8);
      int mj2 = m_orbits[i_m].mj2;
      int j2 = m_orbits[i_m].j2;
      if (mj2>0)
      {
        i_m += (j2-mj2)/2+1;
        continue;
      }
      size_t i_r = i_m - mj2; // time reversed m-state << TODO: This isn't right
//      int iword_r = i_r /(sizeof(mvec_type)*8);
      i_r = i_r %(sizeof(mvec_type)*8);
//      if ( (mvec_out[iword]>>i_m_local)&0x1L ^ (mvec_out[iword_r]>>i_r)&0x1L)
      if ( mvec_out[i_m] ^ mvec_out[i_r])
      {
//        mvec_out[iword] ^= (0x1L << i_m_local);
//        mvec_out[iword_r] ^= (0x1L << i_r);
        mvec_out.flip(i_m);
        mvec_out.flip(i_r);
      }
      i_m++;
   }
   jmout.m_coefs[mvec_out] = m_in_coef;
  }
  return jmout;
}



void JMState::RotateToM( int M_in )
{
  JMState& jm = *this;
//  if (std::abs(jm.M2-M_in)>std::abs(jm.M2+M_in))
//  {
////    jm = jm.TimeReverse();
//  }
//  std::cout << "In " << __func__ << "  rotating from " << M2 << " to " << M_in << std::endl;
  while( jm.M2 > M_in) jm = jm.Jminus();
//  std::cout << "After Jminus loop, M = " << M2 << " note jm.M2 = " << jm.M2 << " is it less than " << M_in << " ? " << (jm.M2<M_in) << std::endl;
  while( jm.M2 < M_in) jm = jm.Jplus();
//  std::cout << "Rotated to M_in. M = " << jm.M2 << std::endl;
}



void JMState::Normalize()
{
  float sum_coef = 0;
  for ( auto& it_mstate : m_coefs)  sum_coef += it_mstate.second * it_mstate.second;
  float norm = sqrt(sum_coef);
  if (norm > 1e-8)
     for ( auto& it_mstate : m_coefs)  it_mstate.second /= norm;
}


// Eliminate std::vectors with zero coefficients
void JMState::EliminateZeros()
{
  std::vector<key_type> zero_keys;
  for ( auto& it_mstate : m_coefs)
  {
    if (std::abs(it_mstate.second)<1e-6)
      zero_keys.push_back(it_mstate.first);
  }
  for (auto& k : zero_keys) m_coefs.erase(k);

}





JMState& JMState::operator=(const JMState& rhs) = default;
JMState& JMState::operator=( JMState&& rhs) = default;

JMState& JMState::operator+=(const JMState& rhs)
{
  for (auto& it_mvec : rhs.m_coefs)
  {
    m_coefs[ it_mvec.first ] += it_mvec.second;
  }
  return *this;
}

JMState JMState::operator+(const JMState& rhs)
{
   JMState jmout = JMState(*this);
   jmout += rhs;
   return jmout;
}

JMState& JMState::operator*=(const double rhs)
{
  for (auto& it_mvec : m_coefs)
  {
    it_mvec.second *= rhs;
  }
  return *this;
}

JMState JMState::operator*(const double rhs)
{
   JMState jmout = JMState(*this);
   jmout *= rhs;
   return jmout;
}



float JMState::Norm() const
{
  float norm = 0;
  for (auto& it_m : m_coefs)
  {
    norm += it_m.second * it_m.second;
  }
  return sqrt(norm);
}



JMState JMState::OuterProduct( const JMState& rhs ) const
{
  JMState jmout;
  jmout.m_coefs.clear();
  for (auto& it_m1 : m_coefs)
  {
    for (auto& it_m2 : rhs.m_coefs)
    {
      if ( (it_m1.first & it_m2.first) != 0) continue;  // cant take an outer product if the two states share an orbit
      double coef = it_m1.second * it_m2.second;
      if (std::abs(coef)>1e-8)
      {
        key_type key = it_m1.first | it_m2.first;
        jmout.m_coefs[ key ] = coef;
      }
    }
  }
  return jmout;
}




// Non-class method, multiply from the left by a scalar,
// which is the same as multiplying from the right.
JMState operator*(const double lhs, const JMState& rhs)
{
  JMState jmout(rhs);
  jmout *= lhs;
  return jmout;
}



// Clebsch-Gordan coefficient, accepts integers equal to twice the corresponding jvalues
double CG(int j2a, int m2a, int j2b, int m2b, int J2, int M2)
{
  return (1-std::abs(j2a-j2b+M2)%4) * sqrt(J2+1) * gsl_sf_coupling_3j(j2a,j2b,J2,m2a,m2b,-M2);
}



JMState TensorProduct( JMState jm1, JMState jm2, int J, int M)
{
//  std::cout << "Calling TensorProduct to combine two jmstates with j2,m2 = " << jm1.J2 << " " << jm1.M2 << "  " << jm2.J2 << " "  << jm2.M2 << "   to J = " << J << "  M = " << M << std::endl;
  int m1 = std::min(jm1.J2, M+jm2.J2);
  int m2 = M-m1;
  int m1_min = std::max(-jm1.J2, M-jm2.J2);

  JMState jmout;
  jmout.SetJ(J);
  jmout.SetM(M);
  for (auto p : jm1.pindx) jmout.pindx.push_back(p);
  for (auto p : jm2.pindx) jmout.pindx.push_back(p);


//  std::cout << "  before rotating, M1,M2 = " << jm1.M2 << " " << jm2.M2 << std::endl;

  // start jm1 and jm2 in the proper m states
  jm1.RotateToM(m1);
  jm2.RotateToM(m2);
//  std::cout << " I thought I was rotating to m1,m2 = " << m1 << " " << m2 << std::endl;
//  std::cout << "  I have " << jm1.M2 << "  " << jm2.M2 << std::endl;

//  std::cout << "In TensorProduct, begin while  " << jm1.M2 << " >= " << m1_min << std::endl;

  while(jm1.M2 >= m1_min)
  {
    double clebsch = CG(jm1.J2, jm1.M2, jm2.J2, jm2.M2, J, M);
//    std::cout << "     jm1.M2 = " << jm1.M2 << "  jm2.M2 = " << jm2.M2 << "   clebsch = " << clebsch  << std::endl;

    if (std::abs(clebsch)>1e-4)
       jmout += clebsch * ( jm1.OuterProduct(jm2) ) ;

//    std::cout << "      size of mcoefs = " << jmout.m_coefs.size() << std::endl;

    jm1 = jm1.Jminus();
    jm2 = jm2.Jplus();
//    std::cout << "after Jminus,Jplus, jm2.M2 = " << jm2.M2 << "  jm2.J2 = " << jm2.J2 << std::endl;
    if (jm2.M2 > jm2.J2) break;
  }

  return jmout;

}






