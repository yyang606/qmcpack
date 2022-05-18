#include "PWOrbital.h"

namespace qmcplusplus
{

PWOrbital::PWOrbital(const std::vector<PosType>& kpts_cart)
 : K(kpts_cart), maxk(kpts_cart.size())
{
#ifdef QMC_COMPLEX
  OrbitalSetSize = maxk;
  mink = 0;
#else
  OrbitalSetSize = 2*maxk-1; // k=0 has no (cos, sin) split
  mink = 1;  // treat k=0 as special case
#endif
  mK2.resize(maxk);
  for (int ik=0; ik<maxk; ik++)
  {
    mK2[ik] = -dot(K[ik], K[ik]);
  }
}

PWOrbital::~PWOrbital(){}

void PWOrbital::evaluateVGL(
  const ParticleSet& P,
  int iat,
  ValueVector_t& pvec,
  GradVector_t& dpvec,
  ValueVector_t& d2pvec)
{
  const PosType& r = P.activeR(iat);
  RealType sinkr, coskr;
  for (int ik=mink; ik<maxk; ik++)
  {
    sincos(dot(K[ik], r), &sinkr, &coskr);
#ifdef QMC_COMPLEX
    pvec[ik]   = ValueType(coskr, sinkr);
    dpvec[ik]  = ValueType(-sinkr, coskr) * K[ik];
    d2pvec[ik] = ValueType(mK2[ik] * coskr, mK2[ik] * sinkr);
#else
    const int j2 = 2*ik;
    const int j1 = j2-1;
    pvec[j1]  = coskr;
    pvec[j2]  = sinkr;
    dpvec[j1]  = -sinkr * K[ik];
    dpvec[j2]  = coskr * K[ik];
    d2pvec[j1] = mK2[ik] * coskr;
    d2pvec[j2] = mK2[ik] * sinkr;
#endif
  }
#ifndef QMC_COMPLEX
  pvec[0]   = 1.0;
  dpvec[0]  = 0.0;
  d2pvec[0] = 0.0;
#endif
}

void PWOrbital::evaluateValue(
  const ParticleSet& P,
  int iat,
  ValueVector_t& pvec)
{
  const PosType& r = P.activeR(iat);
  RealType sinkr, coskr;
  for (int ik=mink; ik<maxk; ik++)
  {
    sincos(dot(K[ik], r), &sinkr, &coskr);
#ifdef QMC_COMPLEX
    pvec[ik] = ValueType(coskr, sinkr);
#else
    const int j2 = 2*ik;
    const int j1 = j2-1;
    pvec[j1] = coskr;
    pvec[j2] = sinkr;
#endif
  }
#ifndef QMC_COMPLEX
  pvec[0] = 1.0;
#endif
}

void PWOrbital::evaluate_notranspose(
  const ParticleSet& P,
  int first,
  int last,
  ValueMatrix_t& phi,
  GradMatrix_t& dphi,
  ValueMatrix_t& d2phi)
{
  for (int iat=first, i=0;iat<last;iat++,i++)
  {
    ValueVector_t p(phi[i], OrbitalSetSize);
    GradVector_t dp(dphi[i], OrbitalSetSize);
    ValueVector_t d2p(d2phi[i], OrbitalSetSize);
    evaluateVGL(P, iat, p, dp, d2p);
  }
}

void PWOrbital::evaluate_notranspose(
  const ParticleSet& P,
  int first,
  int last,
  ValueMatrix_t& phi,
  GradMatrix_t& dphi,
  HessMatrix_t& d2phi_mat)
{
  RealType sinkr, coskr;
  ValueType phi_of_r;
  for (int iat=first,i=0;iat<last;iat++,i++)
  {
    const PosType& r = P.activeR(iat);
    ValueVector_t p(phi[i], OrbitalSetSize);
    GradVector_t dp(dphi[i], OrbitalSetSize);
    HessVector_t hess(d2phi_mat[i], OrbitalSetSize);
    for (int ik=mink;ik<maxk;ik++)
    {
      sincos(dot(K[ik], r), &sinkr, &coskr);
#ifdef QMC_COMPLEX
      // phi(r) = cos(kr)+i*sin(kr)
      phi_of_r = ValueType(coskr, sinkr);
      p[ik] = phi_of_r;
      // i*phi(r) = -sin(kr) + i*cos(kr)
      dp[ik] = ValueType(-sinkr, coskr) * K[ik];
      for (int la = 0; la < OHMMS_DIM; la++)
      {
        (hess[ik])(la, la) = -phi_of_r * (K[ik])[la] * (K[ik])[la];
        for (int lb = la + 1; lb < OHMMS_DIM; lb++)
        {
          (hess[ik])(la, lb) = -phi_of_r * (K[ik])[la] * (K[ik])[lb];
          (hess[ik])(lb, la) = (hess[ik])(la, lb);
        }
      }
#else
      const int j2 = 2*ik;
      const int j1 = j2-1;
      p[j1]  = coskr;
      p[j2]  = sinkr;
      dp[j1] = -sinkr * K[ik];
      dp[j2] = coskr * K[ik];
      for (int la = 0; la < OHMMS_DIM; la++)
      {
        (hess[j1])(la, la) = -coskr * (K[ik])[la] * (K[ik])[la];
        (hess[j2])(la, la) = -sinkr * (K[ik])[la] * (K[ik])[la];
        for (int lb = la + 1; lb < OHMMS_DIM; lb++)
        {
          (hess[j1])(la, lb) = -coskr * (K[ik])[la] * (K[ik])[lb];
          (hess[j2])(la, lb) = -sinkr * (K[ik])[la] * (K[ik])[lb];
          (hess[j1])(lb, la) = (hess[j1])(la, lb);
          (hess[j2])(lb, la) = (hess[j2])(la, lb);
        }
      }
#endif
    }
  }
}

void PWOrbital::report(const std::string& pad) const
{
  app_log() << pad << "PWOrbital report" << std::endl;
  for (int ik=0;ik<K.size();ik++)
  {
    app_log() << pad << ik << " " << K[ik] << std::endl;
  }
  app_log() << pad << "end PWOrbital report" << std::endl;
  app_log().flush();
}
} // qmcplusplus
