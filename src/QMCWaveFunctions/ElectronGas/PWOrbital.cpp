#include "PWOrbital.h"

namespace qmcplusplus
{

PWOrbital::PWOrbital(const std::vector<PosType>& kpts_cart)
 : K(kpts_cart)
{
  OrbitalSetSize = kpts_cart.size();
  mK2.resize(OrbitalSetSize);
  for (int ik=0; ik<OrbitalSetSize; ik++)
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
  for (int ik = 0; ik < OrbitalSetSize; ik++)
  {
    sincos(dot(K[ik], r), &sinkr, &coskr);
    pvec[ik]   = ValueType(coskr, sinkr);
    dpvec[ik]  = ValueType(-sinkr, coskr) * K[ik];
    d2pvec[ik] = ValueType(mK2[ik] * coskr, mK2[ik] * sinkr);
  }
}

void PWOrbital::evaluateValue(
  const ParticleSet& P,
  int iat,
  ValueVector_t& pvec)
{
  const PosType& r = P.activeR(iat);
  RealType sinkr, coskr;
  for (int ik=0; ik<OrbitalSetSize; ik++)
  {
    sincos(dot(K[ik], r), &sinkr, &coskr);
    pvec[ik] = ValueType(coskr, sinkr);
  }
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
  for (int iat=first, ik=0;iat<last;iat++,ik++)
  {
    const PosType& r = P.activeR(iat);
    sincos(dot(K[ik], r), &sinkr, &coskr);
    ValueVector_t p(phi[ik], OrbitalSetSize);
    GradVector_t dp(dphi[ik], OrbitalSetSize);
    HessVector_t hess(d2phi_mat[ik], OrbitalSetSize);
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
