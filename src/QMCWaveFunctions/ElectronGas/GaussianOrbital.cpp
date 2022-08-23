#include "GaussianOrbital.h"
#include "Particle/DistanceTableData.h"

namespace qmcplusplus
{

GaussianOrbital::GaussianOrbital(ParticleSet& els, ParticleSet& ions, RealType cexpo_in, int ndim_in)
 : targetPtcl(els),
   sourcePtcl(ions),
   cexpo(cexpo_in),
   ideitab(els.addTable(ions)),
   ndim(ndim_in)
{
  OrbitalSetSize = ions.getTotalNum();
}

GaussianOrbital::~GaussianOrbital(){}

GaussianOrbital::RealType GaussianOrbital::operator()(const RealType rij)
{
  return std::exp(-cexpo*rij*rij);
}

void GaussianOrbital::gradient_log(GradType& dp, const RealType& rij, const PosType& drij)
{
  ValueType p = (*this)(rij);
  for (int l=0;l<ndim;l++)
    dp[l] = -2.0*cexpo*drij[l];
}

void GaussianOrbital::hessian_log(HessType& h, const GradType& dp)
{
  h = outerProduct(dp, dp);
  for (int la=0;la<ndim;la++)
    h(la, la) -= 2*cexpo;
}

void GaussianOrbital::gradHess_log(GGGType& g3, const GradType& dp)
{
  g3 = outerProduct(dp, dp);
  for (int l=0;l<ndim;l++)
    g3[l] *= dp[l];
  const ValueType g3pre=-2.0*cexpo;
  if (ndim == 2) {
    // x
    g3[0](0, 0) += g3pre*dp[0]*3.0;
    g3[0](0, 1) += g3pre*dp[1];
    g3[0](1, 0) += g3pre*dp[1];
    g3[0](1, 1) += g3pre*dp[0];
    // y
    g3[1](0, 0) += g3pre*dp[1];
    g3[1](0, 1) += g3pre*dp[0];
    g3[1](1, 0) += g3pre*dp[0];
    g3[1](1, 1) += g3pre*dp[1]*3.0;
  } else if (ndim == 3) {
    // x
    g3[0](0, 0) += g3pre*dp[0]*3.0;
    g3[0](0, 1) += g3pre*dp[1];
    g3[0](0, 2) += g3pre*dp[2];
    g3[0](1, 0) += g3pre*dp[1];
    g3[0](1, 1) += g3pre*dp[0];
    g3[0](1, 2) += 0;
    g3[0](2, 0) += g3pre*dp[2];
    g3[0](2, 1) += 0;
    g3[0](2, 2) += g3pre*dp[0];
    // y
    g3[1](0, 0) += g3pre*dp[1];
    g3[1](0, 1) += g3pre*dp[0];
    g3[1](0, 2) += 0;
    g3[1](1, 0) += g3pre*dp[0];
    g3[1](1, 1) += g3pre*dp[1];
    g3[1](1, 2) += g3pre*dp[2];
    g3[1](2, 0) += 0;
    g3[1](2, 1) += g3pre*dp[2];
    g3[1](2, 2) += g3pre*dp[1];
    // z
    g3[2](0, 0) += 0;
    g3[2](0, 1) += g3pre*dp[2];
    g3[2](0, 2) += g3pre*dp[1];
    g3[2](1, 0) += g3pre*dp[2];
    g3[2](1, 1) += 0;
    g3[2](1, 2) += g3pre*dp[0];
    g3[2](1, 0) += g3pre*dp[0];
    g3[2](1, 1) += g3pre*dp[1];
    g3[2](1, 2) += g3pre*dp[2]*3.0;
  } else {
    throw std::runtime_error("GaussianOrbital ndim not implemented");
  }
}

void GaussianOrbital::evaluateValue(
  const ParticleSet& P,
  int i,
  ValueVector_t& pvec)
{
  const auto& d_table = P.getDistTable(ideitab);
  const auto& dist  = (P.activePtcl == i) ? d_table.getTempDists() : d_table.getDistRow(i);
  RealType rij;
  for (int j=0;j<OrbitalSetSize;j++)
  {
    rij = dist[j];
    pvec[j] = (*this)(rij);
  }
}

// three variants of evaluate_notranspose follow:
//   1. value, gradient, laplacian (diagonal of hessian)
//   2. value, gradient, hessian
//   3. value, gradient, hessian, grad(hessian)

void GaussianOrbital::evaluate_notranspose(
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

void GaussianOrbital::evaluateVGL(
  const ParticleSet& P,
  int i,
  ValueVector_t& pvec,
  GradVector_t& dpvec,
  ValueVector_t& d2pvec)
{
  const auto& d_table = P.getDistTable(ideitab);
  const auto& dist  = (P.activePtcl == i) ? d_table.getTempDists() : d_table.getDistRow(i);
  const auto& displ = (P.activePtcl == i) ? d_table.getTempDispls() : d_table.getDisplRow(i);
  RealType rij;
  PosType drij;
  for (int j=0;j<OrbitalSetSize;j++)
  {
    rij = dist[j];
    drij = -displ[j];
    pvec[j] = (*this)(rij);
    gradient_log(dpvec[j], rij, drij);
    d2pvec[j] = (dot(dpvec[j], dpvec[j])-2*ndim*cexpo)*pvec[j];
    dpvec[j] *= pvec[j];
  }
}

void GaussianOrbital::evaluate_notranspose(const ParticleSet& P,
                                           int first,
                                           int last,
                                           ValueMatrix_t& phi,
                                           GradMatrix_t& dphi,
                                           HessMatrix_t& d2phi_mat)
{
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    ValueVector_t p(phi[i], OrbitalSetSize);
    GradVector_t dp(dphi[i], OrbitalSetSize);
    HessVector_t hess(d2phi_mat[i], OrbitalSetSize);

    const auto& d_table = P.getDistTable(ideitab);
    const auto& dist  = (P.activePtcl == i) ? d_table.getTempDists() : d_table.getDistRow(i);
    const auto& displ = (P.activePtcl == i) ? d_table.getTempDispls() : d_table.getDisplRow(i);

    RealType rij;
    PosType drij;
    for (int j=0;j<OrbitalSetSize;j++)
    {
      rij = dist[j];
      drij = -displ[j];
      p[j] = (*this)(rij);
      gradient_log(dp[j], rij, drij);
      // second derivative
      hessian_log(hess[j], dp[j]);
      hess[j] *= p[j];
      // finish first derivative
      dp[j] *= p[j];
    }
  }
}

void GaussianOrbital::evaluate_notranspose(const ParticleSet& P,
                                           int first,
                                           int last,
                                           ValueMatrix_t& phi,
                                           GradMatrix_t& dphi,
                                           HessMatrix_t& d2phi_mat,
                                           GGGMatrix_t& d3phi_mat)
{
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    ValueVector_t p(phi[i], OrbitalSetSize);
    GradVector_t dp(dphi[i], OrbitalSetSize);
    HessVector_t hess(d2phi_mat[i], OrbitalSetSize);
    GGGVector_t ggg(d3phi_mat[i], OrbitalSetSize);
    const auto& d_table = P.getDistTable(ideitab);
    const auto& dist  = (P.activePtcl == i) ? d_table.getTempDists() : d_table.getDistRow(i);
    const auto& displ = (P.activePtcl == i) ? d_table.getTempDispls() : d_table.getDisplRow(i);

    RealType rij;
    PosType drij;
    for (int j=0;j<OrbitalSetSize;j++)
    {
      rij = dist[j];
      drij = -displ[j];
      p[j] = (*this)(rij);
      gradient_log(dp[j], rij, drij);
      // second derivative
      hessian_log(hess[j], dp[j]);
      hess[j] *= p[j];
      // third derivative
      gradHess_log(ggg[j], dp[j]);
      ggg[j] *= p[j];
      // finish first derivative
      dp[j] *= p[j];
    }
  }
}

SPOSet* GaussianOrbital::makeClone() const { return new GaussianOrbital(*this); }

void GaussianOrbital::report(const std::string& pad) const
{
  app_log() << pad << "GaussianOrbital report" << std::endl;
  app_log() << pad << "  source   = " << sourcePtcl.getName() << std::endl;
  app_log() << pad << "  cexpo    = " << cexpo << std::endl;
  app_log() << pad << "  ndim     = " << ndim << std::endl;
  app_log() << pad << "end GaussianOrbital report" << std::endl;
  app_log().flush();
}
} // qmcplusplus
