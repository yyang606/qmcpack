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
    drij = displ[j];
    pvec[j] = (*this)(rij);
    for (int l=0;l<ndim;l++)
      dpvec[j][l] = 2.0*cexpo*drij[l];
    d2pvec[j] = (dot(dpvec[j], dpvec[j])-2*ndim*cexpo)*pvec[j];
    dpvec[j] *= pvec[j];
  }
}

GaussianOrbital::RealType GaussianOrbital::operator()(RealType rij)
{
  return std::exp(-cexpo*rij*rij);
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
      drij = displ[j];
      p[j] = (*this)(rij);
      for (int l=0;l<ndim;l++)
        dp[j][l] = 2.0*cexpo*drij[l];
      // second derivative
      for (int la=0;la<ndim;la++)
      {
        hess[j](la, la) = (dp[j][la]*dp[j][la]-2*cexpo);
        for (int lb=la+1;lb<ndim;lb++)
        {
          hess[j](la, lb) = dp[j][la]*dp[j][lb];
          hess[j](lb, la) = hess[j](la, lb);
        }
      }
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
      drij = displ[j];
      p[j] = (*this)(rij);
      for (int l=0;l<ndim;l++)
        dp[j][l] = 2.0*cexpo*drij[l];
      // second derivative
      for (int la=0;la<ndim;la++)
      {
        hess[j](la, la) = (dp[j][la]*dp[j][la]-2*cexpo);
        for (int lb=la+1;lb<ndim;lb++)
        {
          hess[j](la, lb) = dp[j][la]*dp[j][lb];
          hess[j](lb, la) = hess[j](la, lb);
        }
      }
      hess[j] *= p[j];
      // third derivative
      for (int la=0;la<ndim;la++)
      {
        ggg[j][la] = dp[j][la]*hess[j];
      }
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
