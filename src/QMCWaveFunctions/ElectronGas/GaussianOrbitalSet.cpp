#include "GaussianOrbitalSet.h"
#include "Particle/DistanceTableData.h"

namespace qmcplusplus
{

GaussianOrbitalSet::GaussianOrbitalSet(ParticleSet& els, ParticleSet& ions, RealType cexpo_in, int ndim_in)
 : targetPtcl(els),
   sourcePtcl(ions),
   cexpo(cexpo_in),
   ideitab(els.addTable(ions)),
   ndim(ndim_in)
{
   OrbitalSetSize = ions.getTotalNum();
}

GaussianOrbitalSet::~GaussianOrbitalSet(){}

void GaussianOrbitalSet::evaluate_notranspose(
  const ParticleSet& P,
  int first,
  int last,
  ValueMatrix_t& phi,
  GradMatrix_t& dphi,
  ValueMatrix_t& d2phi)
{
  for (int i=first;i<last;i++)
  {
    ValueVector_t p(phi[i], OrbitalSetSize);
    GradVector_t dp(dphi[i], OrbitalSetSize);
    ValueVector_t d2p(d2phi[i], OrbitalSetSize);
    evaluateVGL(P, i, p, dp, d2p);
  }
}

void GaussianOrbitalSet::evaluateVGL(
    const ParticleSet& P,
    int i,
    ValueVector_t& pvec,
    GradVector_t& dpvec,
    ValueVector_t& d2pvec)
{
  const auto& dei = P.getDistTable(ideitab);
  const auto& dist = dei.getDistRow(i);
  const auto& displ = dei.getDisplRow(i);
  RealType rij;
  PosType drij;
  for (int j=0;j<OrbitalSetSize;j++)
  {
    rij = dist[j];
    drij = displ[j];
    pvec[j] = std::exp(-cexpo*rij);
    for (int l=0;l<ndim;l++)
      dpvec[j][l] = -2.0*cexpo*drij[l];
    d2pvec[j] = 4.0*cexpo*cexpo*rij*rij-2*ndim*cexpo;
  }
}

SPOSet* GaussianOrbitalSet::makeClone() const { return new GaussianOrbitalSet(*this); }

void GaussianOrbitalSet::report(const std::string& pad) const
{
  app_log() << pad << "GaussianOrbital report" << std::endl;
  app_log() << pad << "  source   = " << sourcePtcl.getName() << std::endl;
  app_log() << pad << "  cexpo    = " << cexpo << std::endl;
  app_log() << pad << "  ndim     = " << ndim << std::endl;
  app_log() << pad << "end GaussianOrbital report" << std::endl;
  app_log().flush();
}
} // qmcplusplus
