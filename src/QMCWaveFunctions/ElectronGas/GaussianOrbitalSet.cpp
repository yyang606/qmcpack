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
{}

GaussianOrbitalSet::~GaussianOrbitalSet(){}

void GaussianOrbitalSet::evaluate_notranspose(
  const ParticleSet& P,
  int first,
  int last,
  ValueMatrix_t& phi,
  GradMatrix_t& dphi,
  ValueMatrix_t& d2phi)
{
  for (int i=first;i<=last;i++)
  {
    evaluateVGL(P, i, phi[i], dphi[i], d2phi[i]);
  }
}

void GaussianOrbitalSet::evaluateVGL(
  const ParticleSet& P,
  int i,
  ValueVector_t& phi,
  GradVector_t& dphi,
  ValueVector_t& d2phi)
{
  const DistanceTableData& dei(P.getDistTable(ideitab));
  const auto& dist = dei.getDistRow(i);
  const auto& displ = dei.getDisplRow(i);
  for (int j=0;j<phi.size();j++)
  {
    RealType rij = dist[j];
    PosType drij = displ[j];
    phi[j] = std::exp(-cexpo*rij);
    for (int l=0;l<ndim;l++)
      dphi[j][l] = -2.0*cexpo*drij[l];
    d2phi[j] = 4.0*cexpo*cexpo*rij*rij-2*ndim*cexpo;
  }
}

void GaussianOrbitalSet::report(const std::string& pad) const
{
  app_log() << pad << "GaussianOrbital report" << std::endl;
  app_log() << pad << "  cexpo    = " << cexpo << std::endl;
  app_log() << pad << "  ndim     = " << ndim << std::endl;
  app_log() << pad << "end GaussianOrbital report" << std::endl;
  app_log().flush();
}
} // qmcplusplus
