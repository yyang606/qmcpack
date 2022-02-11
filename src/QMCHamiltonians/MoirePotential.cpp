#include "MoirePotential.h"
#include "OhmmsData/AttributeSet.h"
#include "Particle/DistanceTableData.h"

namespace qmcplusplus
{
MoirePotential::MoirePotential(ParticleSet& elec, ParticleSet& ions)
  : targetPtcl(elec),
    sourcePtcl(ions),
    ideitab(elec.addTable(ions))
{
}

MoirePotential::Return_t MoirePotential::evaluate(ParticleSet& P)
{
  Value = 0.0;
  double arg;
  const size_t Nelec = P.getTotalNum();
  /*
  Return_t arg;
  const DistanceTableData& d_table(P.getDistTable(ideitab));
  const size_t NumIons = sourcePtcl.getTotalNum();
  for (size_t iel = 0; iel < Nelec; iel++)
  {
    const auto& displ = (P.activePtcl == iel) ? d_table.getTempDispls() : d_table.getDisplRow(iel);
    Return_t esum(0);
    for (size_t iat = 0; iat < NumIons; iat++)
    {
      const auto& drij = displ[iat];
      for (size_t m=0;m<3;m++)
      {
        arg = dot(gvecs[m], drij) + phi;
        esum += 2*std::cos(arg);
      }
    }
    Value += esum;
  }
  */
  app_log() << "vmoire = " << vmoire << std::endl;
  for (size_t iel = 0; iel < Nelec; iel++)
  {
    const auto& r = (P.activePtcl == iel) ? P.activeR(iel) : P.R[iel];
    //Return_t esum(0);
    double esum = 0.0;
    for (size_t m=0;m<gvecs.size();m++)
    {
      arg = dot(gvecs[m], r);
      esum += std::cos(arg+phi);
    }
    Value += esum;
  }
  Value *= 2*vmoire;

  //return vmoire*2*Value; // !!!! return value is NOT used
}

bool MoirePotential::put(xmlNodePtr cur)
{
  const RealType BOHR_RADIUS_ANGS = 0.529177210903;
  const RealType AUTOEV = 27.211386245988034;
  const RealType PI = 4. * std::atan(1);
  RealType amoire_in_ang, vmoire_in_mev, phi_in_deg;
  OhmmsAttributeSet attrib;
  attrib.add(amoire_in_ang, "amoire_in_ang");
  attrib.add(vmoire_in_mev, "vmoire_in_mev");
  attrib.add(phi_in_deg, "phi_in_deg");
  attrib.put(cur);
  phi = phi_in_deg/180.0*PI;
  // !!!! HACK hard-code mstar and eps for now
  RealType mstar=0.35, epsmoire=25.0;
  // end HACK !!!!
  // use a_B^* and E_h^* units
  amoire = amoire_in_ang/BOHR_RADIUS_ANGS*mstar/epsmoire;
  vmoire = vmoire_in_mev*1e-3/AUTOEV;
  vmoire *= epsmoire*epsmoire/mstar;
  app_log() << "aM = " << amoire_in_ang << " " << amoire << std::endl;
  app_log() << "vM = " << vmoire_in_mev << " " << vmoire << std::endl;
  app_log() << "phi = " << phi << std::endl;
  // setup the one shell of gvectors
  RealType bmag = 4*PI/std::sqrt(3)/amoire;
  std::vector<TinyVector<int, 3>> gfracs;
  gfracs.resize(3);
  gfracs[0] = {-1,  0, 0};
  gfracs[1] = { 1, -1, 0};
  gfracs[2] = { 0,  1, 0};
  TinyVector<RealType, 3> b1, b2;
  b1 = {std::sqrt(3)/2, 1.0/2, 0.0};
  b2 = {0.0, 1.0, 0.0};
  b1 = bmag*b1;
  b2 = bmag*b2;
  gvecs.resize(3);
  for (unsigned m=0; m<3;m++)
  {
    for (unsigned l=0; l<3; l++)
    {
      gvecs[m][l] = gfracs[m][0]*b1[l] + gfracs[m][1]*b2[l];
    }
    app_log() << gvecs[m] << std::endl;
  }
}

OperatorBase* MoirePotential::makeClone(ParticleSet& P, TrialWaveFunction& psi)
{
  return new MoirePotential(*this);
}

} // qmcplusplus
