//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//
// File created by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//////////////////////////////////////////////////////////////////////////////////////

#include "MoirePotential.h"
#include "OhmmsData/AttributeSet.h"
#include <cmath>

namespace qmcplusplus
{

MoirePotential::Return_t MoirePotential::evaluate(ParticleSet& P)
{
  value_ = 0.0;
  double arg;
  const size_t Nelec = P.getTotalNum();
  for (size_t iel = 0; iel < Nelec; iel++)
  {
    const auto& r = (P.getActivePtcl() == iel) ? P.activeR(iel) : P.R[iel];
    double esum = 0.0;
    for (size_t m=0;m<gvecs.size();m++)
    {
      arg = dot(gvecs[m], r);
      esum += std::cos(arg+phi);
    }
    value_ += esum;
  }
  value_ *= 2*vmoire;
  return value_;
}

bool MoirePotential::put(xmlNodePtr cur)
{
  const RealType BOHR_RADIUS_ANGS = 0.529177210903;
  const RealType AUTOEV = 27.211386245988034;
  const RealType PI = 4. * std::atan(1);
  RealType amoire_in_ang, vmoire_in_mev, phi_in_deg, epsmoire, mstar;
  std::string moire_shape;
  // set defaults
  epsmoire = 1.0;
  mstar = 1.0;
  amoire_in_ang = -1.0;
  vmoire_in_mev = 0.0;
  phi_in_deg = 0.0;
  // read inputs
  OhmmsAttributeSet attrib;
  attrib.add(amoire_in_ang, "amoire_in_ang");
  attrib.add(vmoire_in_mev, "vmoire_in_mev");
  attrib.add(phi_in_deg, "pmoire_in_deg");
  attrib.add(epsmoire, "epsmoire");
  attrib.add(mstar, "mstar");
  attrib.add(moire_shape, "moire_shape");
  attrib.put(cur);
  phi = phi_in_deg/180.0*PI;
  // check inputs
  if (amoire_in_ang < 0)
    throw std::runtime_error("need amoire_in_ang input");
  // use a_B^* and E_h^* units
  amoire = amoire_in_ang/BOHR_RADIUS_ANGS*mstar/epsmoire;
  vmoire = vmoire_in_mev*1e-3/AUTOEV;
  vmoire *= epsmoire*epsmoire/mstar;
  if (moire_shape=="triangular")
  {
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
      }
  } else if (moire_shape=="square")
    {
       RealType bmag =2*PI/amoire; 
       TinyVector<RealType, 2> r1, r2, g1, g2;
       r1={1.0,0.0};
       r2={0.0,1.0};
       RealType fac = 1.0 / (r1[0]*r2[1]-r1[1]*r2[0]);
       g1={r2[1], -r2[0]};
       g2={-r1[1], r1[0]};
       g1=fac*g1*bmag;
       g2=fac*g2*bmag;
       theta = 90.0/180.0*PI;
       std::vector<TinyVector<RealType, 2>> rotation_matrix;
       rotation_matrix.resize(2);
       rotation_matrix[0]={std::cos(theta), -std::sin(theta)};
       rotation_matrix[1]={std::sin(theta), std::cos(theta)};
       gvecs.resize(2);
       for (unsigned i=0; i<2; i++)
       {
         gvecs[i]={g2[0], g2[1], 0.0};
	 TinyVector<RealType, 2> g2_rotated;
	 g2_rotated[0]=rotation_matrix[0][0]*g2[0]+rotation_matrix[0][1]*g2[1];
	 g2_rotated[1]=rotation_matrix[1][0]*g2[0]+rotation_matrix[1][1]*g2[1];
	 g2=g2_rotated;
       }
    } else
       app_error() << " Error: moire shapei: " << moire_shape << "not implemented!!" << std::endl;
       APP_ABORT("");
  return true;
}

bool MoirePotential::get(std::ostream& os) const
{
  os << "External moire potential" << std::endl;
  os << "  aM = " <<  amoire << " bohr*"<< std::endl;
  os << "  vM = " <<  vmoire << " ha*"<< std::endl;
  os << "  phi = " << phi << std::endl;
  os << "  kshell:" << std::endl;
  for (unsigned m=0; m<gvecs.size();m++)
  {
    app_log() << gvecs[m] << std::endl;
  }
  return true;
}

std::unique_ptr<OperatorBase> MoirePotential::makeClone(ParticleSet& P, TrialWaveFunction& psi)
{
  return std::make_unique<MoirePotential>(*this);
}

} // qmcplusplus
