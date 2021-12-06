//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "ElectronGasComplexOrbitalBuilder.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
/** constructor for EGOSet
 * @param norb number of orbitals for the EGOSet
 * @param k list of unique k points in Cartesian coordinate excluding gamma
 * @param k2 k2[i]=dot(k[i],k[i])
 */
EGOSet::EGOSet(const std::vector<PosType>& k, const std::vector<RealType>& k2) : K(k), mK2(k2)
{
  KptMax         = k.size();
  OrbitalSetSize = k.size();
  className      = "EGOSet";
  //assign_energies();
}

EGOSet::EGOSet(const std::vector<PosType>& k, const std::vector<RealType>& k2, const std::vector<int>& d)
    : K(k), mK2(k2)
{
  KptMax         = k.size();
  OrbitalSetSize = k.size();
  className      = "EGOSet";
  //assign_energies();
  //assign_degeneracies(d);
}

ElectronGasComplexOrbitalBuilder::ElectronGasComplexOrbitalBuilder(Communicate* comm, ParticleSet& els)
    : WaveFunctionComponentBuilder(comm, els)
{}


WaveFunctionComponent* ElectronGasComplexOrbitalBuilder::buildComponent(xmlNodePtr cur)
{
  int nc = 0;
  int nup = -1;
  int ndn = -1;
  int ndim = 3;
  PosType twist(0.0);
  OhmmsAttributeSet aAttrib;
  aAttrib.add(nc, "shell");
  aAttrib.add(nup, "nup");
  aAttrib.add(ndn, "ndn");
  aAttrib.add(twist, "twist");
  aAttrib.add(ndim, "ndim");
  aAttrib.put(cur);
  //typedef DiracDeterminant<EGOSet>  Det_t;
  //typedef SlaterDeterminant<EGOSet> SlaterDeterminant_t;
  typedef DiracDeterminant<> Det_t;
  typedef SlaterDet SlaterDeterminant_t;
  int nat = targetPtcl.getTotalNum();
  if (nup < 0) APP_ABORT("nup must be given");
  if (ndn < 0) ndn = nat-nup;
  HEGGrid<RealType, OHMMS_DIM> egGrid(targetPtcl.Lattice);
  if (ndim != 3) app_log() << " setting HEGGrid to " << ndim << " dimensions." << std::endl;
  egGrid.ndim = ndim;
  if (nc == 0)
    nc = egGrid.getShellIndex(nup);
  app_log() << "In ElectronGasComplexOrbitalBuilder:" << std::endl;
  app_log() << "  nup = " << nup << std::endl;
  app_log() << "  ndn = " << ndn << std::endl;
  app_log() << "  ntot = " << nat << std::endl;
  app_log() << "  nc = " << nc << std::endl;
  egGrid.createGrid(nc, nup, twist);
  targetPtcl.setTwist(twist);
  int noffset = 0;
  //create up determinant
  Det_t* updet = new Det_t(std::make_unique<EGOSet>(egGrid.kpt, egGrid.mk2));
  updet->set(noffset, nup);
  noffset += nup;
  SlaterDet* sdet = new SlaterDet(targetPtcl);
  sdet->add(updet, 0);
  if (ndn > 0)
  {
    //create down determinant
    Det_t* downdet = new Det_t(std::make_unique<EGOSet>(egGrid.kpt, egGrid.mk2));
    downdet->set(noffset, ndn);
    noffset += ndn;
    sdet->add(downdet, 1);
  }
  if (noffset < nat)
  {
  app_log() << " !!!! adding more than 2 determinants!" << std::endl;
  //create more up determinant
  Det_t* up1det = new Det_t(std::make_unique<EGOSet>(egGrid.kpt, egGrid.mk2));
  up1det->set(noffset, nup);
  noffset += nup;
  //create more down determinant
  Det_t* dn1det = new Det_t(std::make_unique<EGOSet>(egGrid.kpt, egGrid.mk2));
  dn1det->set(noffset, ndn);
  noffset += ndn;
  //create a Slater determinant
  //SlaterDeterminant_t *sdet  = new SlaterDeterminant_t;
  sdet->add(up1det, 2);
  sdet->add(dn1det, 3);
  }
  if (noffset != nat)
  {
    app_log() << "  noffset = " << noffset << std::endl;
    APP_ABORT("particle number mismatch in ElectronGasComplexOrbitalBuilder");
  }
  return sdet;
}

ElectronGasSPOBuilder::ElectronGasSPOBuilder(ParticleSet& p, Communicate* comm, xmlNodePtr cur)
    : SPOSetBuilder("ElectronGas", comm), has_twist(false), unique_twist(-1.0), egGrid(p.Lattice), spo_node(NULL)
{
  ClassName = "ElectronGasSPOBuilder";
}

SPOSet* ElectronGasSPOBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  app_log() << "ElectronGasSPOBuilder::createSPOSet " << std::endl;
  int nc = 0;
  int ns = 0;
  PosType twist(0.0);
  std::string spo_name("heg");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(ns, "size");
  aAttrib.add(twist, "twist");
  aAttrib.add(spo_name, "name");
  aAttrib.add(spo_name, "id");
  aAttrib.put(cur);
  if (has_twist)
    twist = unique_twist;
  else
  {
    unique_twist = twist;
    has_twist    = true;
    for (int d = 0; d < OHMMS_DIM; ++d)
      has_twist &= (unique_twist[d] + 1.0) > 1e-6;
  }
  if (ns > 0)
    nc = egGrid.getShellIndex(ns);
  if (nc < 0)
  {
    app_error() << "  HEG Invalid Shell." << std::endl;
    APP_ABORT("ElectronGasSPOBuilder::put");
  }
  egGrid.createGrid(nc, ns, twist);
  EGOSet* spo = new EGOSet(egGrid.kpt, egGrid.mk2, egGrid.deg);
  return spo;
}


SPOSet* ElectronGasSPOBuilder::createSPOSetFromIndices(indices_t& indices)
{
  egGrid.createGrid(indices);
  EGOSet* spo = new EGOSet(egGrid.kpt, egGrid.mk2, egGrid.deg);
  return spo;
}


} // namespace qmcplusplus
