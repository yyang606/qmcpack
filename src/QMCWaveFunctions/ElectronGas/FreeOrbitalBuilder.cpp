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

#include "OhmmsData/AttributeSet.h"
#include "LongRange/StructFact.h"
#include "LongRange/KContainer.h"
#include "QMCWaveFunctions/ElectronGas/FreeOrbital.h"
#include "QMCWaveFunctions/ElectronGas/FreeOrbitalBuilder.h"

namespace qmcplusplus
{
FreeOrbitalBuilder::FreeOrbitalBuilder(ParticleSet& els, Communicate* comm, xmlNodePtr cur)
    : SPOSetBuilder("PW", comm), targetPtcl(els)
{}

double FreeOrbitalBuilder::calc_kf(const PtclOnLatticeTraits::ParticleLayout lattice, const size_t nelec)
{
  const size_t ndim = lattice.ndim;
  double rs, kf;
  if (ndim == 2)
  {
    const double sa = M_PI;
    rs = std::pow(lattice.Area/nelec/sa, 1.0/ndim);
    kf = 2.0/rs;
  }
  else if (lattice.ndim == 3)
  {
    const double sa = 4.0*M_PI/3.0;
    rs = std::pow(lattice.Volume/nelec/sa, 1.0/ndim);
    kf = std::pow(9.0*M_PI/2.0, 1.0/ndim)/rs;
  }
  else
  {
    std::ostringstream msg;
    msg << "calc_kf ndim = " << ndim << "\n";
    throw std::runtime_error(msg.str());
  }
  return kf;
}

std::unique_ptr<SPOSet> FreeOrbitalBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  int norb = -1;
  std::string spo_object_name = "PW";
  PosType twist(0.0);
  OhmmsAttributeSet attrib;
  attrib.add(norb, "size");
  attrib.add(twist, "twist");
  attrib.add(spo_object_name, "name");
  attrib.put(cur);

  if (norb < 0)
    throw std::runtime_error("free orbital SPO set require the \"size\" input");

  auto lattice = targetPtcl.getLattice();

  PosType tvec = lattice.k_cart(twist);
#ifdef QMC_COMPLEX
  const int npw = norb;
  targetPtcl.setTwist(twist);
  app_log() << "twist fraction = " << twist << std::endl;
  app_log() << "twist cartesian = " << tvec << std::endl;
#else
  const int npw = std::ceil((norb + 1.0) / 2);
  if (2 * npw - 1 != norb)
  {
    std::ostringstream msg;
    msg << "norb = " << norb << " npw = " << npw;
    msg << " cannot be ran in real PWs (sin, cos)" << std::endl;
    msg << "either use complex build or change the size of SPO set" << std::endl;
    msg << "ideally, set size to a closed shell of PWs." << std::endl;
    throw std::runtime_error(msg.str());
  }
  for (int ldim = 0; ldim < twist.size(); ldim++)
  {
    if (std::abs(twist[ldim]) > 1e-16)
      throw std::runtime_error("no twist for real orbitals");
  }
#endif

  // extract npw k-points from container
  // kpts_cart is sorted by magnitude
  std::vector<PosType> kpts(npw);
  KContainer klists;
  RealType kcut = 1.1*calc_kf(lattice, targetPtcl.getTotalNum());
  klists.updateKLists(lattice, kcut, lattice.ndim, twist);
  if (klists.numk < npw)
  {
    std::ostringstream msg;
    msg << "nklist = " << klists.numk << " npw = " << npw;
    msg << " updateKLists did not create enough PWs" << std::endl;
    throw std::runtime_error(msg.str());
  }

  // k0 is not in kpts_cart
  kpts[0] = tvec;
#ifdef QMC_COMPLEX
  for (int ik = 1; ik < npw; ik++)
  {
    kpts[ik] = klists.kpts_cart[ik - 1];
  }
#else
  const int nktot = klists.kpts.size();
  std::vector<int> mkidx(npw, 0);
  int ik = 1;
  for (int jk = 0; jk < nktot; jk++)
  {
    // check if -k is already chosen
    const int jmk = klists.minusk[jk];
    if (in_list(jk, mkidx))
      continue;
    // if not, then add this kpoint
    kpts[ik]  = klists.kpts_cart[jk];
    mkidx[ik] = jmk; // keep track of its minus
    ik++;
    if (ik >= npw)
      break;
  }
#endif
  auto sposet = std::make_unique<FreeOrbital>(spo_object_name, kpts);
  sposet->report("  ");
  return sposet;
}

bool FreeOrbitalBuilder::in_list(const int j, const std::vector<int> l)
{
  for (int i = 0; i < l.size(); i++)
  {
    if (j == l[i])
      return true;
  }
  return false;
}

} // namespace qmcplusplus
