#include "catch.hpp"

#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/CoulombPBCAA.h"

namespace qmcplusplus
{

TEST_CASE("Coulomb PBC A-A Ewald quasi 2D square", "[hamiltonian]")
{
  LRCoulombSingleton::CoulombHandler = 0; // !!!! crucial if not first test
  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::QUASI2D;
  const double vmad_sq = -1.95013246;
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true;
  lattice.BoxBConds[2] = false; // ppn
  lattice.ndim = 2;
  lattice.R.diagonal(1.0);
  lattice.LR_dim_cutoff = 30.0;
  lattice.reset();

  ParticleSet elec;
  elec.Lattice = lattice;
  elec.setName("e");
  elec.create({1});
  elec.R[0] = {0.0, 0.0, 0.0};

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();
  elec.addTable(elec);
  elec.update();

  CoulombPBCAA caa = CoulombPBCAA(elec, true);

  double val = caa.evaluate(elec);
  CHECK(val == Approx(vmad_sq));
}

TEST_CASE("Coulomb PBC A-A Ewald quasi 2D body center", "[hamiltonian]")
{
  LRCoulombSingleton::CoulombHandler = 0; // !!!! crucial if not first test
  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::QUASI2D;
  const double vmad_bc = -2.7579038;
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true;
  lattice.BoxBConds[2] = false; // ppn
  lattice.ndim = 2;
  lattice.R.diagonal(1.0);
  lattice.LR_dim_cutoff = 30.0;
  lattice.reset();

  ParticleSet elec;
  elec.Lattice = lattice;
  elec.setName("e");
  const int npart = 2;
  elec.create({npart});
  elec.R[0] = {0.0, 0.0, 0.0};
  elec.R[1] = {0.5, 0.5, 0.0};

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();
  elec.addTable(elec);
  elec.update();

  CoulombPBCAA caa = CoulombPBCAA(elec, true);

  double val = caa.evaluate(elec);
  CHECK(val/npart == Approx(vmad_bc));
}

TEST_CASE("Coulomb PBC A-A Ewald quasi 2D triangle", "[hamiltonian]")
{
  LRCoulombSingleton::CoulombHandler = 0; // !!!! crucial if not first test
  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::QUASI2D;
  const double vmad_tri = -1.106102587;
  const double alat = std::sqrt(2.0*M_PI/std::sqrt(3));
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true;
  lattice.BoxBConds[2] = false; // ppn
  lattice.ndim = 2;
  lattice.R = 0.0;
  lattice.R(0, 0) = alat;
  lattice.R(1, 0) = -1.0/2*alat;
  lattice.R(1, 1) = std::sqrt(3)/2*alat;
  lattice.R(2, 2) = 2*alat;
  lattice.LR_dim_cutoff = 30.0;
  lattice.reset();

  ParticleSet elec;
  elec.Lattice = lattice;
  elec.setName("e");
  const int npart = 1;
  elec.create({npart});
  elec.R[0] = {0.0, 0.0, 0.0};

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();
  elec.addTable(elec);
  elec.update();

  CoulombPBCAA caa = CoulombPBCAA(elec, true);

  double val = caa.evaluate(elec);
  CHECK(val/npart == Approx(vmad_tri));
}

TEST_CASE("Coulomb PBC A-A Ewald quasi 2D honeycomb", "[hamiltonian]")
{
  LRCoulombSingleton::CoulombHandler = 0; // !!!! crucial if not first test
  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::QUASI2D;
  const double vmad_hon = -1.510964233;
  const double alat = std::sqrt(2.0*M_PI/std::sqrt(3));
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true;
  lattice.BoxBConds[2] = false; // ppn
  lattice.ndim = 2;
  lattice.R = 0.0;
  lattice.R(0, 0) = alat;
  lattice.R(1, 0) = -1.0/2*alat;
  lattice.R(1, 1) = std::sqrt(3)/2*alat;
  lattice.R(2, 2) = 2*alat;
  lattice.LR_dim_cutoff = 30.0;
  lattice.reset();

  ParticleSet elec;
  elec.Lattice = lattice;
  elec.setName("e");
  const int npart = 2;
  elec.create({npart});
  // initialize fractional coordinates
  elec.R[0] = {0.0, 0.0, 0.0};
  elec.R[1] = {2./3, 1./3, 0.0};
  // convert to Cartesian coordinates
  for (int i=0;i<npart;i++)
    elec.R[i] = dot(elec.R[i], lattice.R);

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();
  elec.addTable(elec);
  elec.update();

  CoulombPBCAA caa = CoulombPBCAA(elec, true);

  double val = caa.evaluate(elec);
  CHECK(val/npart == Approx(vmad_hon));
}

TEST_CASE("Coulomb PBC A-A Ewald staggered triangle", "[hamiltonian]")
{
  LRCoulombSingleton::CoulombHandler = 0; // !!!! crucial if not first test
  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::QUASI2D;
  const double rs = 20.0;
  const double alat = rs*std::sqrt(2.0*M_PI/std::sqrt(3));
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true;
  lattice.BoxBConds[2] = false; // ppn
  lattice.ndim = 2;
  lattice.R = 0.0;
  lattice.R(0, 0) = alat;
  lattice.R(1, 0) = -1.0/2*alat;
  lattice.R(1, 1) = std::sqrt(3)/2*alat;
  lattice.R(2, 2) = 2*alat;
  lattice.LR_dim_cutoff = 30.0;
  lattice.reset();

  ParticleSet elec;
  elec.Lattice = lattice;
  elec.setName("e");
  const int npart = 2;
  elec.create({npart});
  // initialize fractional coordinates
  elec.R[0] = {0.0, 0.0, 0.0};
  elec.R[1] = {2./3, 1./3, 0.0};
  // convert to Cartesian coordinates
  for (int i=0;i<npart;i++)
    elec.R[i] = dot(elec.R[i], lattice.R);

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();
  elec.addTable(elec);
  elec.update();

  CoulombPBCAA caa = CoulombPBCAA(elec, true);

  const int ntest = 4;
  TinyVector<double, ntest> zheight = {0, 0.1, 0.5, 3.0};
  zheight *= rs;
  const double vmad_hon = -1.510964233;
  const double vmad_tri = -1.106102587;
  TinyVector<double, ntest> vmad_ref = {vmad_hon, -1.4193042644, -1.2005504968, vmad_tri};
  vmad_ref /= rs;
  double val;
  for (int itest=0; itest<ntest; itest++)
  {
    for (int i=npart/2;i<npart;i++)
      elec.R[i][2] = zheight[itest];
    elec.update();
    val = caa.evaluate(elec);
    CHECK(val/npart == Approx(vmad_ref[itest]));
  }
}

} // qmcplusplus
