#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "ParticleIO/XMLParticleIO.h"
#include "ParticleIO/ParticleLayoutIO.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "QMCWaveFunctions/ElectronGas/GaussianOrbitalSet.h"

using std::cout;
using std::endl;

namespace qmcplusplus
{

double value(const GaussianOrbitalSet::PosType r, GaussianOrbitalSet* sposet, ParticleSet& P)
{
  const int ielec = 0;
  const int norb = sposet->getOrbitalSetSize();
  GaussianOrbitalSet::ValueVector_t p(norb);
  GaussianOrbitalSet::PosType dr = r-P.R[ielec];
  P.makeMove(ielec, dr);
  sposet->evaluateValue(P, ielec, p);
  return p[ielec].real();
}

GaussianOrbitalSet::PosType fd_gradient(const GaussianOrbitalSet::PosType r, const double h, GaussianOrbitalSet* sposet, ParticleSet& P)
{
  const int ndim = 3;
  GaussianOrbitalSet::PosType g(0);
  GaussianOrbitalSet::PosType r1;
  double fm, fp;
  for (int l=0;l<ndim;l++)
  {
    r1 = r;
    r1[l] -= h;
    fm = value(r1, sposet, P);
    r1[l] += 2*h;
    fp = value(r1, sposet, P);
    g[l] = (fp-fm)/(2*h);
  }
  return g;
}

double fd_laplacian(const GaussianOrbitalSet::PosType r, const double h, GaussianOrbitalSet* sposet, ParticleSet& P)
{
  const int ndim = 3;
  const double f0 = value(r, sposet, P);
  double lap = 0.0;
  double fm, fp;
  GaussianOrbitalSet::PosType r1;
  for (int l=0;l<ndim;l++)
  {
    r1 = r;
    r1[l] -= h;
    fm = value(r1, sposet, P);
    r1[l] += 2*h;
    fp = value(r1, sposet, P);
    lap += (fp+fm-2*f0);
  }
  return lap/(h*h);
}

TEST_CASE("Gaussian orbital for 3D HEG", "[wavefunction]")
{
  const char* particles = "<tmp> \
  <simulationcell>\
     <parameter name=\"lattice\" units=\"bohr\">\
              3.77945227        0.00000000        0.00000000\
              0.00000000        3.77945227        0.00000000\
              0.00000000        0.00000000        3.77945227\
     </parameter>\
     <parameter name=\"bconds\">\
        p p p\
     </parameter>\
     <parameter name=\"LR_dim_cutoff\"       >    15                 </parameter>\
  </simulationcell>\
  <particleset name=\"e\">\
     <group name=\"u\" size=\"1\" mass=\"1.0\">\
        <parameter name=\"charge\"              >    -1                    </parameter>\
        <parameter name=\"mass\"                >    1.0                   </parameter>\
        <attrib name=\"position\" datatype=\"posArray\" condition=\"0\">\
                 0.00000000        0.00000000        0.00000000\
        </attrib>\
     </group>\
     <group name=\"d\" size=\"1\" mass=\"1.0\">\
        <parameter name=\"charge\"              >    -1                    </parameter>\
        <parameter name=\"mass\"                >    1.0                   </parameter>\
        <attrib name=\"position\" datatype=\"posArray\" condition=\"0\">\
                 1.00000000        1.00000000        1.00000000\
        </attrib>\
     </group>\
  </particleset>\
  <particleset name=\"ion0\">\
     <group name=\"H\" size=\"2\" mass=\"1836.15\">\
        <parameter name=\"charge\"              >    1                     </parameter>\
        <parameter name=\"valence\"             >    1                     </parameter>\
        <parameter name=\"atomicnumber\"        >    1                     </parameter>\
        <parameter name=\"mass\"                >    1836.15               </parameter>\
        <attrib name=\"position\" datatype=\"posArray\" condition=\"0\">\
                 0.00000000        0.00000000        0.00000000\
                 1.889726135       1.889726135       1.889726135\
        </attrib>\
     </group>\
  </particleset>\
</tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root  = doc.getRoot();
  xmlNodePtr part1 = xmlFirstElementChild(root);
  xmlNodePtr part2 = xmlNextElementSibling(part1);
  xmlNodePtr part3 = xmlNextElementSibling(part2);

  // read lattice
  ParticleSet::ParticleLayout_t SimulationCell;
  LatticeParser lp(SimulationCell);
  lp.put(part1);
  SimulationCell.print(app_log(), 0);

  // read particle set
  ParticleSet ions, electrons;
  Tensor<int, 3> tmat; // assuming OHMMSDIM==3
  tmat(0, 0) = 1;
  tmat(1, 1) = 1;
  tmat(2, 2) = 1;

  // enforce global Lattice on ions and electrons
  ions.Lattice      = SimulationCell;
  electrons.Lattice = SimulationCell;

  XMLParticleParser parse_electrons(electrons, tmat);
  parse_electrons.put(part2);

  XMLParticleParser parse_ions(ions, tmat);
  parse_ions.put(part3);

  REQUIRE(electrons.getName() == "e");
  REQUIRE(ions.getName() == "ion0");
  REQUIRE(ions.SameMass);
  REQUIRE(electrons.SameMass);

  // calculate particle distances
  const int ei_tid = electrons.addTable(ions);
  electrons.update();

  // initialize orbital set
  const int cexpo = 1.0;
  const int ndim = 3;
  GaussianOrbitalSet* sposet = new GaussianOrbitalSet(electrons, ions, cexpo, ndim);
  electrons.update();
  int norb = sposet->getOrbitalSetSize();
  REQUIRE(norb==2);
  GaussianOrbitalSet::ValueVector_t p(norb);
  GaussianOrbitalSet::PosType dr(0);
  GaussianOrbitalSet::GradVector_t dp(norb);
  GaussianOrbitalSet::ValueVector_t d2p(norb);
  // finite difference
  const int ielec = 0;
  const double h = 0.00001;
  const double rel_tol = 1e-4;
  double v0, v1, lap0, lap1, dlap;
  GaussianOrbitalSet::PosType g0, g1;
  for (int ir=0;ir<10;ir++)
  {
    dr[0] = 0.1*ir;
    // numerical derivative
    v0 = std::exp(-cexpo*dot(dr, dr));
    lap0 = fd_laplacian(dr, h, sposet, electrons);
    g0 = fd_gradient(dr, h, sposet, electrons);
    // analytical derivative
    electrons.makeMove(ielec, dr);
    sposet->evaluateVGL(electrons, ielec, p, dp, d2p);
    v1 = p[0].real();
    for (int l=0;l<ndim;l++)
      g1[l] = (dp[ielec][l]*p[ielec]).real();
    lap1 = d2p[ielec].real()*p[ielec].real();
    // check value
    CHECK(v1 == Approx(v0));
    // check gradient
    for (int l=0;l<ndim;l++)
      CHECK(g1[l] == Approx(g0[l]));
    // check laplacian
    dlap = lap1-lap0;
    CHECK(std::abs(dlap/lap1) < rel_tol);
  }
}
} // qmcplusplus
