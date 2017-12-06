//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"


#include "Utilities/RandomGenerator.h"
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Utilities/OhmmsInfo.h"
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/ConstantOrbital.h"
#include "QMCHamiltonians/BareKineticEnergy.h"
#include "Estimators/EstimatorManager.h"
#include "Estimators/TraceManager.h"
#include "QMCDrivers/VMC/VMCUpdatePbyP.h"
#include "QMCDrivers/DriftOperators.h"


#include <stdio.h>
#include <string>


using std::string;

namespace qmcplusplus
{

TEST_CASE("drift pbyp and node correction", "[drivers][drift]")
{
  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;
  OhmmsInfo("testlogfile");

  MCWalkerConfiguration elec;

  elec.setName("elec");
  elec.setBoundBox(false);
  std::vector<int> agroup(1);
  agroup[0] = 1;
  elec.create(agroup);

  double tau = 0.5;
  double mass= 1.0;
  double drift_max = std::sqrt(tau/mass);
  std::vector<double> massinv(1,1./mass);
  ParticleSet::ParticlePos_t drift(1);

  double xtot = 10.; // check a span of 10 (1/bohr) around 0.0
  double gradx = -xtot/2.;
  int nx = 100;
  double dx   = xtot/nx;

  //app_log() << " begin printing" << std::endl;
  for (int ix=0;ix<nx;ix++)
  {
    elec.G[0][0] = gradx;
    setScaledDriftPbyPandNodeCorr(tau,massinv,elec.G,drift);
    double dval = drift[0][0]; 

    // check Ceperley cap
    if (gradx>3.) REQUIRE( dval == Approx(drift_max) );
    if (gradx<-3.) REQUIRE( dval == Approx(-drift_max) );

    // check Umrigar rescale
    if (std::abs(gradx)<2.)  
    {
      double scale_factor = (-1.+std::sqrt(1.+2.*gradx*gradx*tau))/(gradx*gradx*tau);
      REQUIRE( dval == Approx(scale_factor*gradx*tau) );
    }

    //app_log() << gradx << " " << dval << std::endl;
    gradx += dx;
  }
  //app_log() << " end printing." << std::endl;
}
}

