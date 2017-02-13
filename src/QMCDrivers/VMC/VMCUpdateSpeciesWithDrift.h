//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Yubo Yang, yyang173@illinois.edu, University of Illinois at Urbana-Champaign 
//
// File created by: Yubo Yang, yyang173@illinois.edu, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_VMC_UPDATESPECIESDRIFT_H
#define QMCPLUSPLUS_VMC_UPDATESPECIESDRIFT_H
#include "QMCDrivers/QMCUpdateBase.h"

namespace qmcplusplus
{

/** @ingroup QMCDrivers  SpeciesBySpecies
 *@brief Implements the VMC algorithm using all-particle move for each species with the drift equation.
 */
class VMCUpdateSpeciesWithDrift: public QMCUpdateBase
{
public:
  /// Constructor.
  VMCUpdateSpeciesWithDrift(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                        QMCHamiltonian& h, RandomGenerator_t& rg);

  ~VMCUpdateSpeciesWithDrift();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
  bool put(xmlNodePtr cur);

private:
  /// Copy Constructor (disabled)
  VMCUpdateSpeciesWithDrift(const VMCUpdateSpeciesWithDrift& a): QMCUpdateBase(a),tspecies(a.tspecies) { }
  /// Copy operator (disabled).
  VMCUpdateSpeciesWithDrift& operator=(const VMCUpdateSpeciesWithDrift&)
  {
    return *this;
  }
  std::vector<int> move_interval;
  std::vector<int> nrest;
  bool move_species_all(int ig, Walker_t& thisWalker);
  SpeciesSet& tspecies;
};

}

#endif
/***************************************************************************
 * $RCSfile: VMCUpdateSpeciesWithDrift.h,v $   $Author: yyang173 $
 * $Revision: 7223 $   $Date: 2016/11/14 14:29:40 $
 * $Id: VMCUpdateAll.h,v 1.5 2016/11/14 14:29:40 yyang173 Exp $
 ***************************************************************************/
