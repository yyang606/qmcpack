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

#include "GaussianOrbitalBuilder.h"
#include "GaussianOrbital.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

GaussianOrbitalBuilder::GaussianOrbitalBuilder(ParticleSet& els, const PSetMap& psets, Communicate* comm, xmlNodePtr cur)
  : SPOSetBuilder("Gaussian", comm),
    targetPtcl(els),
    particleSets(psets)
{}

GaussianOrbitalBuilder::~GaussianOrbitalBuilder(){}

std::unique_ptr<SPOSet> GaussianOrbitalBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  std::string spo_object_name="Gaussian";
  RealType cexpo=0.1;
  std::vector<RealType> cexpos;
  std::string sourceName="ion0";
  TinyVector<int, OHMMS_DIM> pbc_images=0;

  OhmmsAttributeSet attrib;
  attrib.add(cexpo, "c");
  attrib.add(sourceName, "source");
  attrib.add(spo_object_name, "name");
  attrib.add(pbc_images, "PBCimages");
  attrib.put(cur);

  ParticleSet* sourcePtcl;
  auto pit(particleSets.find(sourceName));
  if (pit == particleSets.end())
    myComm->barrier_and_abort("GaussianOrbital needs the source particleset");
  else
    sourcePtcl = pit->second.get();
  cexpos.resize(sourcePtcl->getTotalNum());

  xmlNodePtr element = cur->xmlChildrenNode;
  while (element != NULL)
  {
    std::string ename((const char*)element->name);
    if (ename == "parameter")
    {
      const std::string name(getXMLAttributeValue(element, "name"));
      if (name == "cexpos")
      {
        putContent(cexpos, element);
      }
    }
    element = element->next;
  }

  auto sposet = std::make_unique<GaussianOrbital>(spo_object_name, targetPtcl, *sourcePtcl, cexpos, pbc_images);
  sposet->report("  ");
  return sposet;
}

} // qmcplusplus
