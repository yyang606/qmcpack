#include <QMCWaveFunctions/Boltzmannon/HartreeProductBuilder.h>
#include <QMCWaveFunctions/Boltzmannon/HartreeProduct.h>
#include <OhmmsData/AttributeSet.h>

namespace qmcplusplus
{

  // assign targetPtcl, targetPsi
  HartreeProductBuilder::HartreeProductBuilder(ParticleSet& psetIn ,TrialWaveFunction& psiIn):
    OrbitalBuilderBase(psetIn,psiIn) {};

  bool HartreeProductBuilder::put(xmlNodePtr cur)
  {
    // set defaults
    std::string group_name  = "p";
    std::string sposet_name = "spo_p";
    std::string orb_name    = "HartreeProduct";
    // use inputs to over-ride defaults
    std::string nameOpt,sposetOpt,groupOpt;
    OhmmsAttributeSet attrib;
    attrib.add(nameOpt,"name");
    attrib.add(groupOpt,"group");
    attrib.add(sposetOpt,"sposet");
    attrib.put(cur);
    if (nameOpt != "") orb_name = nameOpt;
    if (sposetOpt != "") sposet_name = sposetOpt;
    if (groupOpt != "") group_name = groupOpt;

    // find desired particles
    int first,last; first=last=-1;
    SpeciesSet& specset( targetPtcl.getSpeciesSet() );
    for (int ispec=0;ispec<specset.size();ispec++)
    {
      std::string species_name = specset.speciesName[ispec];
      if (species_name == group_name)
      {
        first = targetPtcl.first(ispec);
        last  = targetPtcl.last(ispec);
      }
    }

    if (first==-1 || last==-1) APP_ABORT("failed to find particle species "+group_name); 

    // find desired SPO set
    SPOSetBasePtr sposet = get_sposet(sposet_name); // get_sposet is from BasisSetFactory

    // make sure #particles = #spos
    if (last-first!=sposet->size())
    {
      std::ostringstream msg;
      msg << "number of particles " << last-first << " != number of orbitals " << sposet->size() << std::endl;
      APP_ABORT(msg.str());
    }

    // initialize wavefunction
    HartreeProduct *wf = new HartreeProduct(sposet,first,last);
    // add to wavefunction
    targetPsi.addOrbital(wf,orb_name);
  }

}
