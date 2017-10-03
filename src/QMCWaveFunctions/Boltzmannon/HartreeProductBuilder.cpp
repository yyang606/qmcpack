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
    // This function reads the xml node 'cur' and builds the HartreeProduct wavefunction.
    //  A pointer to the built wavefunction is passed to targetPsi->Z.
    // Returns:
    //   bool: success
    // Effect:
    //   add a term to targetPsi->Z
    
    // check that the builder is reading the correct xml node
    std::string tag;
    getNodeName(tag,cur);
    if (tag!=hartree_product_tag) APP_ABORT("Dev. Error: HartreeProductBuilder is intended to read <"<<OrbitalBuilderBase::hartree_product_tag<<"> not <"<<tag<<">. Disable abort with caution.");

    // read and validate inputs
    std::string orb_name    = "HartreeProduct";
    std::string group_name  = "";
    std::string sposet_name = "";
    std::string nameOpt,sposetOpt,groupOpt;
    OhmmsAttributeSet attrib;
    attrib.add(nameOpt,"name");
    attrib.add(groupOpt,"group");
    attrib.add(sposetOpt,"sposet");
    attrib.put(cur);
    if (nameOpt != "") orb_name = nameOpt;
    sposet_name = sposetOpt;
    group_name  = groupOpt;
    if (sposet_name=="") APP_ABORT("HartreeProduct requires sposet attribute.");
    if (group_name =="") APP_ABORT("HartreeProduct requires group attribute.");

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
    if (first==-1 || last==-1) APP_ABORT("failed to find particle species: "+group_name); 

    // find desired SPO set
    SPOSetBasePtr sposet = get_sposet(sposet_name); // get_sposet is from BasisSetFactory
    if (sposet==NULL) APP_ABORT("failed to find single-particle orbital set: "+sposet_name);

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

    return true; // this needs to be true for WaveFunctionPool::primaryPsi to be set
  }

}
