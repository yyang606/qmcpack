#include "ShortRangePairPotential.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
bool ShortRangePairPotential::put(xmlNodePtr cur)
{
  using std::sqrt;
  amplitute = 1.0;
  sigma = 0.2;
  OhmmsAttributeSet attrib;
  attrib.add(amplitute, "amplitute");
  attrib.add(sigma, "sigma");
  attrib.put(cur);
  return true;
}


bool ShortRangePairPotential::get(std::ostream& os) const
{
  os << "Short range pair potential" << std::endl;
  os << "  amplitute = " << amplitute << std::endl;
  os << "  sigma = " << sigma << std::endl;
  return true;
}


OperatorBase* ShortRangePairPotential::makeClone(ParticleSet& P, TrialWaveFunction& psi)
{
  return new ShortRangePairPotential(*this);
}


ShortRangePairPotential::Return_t ShortRangePairPotential::evaluate(ParticleSet& P)
{
  Value              = 0.0;
  return Value;
}

} // namespace qmcplusplus
