#include "ShortRangePairPotential.h"
#include "OhmmsData/AttributeSet.h"
#include "Particle/DistanceTableData.h"

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
  const DistanceTableData& d_aa(P.getDistTable(id_daa));
  RealType s2=2*sigma*sigma;
  for (size_t i=1; i<P.getTotalNum(); i++)
  {
    const auto& dist = d_aa.getDistRow(i);
    for (size_t j=0; j<i; j++)
    {
      auto rij = dist[j];
      Value += amplitute*std::exp(-rij*rij/s2);
    }
  }
  return Value;
}

} // namespace qmcplusplus
