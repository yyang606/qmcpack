#include "DisorderPotential.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

DisorderPotential::Return_t DisorderPotential::evaluate(ParticleSet& P)
{
  value_ = 0.0;
  RealType expo;
  const auto& d_table = P.getDistTableAB(itab);
  for (int i=0;i<P.getTotalNum();i++)
  {
    const auto& dist = d_table.getDistRow(i);
    for (int j=0;j<nsite;j++)
    {
      expo = -dist[j]*dist[j]/(2*widths[j]*widths[j]);
      value_ += heights[j]*std::exp(expo);
    }
  }
  return value_;
}

bool DisorderPotential::put(xmlNodePtr cur)
{
  xmlNodePtr node = cur->xmlChildrenNode;
  while (node != NULL)
  {
    std::string cname((const char*)node->name);
    if (cname == "widths") putContent(widths, node);
    if (cname == "heights") putContent(heights, node);
    node = node->next;
  }
  return true;
}

bool DisorderPotential::get(std::ostream& os) const
{
  os << "External disorder potential" << std::endl;
  os << "  widths = ";
  for (int i=0;i<nsite;i++) os << widths[i] << " ";
  os << std::endl;
  os << "  heights = ";
  for (int i=0;i<nsite;i++) os << heights[i] << " ";
  os << std::endl;
  return true;
}

std::unique_ptr<OperatorBase> DisorderPotential::makeClone(ParticleSet& P, TrialWaveFunction& psi)
{
  return std::make_unique<DisorderPotential>(*this);
}

} // qmcplusplus
