#ifndef QMCPLUSPLUS_BREAK_BASE_H
#define QMCPLUSPLUS_BREAK_BASE_H
#include <ostream>
#include <vector>
#include <functional>
#include "LongRange/LRBreakup.h"
namespace qmcplusplus
{
typedef Uniform3DGridLayout                  BoxType;
typedef QMCTraits::RealType                 RealType;

struct BreakSpec
{
  std::string basis;
  int nknot;
  double rc, kc, kcut, kmax;
};
inline std::ostream& operator<<(std::ostream& os, const BreakSpec& params)
{
  os << "  basis  = " << params.basis << std::endl;
  os << "  nknot  = " << params.nknot << std::endl;
  os << "  rc     = " << params.rc    << std::endl;
  os << "  kc     = " << params.kc    << std::endl;
  os << "  kcut   = " << params.kcut  << std::endl;
  os << "  kmax   = " << params.kmax  << std::endl;
  return os;
}

// Prototype of Break classes, which store short-range long-range breakup
//  information of periodic potential over the crystal lattice e.g. Ewald.
// This base class provides the interface of a "breaker". 
// Instantiate using BreakPrototypeFactory::create_break
//
// Example:
//   BreakBase* breaker = break_factory.create_break(doc);
//   double fklr1 = breaker->evaluate_fklr(1.0);
class BreakBase
{
 public:
  BreakBase(BreakSpec params):params_(params), chisq_(-1){};
  ~BreakBase(){};
  RealType get_chisq(){return chisq_;};
  RealType get_rc(){return params_.rc;};
  RealType get_kc(){return params_.kc;};
  virtual void report(std::ostream& os) = 0;

  // goal in life: evaluate long-range potential at k
  virtual RealType evaluate_fklr(RealType k) = 0;
 protected:
  BreakSpec params_;
  RealType chisq_;
};
}
#endif
