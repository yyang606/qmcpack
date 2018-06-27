#ifndef QMCPLUSPLUS_BREAK_BASE_H
#define QMCPLUSPLUS_BREAK_BASE_H
#include <ostream>
#include <vector>
#include <functional>
#include "LongRange/LRBreakup.h"
#include "coulomb_types.h"
namespace qmcplusplus
{
DECLARE_COULOMB_TYPES
typedef Uniform3DGridLayout BoxType;
typedef mRealType           RealType;

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

// Base class of Break classes, which store short-range long-range breakup
//  information of periodic potential over the crystal lattice e.g. Ewald.
// This base class provides the interface of a "breaker". 
// Instantiate using BreakFactory::create_break
//
// Example:
//   Libxml2Document fxml;
//   fxml.parse("input.xml");
//   xmlXPathContextPtr doc = fxml.getXPathContext();
//   Uniform3DGridLayout box = create_box(doc);
//   BreakBase* breaker = break_factory.create_break(box, doc);
//   double fklr1 = breaker->evaluate_fklr(1.0);
class BreakBase
{
 public:
  BreakBase(BreakSpec params):params_(params), chisq_(-1){};
  ~BreakBase(){};
  RealType get_chisq(){return chisq_;};
  RealType get_rc(){return params_.rc;};
  RealType get_kc(){return params_.kc;};
  // do NOT overload operator<<, because instances will have pointer type
  void report(std::ostream& os)
  {
    os << std::endl;
    os << " long-range breakup" << std::endl;
    os << " ------------------" << std::endl;
    os << params_;
    std::streamsize ss = std::cout.precision();
    os << "  chi^2  = " << std::scientific << chisq_ << std::endl;
    os << std::setprecision(ss);
  }

  // goal in life: evaluate long-range potential at k
  virtual RealType evaluate_fklr(RealType k) = 0;

 protected:
  BreakSpec params_;
  RealType chisq_;
};
}
#endif
