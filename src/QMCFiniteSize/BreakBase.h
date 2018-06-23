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
}
#endif
