#ifndef QMCPLUSPLUS_BREAK_FACTORY_H
#define QMCPLUSPLUS_BREAK_FACTORY_H
#include "QMCFiniteSize/KspaceFunctions.h"
#include "BreakBase.h"
#include "NatoliBreak.h"
#include "EslerBreak.h"
#include "QMCFiniteSize/fsc_routines.h"
namespace qmcplusplus
{
class BreakFactory
{
 public:
  BreakBase* create_break(
   Uniform3DGridLayout box,
   xmlXPathContextPtr doc
  );
 private:
  RealType default_kc(Uniform3DGridLayout box);
};
} // qmcplusplus
#endif
