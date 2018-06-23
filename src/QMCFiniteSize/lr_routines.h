#ifndef QMCPLUSPLUS_LR_COMMON_H
#define QMCPLUSPLUS_LR_COMMON_H
#include "QMCFiniteSize/KspaceFunctions.hpp"
#include "QMCFiniteSize/EslerBreak.h"
#include "QMCFiniteSize/NatoliBreak.h"
#include "QMCFiniteSize/fsc_routines.h"
namespace qmcplusplus
{
BreakBase* create_break(
  Uniform3DGridLayout box,
  xmlXPathContextPtr doc
);
}
#endif
