#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Yubo "Paul" Yang, University of Illinois at Urbana-Champaign
#//
#// File created by: Yubo "Paul" Yang, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////


PROJECT(qmcfinitesize)

set(QFS_TOOLS spline3d)

FOREACH(p ${QFS_TOOLS})

  ADD_EXECUTABLE(${p} ${p}.cpp)
  TARGET_LINK_LIBRARIES(${p} qmcbase qmcutil)

  FOREACH(l ${QMC_UTIL_LIBS})
    TARGET_LINK_LIBRARIES(${p} ${l})
  ENDFOREACH(l ${QMC_UTIL_LIBS})

  if(ENABLE_TAU_PROFILE)
    TARGET_LINK_LIBRARIES(${p} tau)
  endif(ENABLE_TAU_PROFILE)

  IF(MPI_LIBRARY)
    TARGET_LINK_LIBRARIES(${p} ${MPI_LIBRARY})
  ENDIF(MPI_LIBRARY)

ENDFOREACH(p ${QFS_TOOLS})
