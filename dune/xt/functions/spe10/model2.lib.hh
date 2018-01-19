// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_FUNCTIONS_SPE10_MODEL2_LIB_HH
#define DUNE_XT_FUNCTIONS_SPE10_MODEL2_LIB_HH

#include "model2.hh"

#if DUNE_XT_WITH_PYTHON_BINDINGS


#define _DUNE_XT_FUNCTIONS_SPE10_MODEL2_LIB(_p, _G, _R, _r, _rC)                                                       \
  _p class Dune::XT::Functions::Spe10::                                                                                \
      Model2Function<typename _G::template Codim<0>::Entity, typename _G::ctype, _G::dimension, _R, _r, _rC>

#define DUNE_XT_FUNCTIONS_SPE10_MODEL2_LIB_YASP_3D(_p)                                                                 \
  _DUNE_XT_FUNCTIONS_SPE10_MODEL2_LIB(_p, YASP_3D_EQUIDISTANT_OFFSET, double, 3, 3)

DUNE_XT_FUNCTIONS_SPE10_MODEL2_LIB_YASP_3D(extern template);


#endif // DUNE_XT_WITH_PYTHON_BINDINGS


#endif // DUNE_XT_FUNCTIONS_SPE10_MODEL2_LIB_HH
