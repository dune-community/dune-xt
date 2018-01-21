// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_FUNCTIONS_EXPRESSION_LIB_HH
#define DUNE_XT_FUNCTIONS_EXPRESSION_LIB_HH

#include "expression/default.hh"

#if DUNE_XT_WITH_PYTHON_BINDINGS


#define _DUNE_XT_FUNCTIONS_EXPRESSION_LIB(_p, _G, _R, _r, _rC)                                                         \
  _p class Dune::XT::Functions::                                                                                       \
      ExpressionFunction<typename _G::template Codim<0>::Entity, typename _G::ctype, _G::dimension, _R, _r, _rC>

#if HAVE_DUNE_ALUGRID
#define DUNE_XT_FUNCTIONS_EXPRESSION_LIB_ALU_2D(_p)                                                                    \
  _DUNE_XT_FUNCTIONS_EXPRESSION_LIB(_p, ALU_2D_SIMPLEX_CONFORMING, double, 1, 1);                                      \
  _DUNE_XT_FUNCTIONS_EXPRESSION_LIB(_p, ALU_2D_SIMPLEX_CONFORMING, double, 2, 2)
#else
#define DUNE_XT_FUNCTIONS_EXPRESSION_LIB_ALU_2D(_p)
#endif

#define DUNE_XT_FUNCTIONS_EXPRESSION_LIB_YASP_1D(_p)                                                                   \
  _DUNE_XT_FUNCTIONS_EXPRESSION_LIB(_p, YASP_1D_EQUIDISTANT_OFFSET, double, 1, 1)

#define DUNE_XT_FUNCTIONS_EXPRESSION_LIB_YASP_2D(_p)                                                                   \
  _DUNE_XT_FUNCTIONS_EXPRESSION_LIB(_p, YASP_2D_EQUIDISTANT_OFFSET, double, 1, 1);                                     \
  _DUNE_XT_FUNCTIONS_EXPRESSION_LIB(_p, YASP_2D_EQUIDISTANT_OFFSET, double, 2, 2)

#define DUNE_XT_FUNCTIONS_EXPRESSION_LIB_YASP_3D(_p)                                                                   \
  _DUNE_XT_FUNCTIONS_EXPRESSION_LIB(_p, YASP_3D_EQUIDISTANT_OFFSET, double, 1, 1);                                     \
  _DUNE_XT_FUNCTIONS_EXPRESSION_LIB(_p, YASP_3D_EQUIDISTANT_OFFSET, double, 3, 3)

#define DUNE_XT_FUNCTIONS_EXPRESSION_LIB_ALL_GRIDS(_p)                                                                 \
  DUNE_XT_FUNCTIONS_EXPRESSION_LIB_ALU_2D(_p);                                                                         \
  DUNE_XT_FUNCTIONS_EXPRESSION_LIB_YASP_1D(_p);                                                                        \
  DUNE_XT_FUNCTIONS_EXPRESSION_LIB_YASP_2D(_p);                                                                        \
  DUNE_XT_FUNCTIONS_EXPRESSION_LIB_YASP_3D(_p)

DUNE_XT_FUNCTIONS_EXPRESSION_LIB_ALL_GRIDS(extern template);


#endif // DUNE_XT_WITH_PYTHON_BINDINGS


#endif // DUNE_XT_FUNCTIONS_EXPRESSION_LIB_HH
