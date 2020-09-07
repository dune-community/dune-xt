// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <dune/xt/grid/grids.hh>

#include "grid-function_for_all_grids.hh"


PYBIND11_MODULE(_functions_interfaces_grid_function_2d, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.la");

  // All of these need to be there ...
  GridFunctionInterface_for_all_grids<boost::tuple<YASP_2D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                                   ,
                                                   ALU_2D_SIMPLEX_CONFORMING
#endif
                                                   >>::bind_interface(m);
  // ... before we start binding those.
  GridFunctionInterface_for_all_grids<boost::tuple<YASP_2D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                                   ,
                                                   ALU_2D_SIMPLEX_CONFORMING
#endif
                                                   >>::bind_combined(m);
} // PYBIND11_MODULE(...)
