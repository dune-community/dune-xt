// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <python/dune/xt/grid/grids.bindings.hh>

#include "element-function_for_all_grids.hh"


PYBIND11_MODULE(_functions_interfaces_element_function_2d, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.la");

  // All of these need to be there ...
  ElementFunctionSetInterface_for_all_grids<Dune::XT::Grid::bindings::Available2dGridTypes>::bind_interface(m);
  // ... before we start binding those.
  ElementFunctionSetInterface_for_all_grids<Dune::XT::Grid::bindings::Available2dGridTypes>::bind_combined(m);
} // PYBIND11_MODULE(...)
