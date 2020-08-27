// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2019)
//   Tobias Leibner  (2020)

#include "config.h"

#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>

#include "function-interface.hh"


PYBIND11_MODULE(_functions_function_interface_3d, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");

  using namespace Dune::XT::Functions;

  using diff = CombinationType::difference;
  using sum = CombinationType::sum;
  using prod = CombinationType::product;

  auto i_3_1_1 = bind_FunctionInterface<3, 1, 1>(m);
  auto i_3_2_1 = bind_FunctionInterface<3, 2, 1>(m);
  auto i_3_2_2 = bind_FunctionInterface<3, 2, 2>(m);
  auto i_3_3_1 = bind_FunctionInterface<3, 3, 1>(m);
  auto i_3_3_3 = bind_FunctionInterface<3, 3, 3>(m);

  bind_combined_Function<3, diff, 1, 1, 1, 1>(m);
  addbind_FunctionInterface_combined_op<3, diff, 1, 1, 1, 1>(i_3_1_1);
  bind_combined_Function<3, diff, 2, 1, 2, 1>(m);
  addbind_FunctionInterface_combined_op<3, diff, 2, 1, 2, 1>(i_3_2_1);
  bind_combined_Function<3, diff, 2, 2, 2, 2>(m);
  addbind_FunctionInterface_combined_op<3, diff, 2, 2, 2, 2>(i_3_2_2);
  bind_combined_Function<3, diff, 3, 1, 3, 1>(m);
  addbind_FunctionInterface_combined_op<3, diff, 3, 1, 3, 1>(i_3_3_1);
  bind_combined_Function<3, diff, 3, 3, 3, 3>(m);
  addbind_FunctionInterface_combined_op<3, diff, 3, 3, 3, 3>(i_3_3_3);

  bind_combined_Function<3, sum, 1, 1, 1, 1>(m);
  addbind_FunctionInterface_combined_op<3, sum, 1, 1, 1, 1>(i_3_1_1);
  bind_combined_Function<3, sum, 2, 1, 2, 1>(m);
  addbind_FunctionInterface_combined_op<3, sum, 2, 1, 2, 1>(i_3_2_1);
  bind_combined_Function<3, sum, 2, 2, 2, 2>(m);
  addbind_FunctionInterface_combined_op<3, sum, 2, 2, 2, 2>(i_3_2_2);
  bind_combined_Function<3, sum, 3, 1, 3, 1>(m);
  addbind_FunctionInterface_combined_op<3, sum, 3, 1, 3, 1>(i_3_3_1);
  bind_combined_Function<3, sum, 3, 3, 3, 3>(m);
  addbind_FunctionInterface_combined_op<3, sum, 3, 3, 3, 3>(i_3_3_3);

  bind_combined_Function<3, prod, 1, 1, 1, 1>(m);
  addbind_FunctionInterface_combined_op<3, prod, 1, 1, 1, 1>(i_3_1_1);
  bind_combined_Function<3, prod, 1, 1, 2, 1>(m);
  addbind_FunctionInterface_combined_op<3, prod, 1, 1, 2, 1>(i_3_1_1);
  bind_combined_Function<3, prod, 1, 1, 2, 2>(m);
  addbind_FunctionInterface_combined_op<3, prod, 1, 1, 2, 2>(i_3_1_1);
  bind_combined_Function<3, prod, 1, 1, 3, 1>(m);
  addbind_FunctionInterface_combined_op<3, prod, 1, 1, 3, 1>(i_3_1_1);
  bind_combined_Function<3, prod, 1, 1, 3, 3>(m);
  addbind_FunctionInterface_combined_op<3, prod, 1, 1, 3, 3>(i_3_1_1);
} // PYBIND11_MODULE(...)
