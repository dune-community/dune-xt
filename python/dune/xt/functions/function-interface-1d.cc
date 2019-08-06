// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2018)
//   Tim Keil        (2018)

#include "config.h"

#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>

#include "function-interface.hh"


PYBIND11_MODULE(_function_interface_1d, m)
{
  namespace py = pybind11;

  Dune::XT::Common::bindings::add_initialization(m, "dune.xt.functions", "_function_interface_1d");

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");

  using namespace Dune::XT::Functions;

  const auto diff = CombinationType::difference;
  const auto sum = CombinationType::sum;
  const auto prod = CombinationType::product;

  auto i_1_1_1 = bind_FunctionInterface<1, 1, 1>(m);
  auto i_1_2_1 = bind_FunctionInterface<1, 2, 1>(m);
  auto i_1_2_2 = bind_FunctionInterface<1, 2, 2>(m);
  auto i_1_3_1 = bind_FunctionInterface<1, 3, 1>(m);
  auto i_1_3_3 = bind_FunctionInterface<1, 3, 3>(m);

  bind_combined_Function<1, diff, 1, 1, 1, 1>(m);
  addbind_FunctionInterface_combined_op<1, diff, 1, 1, 1, 1>(i_1_1_1);
  bind_combined_Function<1, diff, 2, 1, 2, 1>(m);
  addbind_FunctionInterface_combined_op<1, diff, 2, 1, 2, 1>(i_1_2_1);
  bind_combined_Function<1, diff, 2, 2, 2, 2>(m);
  addbind_FunctionInterface_combined_op<1, diff, 2, 2, 2, 2>(i_1_2_2);
  bind_combined_Function<1, diff, 3, 1, 3, 1>(m);
  addbind_FunctionInterface_combined_op<1, diff, 3, 1, 3, 1>(i_1_3_1);
  bind_combined_Function<1, diff, 3, 3, 3, 3>(m);
  addbind_FunctionInterface_combined_op<1, diff, 3, 3, 3, 3>(i_1_3_3);

  bind_combined_Function<1, sum, 1, 1, 1, 1>(m);
  addbind_FunctionInterface_combined_op<1, sum, 1, 1, 1, 1>(i_1_1_1);
  bind_combined_Function<1, sum, 2, 1, 2, 1>(m);
  addbind_FunctionInterface_combined_op<1, sum, 2, 1, 2, 1>(i_1_2_1);
  bind_combined_Function<1, sum, 2, 2, 2, 2>(m);
  addbind_FunctionInterface_combined_op<1, sum, 2, 2, 2, 2>(i_1_2_2);
  bind_combined_Function<1, sum, 3, 1, 3, 1>(m);
  addbind_FunctionInterface_combined_op<1, sum, 3, 1, 3, 1>(i_1_3_1);
  bind_combined_Function<1, sum, 3, 3, 3, 3>(m);
  addbind_FunctionInterface_combined_op<1, sum, 3, 3, 3, 3>(i_1_3_3);

  bind_combined_Function<1, prod, 1, 1, 1, 1>(m);
  addbind_FunctionInterface_combined_op<1, prod, 1, 1, 1, 1>(i_1_1_1);
  bind_combined_Function<1, prod, 1, 1, 2, 1>(m);
  addbind_FunctionInterface_combined_op<1, prod, 1, 1, 2, 1>(i_1_1_1);
  bind_combined_Function<1, prod, 1, 1, 2, 2>(m);
  addbind_FunctionInterface_combined_op<1, prod, 1, 1, 2, 2>(i_1_1_1);
  bind_combined_Function<1, prod, 1, 1, 3, 1>(m);
  addbind_FunctionInterface_combined_op<1, prod, 1, 1, 3, 1>(i_1_1_1);
  bind_combined_Function<1, prod, 1, 1, 3, 3>(m);
  addbind_FunctionInterface_combined_op<1, prod, 1, 1, 3, 3>(i_1_1_1);
} // PYBIND11_MODULE(...)
