// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2018)
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

#include "expression.hh"


PYBIND11_MODULE(_expression, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.functions._function_interface_1d");
  py::module::import("dune.xt.functions._function_interface_2d");
  py::module::import("dune.xt.functions._function_interface_3d");

  using namespace Dune::XT::Functions;

  bind_ExpressionFunction<1, 1, 1>(m);
  bind_ExpressionFunction<1, 2, 1>(m);
  bind_ExpressionFunction<1, 2, 2>(m);
  bind_ExpressionFunction<1, 3, 1>(m);
  bind_ExpressionFunction<1, 3, 3>(m);

  bind_ExpressionFunction<2, 1, 1>(m);
  bind_ExpressionFunction<2, 2, 1>(m);
  bind_ExpressionFunction<2, 2, 2>(m);
  bind_ExpressionFunction<2, 3, 1>(m);
  bind_ExpressionFunction<2, 3, 3>(m);

  bind_ExpressionFunction<3, 1, 1>(m);
  bind_ExpressionFunction<3, 2, 1>(m);
  bind_ExpressionFunction<3, 2, 2>(m);
  bind_ExpressionFunction<3, 3, 1>(m);
  bind_ExpressionFunction<3, 3, 3>(m);
}
