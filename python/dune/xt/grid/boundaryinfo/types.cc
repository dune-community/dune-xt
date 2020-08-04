// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2018)
//   Tobias Leibner  (2018 - 2019)

#include "config.h"

#include <sstream>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/operators.h>
#include <dune/xt/grid/boundaryinfo/interfaces.hh>
#include <dune/xt/grid/boundaryinfo/types.hh>

using namespace Dune::XT::Grid;


PYBIND11_MODULE(_grid_boundaryinfo_types, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");

  py::class_<BoundaryType>(m, "BoundaryType", "BoundaryType")
      .def(
          "__eq__",
          [](const BoundaryType& self, const BoundaryType& other) { return self == other; },
          py::is_operator())
      .def(
          "__neq__",
          [](const BoundaryType& self, const BoundaryType& other) { return self != other; },
          py::is_operator())
      .def_property_readonly("id", &BoundaryType::id);

#define BIND_(NAME)                                                                                                    \
  py::class_<NAME, BoundaryType>(m, #NAME, #NAME)                                                                      \
      .def(py::init([]() { return std::make_unique<NAME>(); }))                                                        \
      .def("__repr__", [](const NAME&) { return std::string(#NAME) + "()"; })

  BIND_(NoBoundary);
  BIND_(UnknownBoundary);
  BIND_(DirichletBoundary);
  BIND_(NeumannBoundary);
  BIND_(RobinBoundary);
  BIND_(ReflectingBoundary);
  BIND_(AbsorbingBoundary);
  BIND_(InflowBoundary);
  BIND_(OutflowBoundary);
  BIND_(InflowOutflowBoundary);
  BIND_(ImpermeableBoundary);

#undef BIND_
}
