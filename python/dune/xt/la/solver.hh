// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_LA_SOLVER_PBH
#define DUNE_XT_LA_SOLVER_PBH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/operators.h>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/la/container.bindings.hh>

#include <dune/xt/la/container.hh>
#include <dune/xt/la/type_traits.hh>
#include <dune/xt/la/solver.hh>

namespace Dune {
namespace XT {
namespace LA {


template <class M, class V = typename Container<typename M::ScalarType, M::vector_type>::VectorType>
auto bind_Solver(pybind11::module& m)
{
  static_assert(is_matrix<M>::value);
  using C = Solver<M>;

  namespace py = pybind11;
  using namespace pybind11::literals;

  const auto ClassName = Common::to_camel_case(bindings::container_name<M>::value() + "_solver");

  py::class_<C> c(m, ClassName.c_str(), ClassName.c_str());

  c.def_static("types", &C::types);
  c.def_static("options", &C::options);

  c.def(py::init<M>());

  c.def(
      "apply", [](const C& self, const V& rhs, V& solution) { self.apply(rhs, solution); }, "rhs"_a, "solution"_a);
  c.def(
      "apply",
      [](const C& self, const V& rhs, V& solution, const std::string& type) { self.apply(rhs, solution, type); },
      "rhs"_a,
      "solution"_a,
      "type"_a);
  c.def(
      "apply",
      [](const C& self, const V& rhs, V& solution, const Common::Configuration& options) {
        self.apply(rhs, solution, options);
      },
      "rhs"_a,
      "solution"_a,
      "options"_a);

  m.def(
      "make_solver", [](const M& matrix) { return C(matrix); }, pybind11::keep_alive<0, 1>());

  return c;
} // ... bind_Solver(...)


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_SOLVER_PBH
