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

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/expression.hh>

#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/fmatrix.hh>
#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/grid/traits.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>

namespace Dune::XT::Functions::bindings {


template <size_t d, size_t r = 1, size_t rC = 1, class R = double>
class ExpressionFunction
{
  using type = Functions::ExpressionFunction<d, r, rC, R>;
  using base_type = Functions::FunctionInterface<d, r, rC, R>;
  using bound_type = pybind11::class_<type, base_type>;

  template <size_t r_ = r, size_t rC_ = rC, bool anything = true> // the matrix-valued case
  struct addbind
  {
    static void ctor(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      c.def(py::init<const std::string&,
                     const Common::FieldMatrix<std::string, r, rC>&,
                     const Common::FieldVector<Common::FieldMatrix<std::string, rC, d>, r>&,
                     const size_t,
                     const std::string>(),
            "variable"_a,
            "expressions"_a,
            "gradient_expressions"_a,
            "order"_a,
            "name"_a = type::static_id());
      c.def(py::init<const std::string&,
                     const Common::FieldMatrix<std::string, r, rC>&,
                     const size_t,
                     const std::string>(),
            "variable"_a,
            "expressions"_a,
            "order"_a,
            "name"_a = type::static_id());
    } // ... ctor(...)

    static void factory(pybind11::module& m, const std::string& ClassId)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      m.def(
          ClassId.c_str(),
          [](Grid::bindings::Dimension<d> /*dim_domain*/,
             const std::string& variable,
             const Common::FieldMatrix<std::string, r, rC>& expressions,
             const Common::FieldVector<Common::FieldMatrix<std::string, rC, d>, r>& gradient_expressions,
             const size_t order,
             const std::string& name) { return type(variable, expressions, gradient_expressions, order, name); },
          "dim_domain"_a,
          "variable"_a,
          "expressions"_a,
          "gradient_expressions"_a,
          "order"_a,
          "name"_a = type::static_id());
      m.def(
          ClassId.c_str(),
          [](Grid::bindings::Dimension<d> /*dim_domain*/,
             const std::string& variable,
             const Common::FieldMatrix<std::string, r, rC>& expressions,
             const size_t order,
             const std::string& name) { return type(variable, expressions, order, name); },
          "dim_domain"_a,
          "variable"_a,
          "expressions"_a,
          "order"_a,
          "name"_a = type::static_id());
    } // ... factory(...)
  }; // struct addbind, the matrix-valued case

  template <size_t r_, bool anything> // the vector-valued case
  struct addbind<r_, 1, anything>
  {
    static void ctor(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      c.def(py::init<const std::string&,
                     const Common::FieldVector<std::string, r>&,
                     const Common::FieldMatrix<std::string, r, d>&,
                     const size_t,
                     const std::string>(),
            "variable"_a,
            "expressions"_a,
            "gradient_expressions"_a,
            "order"_a,
            "name"_a = type::static_id());
      c.def(py::init<const std::string&, const Common::FieldVector<std::string, r>&, const size_t, const std::string>(),
            "variable"_a,
            "expressions"_a,
            "order"_a,
            "name"_a = type::static_id());
    } // ... ctor(...)

    static void factory(pybind11::module& m, const std::string& ClassId)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      m.def(
          ClassId.c_str(),
          [](Grid::bindings::Dimension<d> /*dim_domain*/,
             const std::string& variable,
             const Common::FieldVector<std::string, r>& expressions,
             const Common::FieldMatrix<std::string, r, d>& gradient_expressions,
             const size_t order,
             const std::string& name) { return type(variable, expressions, gradient_expressions, order, name); },
          "dim_domain"_a,
          "variable"_a,
          "expressions"_a,
          "gradient_expressions"_a,
          "order"_a,
          "name"_a = type::static_id());
      m.def(
          ClassId.c_str(),
          [](Grid::bindings::Dimension<d> /*dim_domain*/,
             const std::string& variable,
             const Common::FieldVector<std::string, r>& expressions,
             const size_t order,
             const std::string& name) { return type(variable, expressions, order, name); },
          "dim_domain"_a,
          "variable"_a,
          "expressions"_a,
          "order"_a,
          "name"_a = type::static_id());
    } // ... factory(...)
  }; // struct addbind, the vector-valued case

  template <bool anything>
  struct addbind<1, 1, anything> // the scalar case
  {
    static void ctor(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      c.def(py::init([](const std::string& variable,
                        const std::string& expression,
                        const Common::FieldVector<std::string, d>& gradient_expressions,
                        const size_t order,
                        const std::string& name) {
              Common::FieldMatrix<std::string, 1, d> gradient_expressions_mat;
              gradient_expressions_mat[0] = gradient_expressions;
              return new type(
                  variable, Common::FieldVector<std::string, 1>(expression), gradient_expressions_mat, order, name);
            }),
            "variable"_a,
            "expression"_a,
            "gradient_expressions"_a,
            "order"_a,
            "name"_a = type::static_id());
      c.def(py::init([](const std::string& variable,
                        const std::string& expression,
                        const size_t order,
                        const std::string& name) {
              return new type(variable, Common::FieldVector<std::string, 1>(expression), order, name);
            }),
            "variable"_a,
            "expression"_a,
            "order"_a,
            "name"_a = type::static_id());
    } // ... ctor(...)

    static void factory(pybind11::module& m, const std::string& ClassId)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      m.def(
          ClassId.c_str(),
          [](Grid::bindings::Dimension<d> /*dim_domain*/,
             const std::string& variable,
             const std::string& expression,
             const Common::FieldVector<std::string, d>& gradient_expressions,
             const size_t order,
             const std::string& name) {
            Common::FieldMatrix<std::string, 1, d> gradient_expressions_mat;
            gradient_expressions_mat[0] = gradient_expressions;
            return new type(
                variable, Common::FieldVector<std::string, 1>(expression), gradient_expressions_mat, order, name);
          },
          "dim_domain"_a,
          "variable"_a,
          "expression"_a,
          "gradient_expressions"_a,
          "order"_a,
          "name"_a = type::static_id());
      m.def(
          ClassId.c_str(),
          [](Grid::bindings::Dimension<d> /*dim_domain*/,
             const std::string& variable,
             const std::string& expression,
             const size_t order,
             const std::string& name) { return new type(variable, expression, order, name); },
          "dim_domain"_a,
          "variable"_a,
          "expression"_a,
          "order"_a,
          "name"_a = type::static_id());
    } // ... factory(...)
  }; // struct addbind, the scalar case

public:
  static bound_type bind(pybind11::module& m, const std::string& class_id = "expression_function")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id + "_" + Common::to_string(d) + "_to_" + Common::to_string(r);
    if (rC > 1)
      class_name += "x" + Common::to_string(rC);
    class_name += "d";
    const auto ClassName = Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), XT::Common::to_camel_case(class_id).c_str());
    addbind<>::ctor(c);

    addbind<>::factory(m, Common::to_camel_case(class_id));

    return c;
  }
}; // class ExpressionFunction


} // namespace Dune::XT::Functions::bindings


PYBIND11_MODULE(_functions_expression, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.functions._functions_function_interface_1d");
  py::module::import("dune.xt.functions._functions_function_interface_2d");
  py::module::import("dune.xt.functions._functions_function_interface_3d");

  Dune::XT::Functions::bindings::ExpressionFunction<1, 1, 1>::bind(m);
  Dune::XT::Functions::bindings::ExpressionFunction<1, 2, 1>::bind(m);
  Dune::XT::Functions::bindings::ExpressionFunction<1, 2, 2>::bind(m);
  Dune::XT::Functions::bindings::ExpressionFunction<1, 3, 1>::bind(m);
  Dune::XT::Functions::bindings::ExpressionFunction<1, 3, 3>::bind(m);

  Dune::XT::Functions::bindings::ExpressionFunction<2, 1, 1>::bind(m);
  Dune::XT::Functions::bindings::ExpressionFunction<2, 2, 1>::bind(m);
  Dune::XT::Functions::bindings::ExpressionFunction<2, 2, 2>::bind(m);
  Dune::XT::Functions::bindings::ExpressionFunction<2, 3, 1>::bind(m);
  Dune::XT::Functions::bindings::ExpressionFunction<2, 3, 3>::bind(m);

  Dune::XT::Functions::bindings::ExpressionFunction<3, 1, 1>::bind(m);
  Dune::XT::Functions::bindings::ExpressionFunction<3, 2, 1>::bind(m);
  Dune::XT::Functions::bindings::ExpressionFunction<3, 2, 2>::bind(m);
  Dune::XT::Functions::bindings::ExpressionFunction<3, 3, 1>::bind(m);
  Dune::XT::Functions::bindings::ExpressionFunction<3, 3, 3>::bind(m);
} // PYBIND11_MODULE(...)
