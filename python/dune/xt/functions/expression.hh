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
//   Tobias Leibner  (2018)

#ifndef DUNE_XT_FUNCTIONS_EXPRESSION_PBH
#define DUNE_XT_FUNCTIONS_EXPRESSION_PBH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/functions/expression.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace internal {


template <size_t d, size_t r, size_t rC>
struct addbind_ExpressionFunction_scalar_ctor
{
  template <class C>
  void operator()(C& /*c*/, const std::string& /*static_id*/)
  {}
};

template <size_t d, size_t r>
struct addbind_ExpressionFunction_scalar_ctor<d, r, 1>
{
  template <class C>
  void operator()(C& c, const std::string& static_id)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    c.def(py::init<const std::string,
                   const Common::FieldVector<std::string, r>,
                   const Common::FieldMatrix<std::string, r, d>,
                   const size_t,
                   const std::string>(),
          "variable"_a,
          "expressions"_a,
          "gradient_expressions"_a,
          "order"_a,
          "name"_a = static_id);

    c.def(py::init<const std::string, const Common::FieldVector<std::string, r>, const size_t, const std::string>(),
          "variable"_a,
          "expressions"_a,
          "order"_a,
          "name"_a = static_id);
  }
}; // struct addbind_ExpressionFunction_scalar_ctor<1, 1>


} // namespace internal


/**
 * \note We would like to drop the d template paremter and use either of
\code
static const           size_t d = G::dimension;
static const constexpr size_t d = G::dimension;
\endcode
 *       but this triggers a bug in gcc-4.9 and we thus need to use G::dimension
 *       everywhere: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=59937
 */
template <size_t d, size_t r>
pybind11::class_<ExpressionFunction<d, r, 3, double>> bind_ExpressionFunction(pybind11::module& m,
                                                                              std::integral_constant<int, 3>)
{
  namespace py = pybind11;
  using namespace pybind11::literals;
  const size_t rC = 3;
  using R = double;

  using I = FunctionInterface<d, r, rC, R>;
  using C = ExpressionFunction<d, r, rC, R>;

  const std::string c_name =
      "ExpressionFunction__" + Common::to_string(d) + "d_to_" + Common::to_string(r) + "x" + Common::to_string(rC);
  py::class_<C, I> c(m, std::string(c_name).c_str(), std::string(c_name).c_str());

  internal::addbind_ExpressionFunction_scalar_ctor<d, r, rC>()(c, C::static_id());
  c.def(py::init<const std::string,
                 const Common::FieldMatrix<std::string, r, rC>,
                 const Common::FieldVector<Common::FieldMatrix<std::string, rC, d>, r>,
                 const size_t,
                 const std::string>(),
        "variable"_a,
        "expressions"_a,
        "gradient_expressions"_a,
        "order"_a,
        "name"_a = C::static_id());

  c.def(py::init<const std::string, const Common::FieldMatrix<std::string, r, rC>, const size_t, const std::string>(),
        "variable"_a,
        "expressions"_a,
        "order"_a,
        "name"_a = C::static_id());

  c.def_property_readonly("static_id", [](const C& /*self*/) { return C::static_id(); });

  const std::string make_name = "make_expression_function_" + Common::to_string(r) + "x" + Common::to_string(rC);
  m.def(std::string(make_name).c_str(),
        [](const std::string& variable,
           const Common::FieldMatrix<std::string, r, rC>& expression,
           const Common::FieldVector<Common::FieldMatrix<std::string, rC, d>, r>& gradient_expressions,
           const size_t& order,
           const std::string& name) { return C(variable, expression, gradient_expressions, order, name); },
        "variable"_a,
        "expression"_a,
        "gradient_expressions"_a,
        "order"_a,
        "name"_a = C::static_id());

  m.def(std::string(make_name).c_str(),
        [](const std::string& variable,
           const Common::FieldMatrix<std::string, r, rC>& expression,
           const size_t& order,
           const std::string& name) { return C(variable, expression, order, name); },
        "variable"_a,
        "expression"_a,
        "order"_a,
        "name"_a = C::static_id());

  return c;
} // ... bind_ExpressionFunction(...)

template <size_t d, size_t r>
pybind11::class_<ExpressionFunction<d, r, 2, double>> bind_ExpressionFunction(pybind11::module& m,
                                                                              std::integral_constant<int, 2>)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  const size_t rC = 2;
  using R = double;

  using I = FunctionInterface<d, r, rC, R>;
  using C = ExpressionFunction<d, r, rC, R>;

  const std::string c_name =
      "ExpressionFunction__" + Common::to_string(d) + "d_to_" + Common::to_string(r) + "x" + Common::to_string(rC);
  py::class_<C, I> c(m, std::string(c_name).c_str(), std::string(c_name).c_str());

  internal::addbind_ExpressionFunction_scalar_ctor<d, r, rC>()(c, C::static_id());
  c.def(py::init<const std::string,
                 const Common::FieldMatrix<std::string, r, rC>,
                 const Common::FieldVector<Common::FieldMatrix<std::string, rC, d>, r>,
                 const size_t,
                 const std::string>(),
        "variable"_a,
        "expressions"_a,
        "gradient_expressions"_a,
        "order"_a,
        "name"_a = C::static_id());

  c.def(py::init<const std::string, const Common::FieldMatrix<std::string, r, rC>, const size_t, const std::string>(),
        "variable"_a,
        "expressions"_a,
        "order"_a,
        "name"_a = C::static_id());

  c.def_property_readonly("static_id", [](const C& /*self*/) { return C::static_id(); });

  const std::string make_name = "make_expression_function_" + Common::to_string(r) + "x" + Common::to_string(rC);
  m.def(std::string(make_name).c_str(),
        [](const std::string& variable,
           const Common::FieldMatrix<std::string, r, rC>& expression,
           const Common::FieldVector<Common::FieldMatrix<std::string, rC, d>, r>& gradient_expressions,
           const size_t& order,
           const std::string& name) { return C(variable, expression, gradient_expressions, order, name); },
        "variable"_a,
        "expression"_a,
        "gradient_expressions"_a,
        "order"_a,
        "name"_a = C::static_id());

  m.def(std::string(make_name).c_str(),
        [](const std::string& variable,
           const Common::FieldMatrix<std::string, r, rC>& expression,
           const size_t& order,
           const std::string& name) { return C(variable, expression, order, name); },
        "variable"_a,
        "expression"_a,
        "order"_a,
        "name"_a = C::static_id());

  return c;
} // ... bind_ExpressionFunction(...)

template <size_t d, size_t r>
pybind11::class_<ExpressionFunction<d, r, 1, double>> bind_ExpressionFunction(pybind11::module& m,
                                                                              std::integral_constant<int, 1>)
{
  namespace py = pybind11;
  using namespace pybind11::literals;
  const size_t rC = 1;
  using R = double;

  using I = FunctionInterface<d, r, rC, R>;
  using C = ExpressionFunction<d, r, rC, R>;

  const std::string c_name =
      "ExpressionFunction__" + Common::to_string(d) + "d_to_" + Common::to_string(r) + "x" + Common::to_string(1);
  py::class_<C, I> c(m, std::string(c_name).c_str(), std::string(c_name).c_str());

  internal::addbind_ExpressionFunction_scalar_ctor<d, r, 1>()(c, C::static_id());
  c.def(py::init<const std::string,
                 const Common::FieldVector<std::string, r>,
                 const Common::FieldMatrix<std::string, r, d>,
                 const size_t,
                 const std::string>(),
        "variable"_a,
        "expressions"_a,
        "gradient_expressions"_a,
        "order"_a,
        "name"_a = C::static_id());

  c.def(py::init<const std::string, const Common::FieldVector<std::string, r>, const size_t, const std::string>(),
        "variable"_a,
        "expressions"_a,
        "order"_a,
        "name"_a = C::static_id());

  c.def_property_readonly("static_id", [](const C& /*self*/) { return C::static_id(); });

  const std::string make_name = "make_expression_function_" + Common::to_string(r) + "x" + Common::to_string(1);
  m.def(std::string(make_name).c_str(),
        [](const std::string& variable,
           const Common::FieldVector<std::string, r>& expression,
           const Common::FieldMatrix<std::string, r, d>& gradient_expressions,
           const size_t& order,
           const std::string& name) { return C(variable, expression, gradient_expressions, order, name); },
        "variable"_a,
        "expression"_a,
        "gradient_expressions"_a,
        "order"_a,
        "name"_a = C::static_id());

  m.def(std::string(make_name).c_str(),
        [](const std::string& variable,
           const Common::FieldVector<std::string, r>& expression,
           const size_t& order,
           const std::string& name) { return C(variable, expression, order, name); },
        "variable"_a,
        "expression"_a,
        "order"_a,
        "name"_a = C::static_id());

  return c;
}


template <size_t d, size_t r, size_t rC>
pybind11::class_<ExpressionFunction<d, r, rC, double>> bind_ExpressionFunction(pybind11::module& m)
{
  return bind_ExpressionFunction<d, r>(m, std::integral_constant<int, rC>());
} // ... bind_ExpressionFunction(...)


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_EXPRESSION_PBH
