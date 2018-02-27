// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_FUNCTIONS_EXPRESSION_PBH
#define DUNE_XT_FUNCTIONS_EXPRESSION_PBH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/functions/expression.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace internal {


template <size_t r, size_t rC>
struct addbind_ExpressionFunction_scalar_ctor
{
  template <class C>
  void operator()(C& /*c*/, const std::string& /*static_id*/)
  {
  }
};

template <size_t r>
struct addbind_ExpressionFunction_scalar_ctor<r, 1>
{
  template <class C>
  void operator()(C& c, const std::string& static_id)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    c.def(py::init<const std::string,
                   const std::vector<std::string>,
                   const size_t,
                   const std::string,
                   const std::vector<std::vector<std::string>>>(),
          "variable"_a,
          "expressions"_a,
          "order"_a,
          "name"_a = static_id,
          "gradient_expressions"_a = std::vector<std::vector<std::string>>());
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
template <class G, size_t d, size_t r, size_t rC>
typename std::enable_if<Grid::is_grid<G>::value,
                        pybind11::class_<ExpressionFunction<typename G::template Codim<0>::Entity,
                                                            typename G::ctype,
                                                            G::dimension,
                                                            double,
                                                            r,
                                                            rC>>>::type
bind_ExpressionFunction(pybind11::module& m, const std::string& grid_id)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  typedef double R;

  typedef LocalizableFunctionInterface<E, D, d, R, r, rC> I;
  typedef ExpressionFunction<E, D, d, R, r, rC> C;

  const std::string c_name =
      "ExpressionFunction__" + grid_id + "_to_" + Common::to_string(r) + "x" + Common::to_string(rC);
  py::class_<C, I> c(m, std::string(c_name).c_str(), std::string(c_name).c_str());

  internal::addbind_ExpressionFunction_scalar_ctor<r, rC>()(c, C::static_id());
  c.def(
      py::init<const std::string, const std::string, const size_t, const std::string, const std::vector<std::string>>(),
      "variable"_a,
      "expression"_a,
      "order"_a,
      "name"_a = C::static_id(),
      "gradient_expressions"_a = std::vector<std::string>());
  c.def(py::init<const std::string,
                 const std::vector<std::vector<std::string>>,
                 const size_t,
                 const std::string,
                 const std::vector<std::vector<std::vector<std::string>>>>(),
        "variable"_a,
        "expressions"_a,
        "order"_a,
        "name"_a = C::static_id(),
        "gradient_expressions"_a = std::vector<std::vector<std::string>>());

  c.def_property_readonly("static_id", [](const C& /*self*/) { return C::static_id(); });

  const std::string make_name = "make_expression_function_" + Common::to_string(r) + "x" + Common::to_string(rC);
  m.def(std::string(make_name).c_str(),
        [](const Grid::GridProvider<G>& /*grid*/,
           const std::string& variable,
           const std::string& expression,
           const size_t& order,
           const std::string& name,
           const std::vector<std::string>& gradient_expressions) {
          return C(variable, expression, order, name, gradient_expressions);
        },
        "grid_provider"_a,
        "variable"_a,
        "expression"_a,
        "order"_a,
        "name"_a = C::static_id(),
        "gradient_expressions"_a = std::vector<std::string>());
  m.def(std::string(make_name).c_str(),
        [](const Grid::GridProvider<G, Grid::DD::SubdomainGrid<G>>& /*grid*/,
           const std::string& variable,
           const std::string& expression,
           const size_t& order,
           const std::string& name,
           const std::vector<std::string>& gradient_expressions) {
          return C(variable, expression, order, name, gradient_expressions);
        },
        "grid_provider"_a,
        "variable"_a,
        "expression"_a,
        "order"_a,
        "name"_a = C::static_id(),
        "gradient_expressions"_a = std::vector<std::string>());
  m.def(std::string(make_name).c_str(),
        [](const Grid::GridProvider<G>& /*grid*/,
           const std::string& variable,
           const std::vector<std::vector<std::string>>& expression,
           const size_t& order,
           const std::string& name,
           const std::vector<std::vector<std::vector<std::string>>>& gradient_expressions) {
          return C(variable, expression, order, name, gradient_expressions);
        },
        "grid_provider"_a,
        "variable"_a,
        "expression"_a,
        "order"_a,
        "name"_a = C::static_id(),
        "gradient_expressions"_a = std::vector<std::vector<std::vector<std::string>>>());
  m.def(std::string(make_name).c_str(),
        [](const Grid::GridProvider<G, Grid::DD::SubdomainGrid<G>>& /*grid*/,
           const std::string& variable,
           const std::vector<std::vector<std::string>>& expression,
           const size_t& order,
           const std::string& name,
           const std::vector<std::vector<std::vector<std::string>>>& gradient_expressions) {
          return C(variable, expression, order, name, gradient_expressions);
        },
        "grid_provider"_a,
        "variable"_a,
        "expression"_a,
        "order"_a,
        "name"_a = C::static_id(),
        "gradient_expressions"_a = std::vector<std::vector<std::vector<std::string>>>());

  return c;
} // ... bind_ExpressionFunction(...)


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_EXPRESSION_PBH
