// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_FUNCTIONS_CONSTANT_PBH
#define DUNE_XT_FUNCTIONS_CONSTANT_PBH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/type_traits.hh>

#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/fmatrix.hh>

#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/functions/constant.hh>

namespace Dune {
namespace XT {
namespace Functions {


/**
 * \note We would like to drop the d template paremter and use either of
\code
static const           size_t d = G::dimension;
static const constexpr size_t d = G::dimension;
\endcode
 *       but this triggers a bug in gcc-4.9, see e.g.: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=59937
 */
template <class G, size_t d, size_t r, size_t rC>
typename std::enable_if<Grid::is_grid<G>::value,
                        pybind11::class_<ConstantFunction<typename G::template Codim<0>::Entity,
                                                          typename G::ctype,
                                                          d,
                                                          double,
                                                          r,
                                                          rC>>>::type
bind_ConstantFunction(pybind11::module& m, const std::string& grid_id)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  typedef double R;

  typedef LocalizableFunctionInterface<E, D, d, R, r, rC> I;
  typedef ConstantFunction<E, D, d, R, r, rC> C;

  py::class_<C, I> c(
      m,
      std::string("ConstantFunction__" + grid_id + "_to_" + Common::to_string(r) + "x" + Common::to_string(rC)).c_str(),
      std::string("ConstantFunction__" + grid_id + "_to_" + Common::to_string(r) + "x" + Common::to_string(rC))
          .c_str());

  c.def(py::init<typename C::RangeType, std::string>(), "value"_a, "name"_a = C::static_id());
  c.def(py::init<typename C::RangeFieldType, std::string>(), "value"_a, "name"_a = C::static_id());

  c.def_property_readonly("static_id", [](const C& /*self*/) { return C::static_id(); });

  const std::string make_name = "make_constant_function_" + Common::to_string(r) + "x" + Common::to_string(rC);
  m.def(std::string(make_name).c_str(),
        [](const Grid::GridProvider<G, Grid::none_t>& /*grid*/,
           const typename C::RangeType& value,
           const std::string& name) { return C(value, name); },
        "grid_provider"_a,
        "value"_a,
        "name"_a = C::static_id());
  m.def(std::string(make_name).c_str(),
        [](const Grid::GridProvider<G, Grid::DD::SubdomainGrid<G>>& /*grid*/,
           const typename C::RangeType& value,
           const std::string& name) { return C(value, name); },
        "grid_provider"_a,
        "value"_a,
        "name"_a = C::static_id());
  m.def(std::string(make_name).c_str(),
        [](const Grid::GridProvider<G, Grid::none_t>& /*grid*/,
           const typename C::RangeFieldType& value,
           const std::string& name) { return C(value, name); },
        "grid_provider"_a,
        "value"_a,
        "name"_a = C::static_id());
  m.def(std::string(make_name).c_str(),
        [](const Grid::GridProvider<G, Grid::DD::SubdomainGrid<G>>& /*grid*/,
           const typename C::RangeFieldType& value,
           const std::string& name) { return C(value, name); },
        "grid_provider"_a,
        "value"_a,
        "name"_a = C::static_id());

  return c;
} // ... bind_ConstantFunction(...)


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_CONSTANT_PBH
