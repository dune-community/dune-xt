// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2018)

#ifndef PYTHON_DUNE_XT_FUNCTIONS_INDICATOR_HH
#define PYTHON_DUNE_XT_FUNCTIONS_INDICATOR_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/type_traits.hh>

#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/fmatrix.hh>

#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/functions/indicator.hh>

namespace Dune {
namespace XT {
namespace Functions {


template <class G, size_t d, size_t r, size_t rC>
typename std::enable_if<Grid::is_grid<G>::value,
                        pybind11::class_<IndicatorGridFunction<typename G::template Codim<0>::Entity, r, rC, double>>>::
    type
    bind_IndicatorGridFunction(pybind11::module& m, const std::string& grid_id)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  typedef double R;
  typedef GridFunctionInterface<E, r, rC, R> I;
  typedef IndicatorGridFunction<E, r, rC, R> C;
  using DomainType = typename C::DomainType;
  using RangeType = typename C::RangeType;
  using CornerVector = std::vector<std::tuple<DomainType, DomainType, RangeType>>;
  using IntervalVector = std::vector<std::pair<Common::FieldMatrix<D, d, 2>, RangeType>>;

  const std::string classname =
      std::string("IndicatorGridFunction__" + grid_id + "_to_" + Common::to_string(r) + "x" + Common::to_string(rC));
  py::class_<C, I> c(m, classname.c_str(), classname.c_str());
  c.def(py::init<CornerVector, std::string>(), "corner_and_value_pairs"_a, "name"_a = C::static_id());
  c.def(py::init<IntervalVector, std::string>(), "interval_and_value_pairs"_a, "name"_a = C::static_id());
  c.def_property_readonly("static_id", [](const C& /*self*/) { return C::static_id(); });

  const std::string make_name = "make_indicator_function_" + Common::to_string(r) + "x" + Common::to_string(rC);
  m.def(std::string(make_name).c_str(),
        [](const Grid::GridProvider<G, Grid::none_t>& /*grid*/, const CornerVector& values, const std::string& name) {
          return C(values, name);
        },
        "grid_provider"_a,
        "values"_a,
        "name"_a = C::static_id());
  m.def(std::string(make_name).c_str(),
        [](const Grid::GridProvider<G, Grid::DD::SubdomainGrid<G>>& /*grid*/,
           const IntervalVector& value,
           const std::string& name) { return C(value, name); },
        "grid_provider"_a,
        "values"_a,
        "name"_a = C::static_id());
  m.def(std::string(make_name).c_str(),
        [](const Grid::GridProvider<G, Grid::none_t>& /*grid*/, const CornerVector& values, const std::string& name) {
          return C(values, name);
        },
        "grid_provider"_a,
        "values"_a,
        "name"_a = C::static_id());
  m.def(std::string(make_name).c_str(),
        [](const Grid::GridProvider<G, Grid::DD::SubdomainGrid<G>>& /*grid*/,
           const IntervalVector& value,
           const std::string& name) { return C(value, name); },
        "grid_provider"_a,
        "values"_a,
        "name"_a = C::static_id());
  return c;
} // ... bind_IndicatorGridFunction(...)


template <size_t d, size_t r, size_t rC>
typename pybind11::class_<IndicatorFunction<d, r, rC, double>> bind_IndicatorFunction(pybind11::module& m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  typedef double D;
  typedef double R;
  typedef FunctionInterface<d, r, rC, R> I;
  typedef IndicatorFunction<d, r, rC, R> C;
  using DomainType = typename C::DomainType;
  using RangeReturnType = typename C::RangeReturnType;
  using CornerVector = std::vector<std::tuple<DomainType, DomainType, RangeReturnType>>;
  using IntervalVector = std::vector<std::pair<Common::FieldMatrix<D, d, 2>, RangeReturnType>>;

  const std::string classname = std::string(
      "IndicatorFunction__" + Common::to_string(d) + "d_to_" + Common::to_string(r) + "x" + Common::to_string(rC));
  py::class_<C, I> c(m, classname.c_str(), classname.c_str());
  c.def(py::init<CornerVector, std::string>(), "values"_a, "name"_a = C::static_id());
  c.def(py::init<IntervalVector, std::string>(), "values"_a, "name"_a = C::static_id());

  return c;
} // ... bind_IndicatorFunction(...)


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // PYTHON_DUNE_XT_FUNCTIONS_INDICATOR_HH
