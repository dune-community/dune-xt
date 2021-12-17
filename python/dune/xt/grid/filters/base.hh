// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)
//   René Fritze     (2020)
//   Tim Keil        (2021)
//   Tobias Leibner  (2021)

#ifndef PYTHON_DUNE_XT_GRID_FILTERS_BASE_HH
#define PYTHON_DUNE_XT_GRID_FILTERS_BASE_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/filters/base.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune::XT::Grid::bindings {


template <class GV>
class ElementFilter
{
  using G = typename GV::Grid;
  static_assert(is_grid<G>::value);

public:
  using type = Grid::ElementFilter<GV>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id = grid_name<G>::value(),
                         const std::string& layer_id = "",
                         const std::string& class_id = "element_filter")
  {
    using namespace pybind11::literals;

    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
    if (!layer_id.empty())
      ClassName += "_" + layer_id;
    bound_type c(m, ClassName.c_str(), ClassName.c_str());

    c.def("not", [](type& self) { return !self; });
    c.def(
        "and", [](type& self, const type& other) { return self && other; }, "other"_a);
    c.def(
        "or", [](type& self, const type& other) { return self || other; }, "other_a");

    return c;
  } // ... bind(...)
}; // class ElementFilter


template <class GV>
class IntersectionFilter
{
  using G = typename GV::Grid;
  static_assert(is_grid<G>::value);

public:
  using type = Grid::IntersectionFilter<GV>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id = grid_name<G>::value(),
                         const std::string& layer_id = "",
                         const std::string& class_id = "intersection_filter")
  {
    using namespace pybind11::literals;

    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
    if (!layer_id.empty())
      ClassName += "_" + layer_id;
    bound_type c(m, ClassName.c_str(), ClassName.c_str());

    c.def("not", [](type& self) { return !self; });
    c.def(
        "and", [](type& self, const type& other) { return self && other; }, "other"_a);
    c.def(
        "or", [](type& self, const type& other) { return self || other; }, "other_a");

    return c;
  } // ... bind(...)
}; // class IntersectionFilter


} // namespace Dune::XT::Grid::bindings

#endif // PYTHON_DUNE_XT_GRID_FILTERS_BASE_HH
