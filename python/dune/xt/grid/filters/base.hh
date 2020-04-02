// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_XT_GRID_FILTERS_BASE_HH
#define PYTHON_DUNE_XT_GRID_FILTERS_BASE_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/filters/base.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace bindings {


template <class G>
class ElementFilter
{
  static_assert(is_grid<G>::value, "");
  using GV = typename G::LeafGridView;
public:

  using type = Grid::ElementFilter<GV>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "element_filter",
                         const std::string& grid_id = grid_name<G>::value())
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());

    c.def("not", [](type& self){return !self;});
    c.def("and", [](type& self, const type& other){return self && other;}, "other"_a);
    c.def("or", [](type& self, const type& other){return self || other;}, "other_a");

    return c;
  } // ... bind(...)
}; // class ElementFilter


template <class G>
class IntersectionFilter
{
  static_assert(is_grid<G>::value, "");
  using GV = typename G::LeafGridView;
public:

  using type = Grid::IntersectionFilter<GV>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "intersection_filter",
                         const std::string& grid_id = grid_name<G>::value())
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());

    c.def("not", [](type& self){return !self;});
    c.def("and", [](type& self, const type& other){return self && other;}, "other"_a);
    c.def("or", [](type& self, const type& other){return self || other;}, "other_a");

    return c;
  } // ... bind(...)
}; // class IntersectionFilter


} // namespace bindings
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // PYTHON_DUNE_XT_GRID_FILTERS_BASE_HH
