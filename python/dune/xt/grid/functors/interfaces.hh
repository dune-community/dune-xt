// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)
//   René Fritze     (2020)

#ifndef PYTHON_DUNE_XT_GRID_FUNCTORS_INTERFACES_HH
#define PYTHON_DUNE_XT_GRID_FUNCTORS_INTERFACES_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <python/dune/xt/common/timedlogging.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune::XT::Grid::bindings {


template <class G>
class ElementFunctor
{
  static_assert(is_grid<G>::value);
  using GV = typename G::LeafGridView;

public:
  using type = Grid::ElementFunctor<GV>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "element_functor",
                         const std::string& grid_id = grid_name<G>::value())
  {
    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def_readwrite("logger", &type::logger);
    c.def("prepare", [](type& self) { self.prepare(); });
    c.def("finalize", [](type& self) { self.finalize(); });
    return c;
  }
}; // class ElementFunctor


template <class G>
class IntersectionFunctor
{
  static_assert(is_grid<G>::value);
  using GV = typename G::LeafGridView;

public:
  using type = Grid::IntersectionFunctor<GV>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "intersection_functor",
                         const std::string& grid_id = grid_name<G>::value())
  {
    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def_readwrite("logger", &type::logger);
    c.def("prepare", [](type& self) { self.prepare(); });
    c.def("finalize", [](type& self) { self.finalize(); });
    return c;
  }
}; // class IntersectionFunctor


template <class G>
class ElementAndIntersectionFunctor
{
  static_assert(is_grid<G>::value);
  using GV = typename G::LeafGridView;

public:
  using type = Grid::ElementAndIntersectionFunctor<GV>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "element_and_intersection_functor",
                         const std::string& grid_id = grid_name<G>::value())
  {
    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def_readwrite("logger", &type::logger);
    c.def("prepare", [](type& self) { self.prepare(); });
    c.def("finalize", [](type& self) { self.finalize(); });
    return c;
  }
}; // class ElementAndIntersectionFunctor


} // namespace Dune::XT::Grid::bindings

#endif // PYTHON_DUNE_XT_GRID_FUNCTORS_INTERFACES_HH
