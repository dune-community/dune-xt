// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <python/dune/xt/grid/grids.bindings.hh>

#include "base.hh"


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct ElementFilter_for_all_grids
{
  static void bind(pybind11::module& m)
  {
    Dune::XT::Grid::bindings::ElementFilter<Dune::XT::Common::tuple_head_t<GridTypes>>::bind(m);
    ElementFilter_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct ElementFilter_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct IntersectionFilter_for_all_grids
{
  static void bind(pybind11::module& m)
  {
    Dune::XT::Grid::bindings::IntersectionFilter<Dune::XT::Common::tuple_head_t<GridTypes>>::bind(m);
    IntersectionFilter_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct IntersectionFilter_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_grid_filters_base, m)
{
  ElementFilter_for_all_grids<>::bind(m);
  IntersectionFilter_for_all_grids<>::bind(m);
}
