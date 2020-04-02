// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <dune/xt/grid/grids.hh>

#include "base.hh"


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct ElementFilter_for_all_grids
{
  static void bind(pybind11::module& m)
  {
    Dune::XT::Grid::bindings::ElementFilter<typename GridTypes::head_type>::bind(m);
    ElementFilter_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct ElementFilter_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct IntersectionFilter_for_all_grids
{
  static void bind(pybind11::module& m)
  {
    Dune::XT::Grid::bindings::IntersectionFilter<typename GridTypes::head_type>::bind(m);
    IntersectionFilter_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct IntersectionFilter_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_grid_filters_base, m)
{
  ElementFilter_for_all_grids<>::bind(m);
  IntersectionFilter_for_all_grids<>::bind(m);
}
