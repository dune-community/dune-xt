// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2018)
//   Tobias Leibner  (2018 - 2019)

#include "config.h"

#include <dune/xt/grid/grids.hh>

#include <python/dune/xt/grid/gridprovider.hh>


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct GridProvider_for_all_grids
{
  static void bind(pybind11::module& m)
  {
    Dune::XT::Grid::bindings::GridProvider<typename GridTypes::head_type>::bind(m);
    GridProvider_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct GridProvider_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_grid_gridprovider_provider, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");

  GridProvider_for_all_grids<>::bind(m);
}
