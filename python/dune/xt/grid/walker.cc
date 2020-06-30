// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017, 2019)
//   René Fritze     (2018 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2018 - 2020)

#include "config.h"

#include <dune/pybindxi/pybind11.h>

#include "walker.hh"


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct Walker_for_all_grids
{
  static void bind(pybind11::module& m)
  {
    Dune::XT::Grid::bindings::Walker<typename GridTypes::head_type>::bind(m);
    Walker_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct Walker_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_grid_walker, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.grid._grid_gridprovider_provider");
  py::module::import("dune.xt.grid._grid_filters_base");
  py::module::import("dune.xt.grid._grid_functors_interfaces");

  Walker_for_all_grids<>::bind(m);
}
