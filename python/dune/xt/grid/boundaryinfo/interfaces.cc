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

#include "interfaces.hh"


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct BoundaryInfo_for_all_grids
{
  static void bind(pybind11::module& m)
  {
    Dune::XT::Grid::bindings::BoundaryInfo<typename GridTypes::head_type>::bind(m);
    BoundaryInfo_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct BoundaryInfo_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_grid_boundaryinfo_interfaces, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common.timedlogging");

  BoundaryInfo_for_all_grids<>::bind(m);
}
