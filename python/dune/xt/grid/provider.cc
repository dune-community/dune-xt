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

#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/layers.hh>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/gridprovider.hh>


template <class Tuple = Dune::XT::Grid::AvailableGridTypes>
void addbind_for_Grid(pybind11::module& m)
{
  using namespace Dune::XT::Grid;
  using G = typename Tuple::head_type;
  const auto grid_id = Dune::XT::Grid::bindings::grid_name<G>::value();

  // bind_GridProvider<G>(m, grid_id);
  bind_make_cube_grid<G>(m, grid_id);

  addbind_for_Grid<typename Tuple::tail_type>(m);
} // ... addbind_for_Grid(...)

template <>
void addbind_for_Grid<boost::tuples::null_type>(pybind11::module&)
{}


PYBIND11_MODULE(_provider, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");

  addbind_for_Grid(m);
}
