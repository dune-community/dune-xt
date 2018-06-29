// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2018)

#include "config.h"

#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/exceptions.bindings.hh>

#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/layers.hh>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/gridprovider.hh>

#include <python/dune/xt/grid/available_types.hh>

template <class Tuple = Dune::XT::Grid::bindings::AvailableTypes>
void addbind_for_Grid(pybind11::module& m)
{
  using namespace Dune::XT::Grid;
  using G = typename Tuple::head_type;
  const auto grid_id = Dune::XT::Grid::bindings::grid_name<G>::value();
  typedef typename Layer<G, Layers::dd_subdomain, Backends::view, DD::SubdomainGrid<G>>::type DdSubdomainPart;

  bind_GridProvider<G>(m, grid_id);
  bind_make_cube_grid<G>(m, grid_id);

  bind_DdSubdomainsGridProvider<G>(m, grid_id);
  bind_make_cube_dd_subdomains_grid<G>(m, grid_id);

  addbind_for_Grid<typename Tuple::tail_type>(m);
} // ... addbind_for_Grid(...)

template <>
void addbind_for_Grid<boost::tuples::null_type>(pybind11::module& m)
{
}

PYBIND11_MODULE(_provider, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  addbind_for_Grid(m);
}
