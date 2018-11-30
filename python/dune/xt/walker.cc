// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2018)

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
#include <python/dune/xt/grid/boundaryinfo.bindings.hh>
#include <python/dune/xt/grid/walker.bindings.hh>
#include <python/dune/xt/grid/walker/apply-on.bindings.hh>

#include <python/dune/xt/grid/available_types.hh>

template <class G, Dune::XT::Grid::Layers layer, Dune::XT::Grid::Backends backend>
void bind_walker(pybind11::module& m)
{
  try {
    Dune::XT::Grid::bindings::Walker<G, layer, backend>::bind(m);
  } catch (std::runtime_error&) {
  }
} // ... bind_walker(...)


template <class Tuple = Dune::XT::Grid::bindings::AvailableTypes>
void addbind_for_Grid(pybind11::module& m)
{
  using namespace Dune::XT::Grid;
  using G = typename Tuple::head_type;

  bind_walker<G, Layers::dd_subdomain, Backends::view>(m);
  bind_walker<G, Layers::dd_subdomain_boundary, Backends::view>(m);
  bind_walker<G, Layers::dd_subdomain_coupling, Backends::view>(m);
  bind_walker<G, Layers::dd_subdomain_oversampled, Backends::view>(m);
  bind_walker<G, Layers::leaf, Backends::view>(m);
  bind_walker<G, Layers::level, Backends::view>(m);

  addbind_for_Grid<typename Tuple::tail_type>(m);
} // ... addbind_for_Grid(...)

template <>
void addbind_for_Grid<boost::tuples::null_type>(pybind11::module&)
{
}


PYBIND11_MODULE(_walker, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;
  addbind_for_Grid(m);
  DUNE_XT_GRID_WALKER_APPLYON_BIND(m);
}
