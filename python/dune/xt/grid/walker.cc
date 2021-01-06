// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017, 2019 - 2020)
//   René Fritze     (2018 - 2020)
//   Tim Keil        (2018)
//   Tobias Leibner  (2018 - 2020)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/grid/view/coupling.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

#include "walker.hh"


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct Walker_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;

  static void bind(pybind11::module& m)
  {
    using Dune::XT::Grid::bindings::grid_name;
    Dune::XT::Grid::bindings::Walker<GV>::bind(m, grid_name<G>::value(), "leaf");
//    Dune::XT::Grid::bindings::Walker<GV>::bind_leaf_factory(m);
    Walker_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct Walker_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};

template <class GridTypes = Dune::XT::Grid::bindings::Available2dGridTypes>
struct Walker_for_all_coupling_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GridGlueType = Dune::XT::Grid::DD::Glued<G,G,Dune::XT::Grid::Layers::leaf>;
  using CGV = Dune::XT::Grid::CouplingGridView<GridGlueType>;

  static void bind(pybind11::module& m)
  {
    using Dune::XT::Grid::bindings::grid_name;
    Dune::XT::Grid::bindings::Walker<CGV>::bind(m, grid_name<G>::value(), "coupling");
//    Dune::XT::Grid::bindings::Walker<CGV>::bind_coupling_factory(m);
    Walker_for_all_coupling_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct Walker_for_all_coupling_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_grid_walker, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.grid._grid_element");
  py::module::import("dune.xt.grid._grid_filters_base");
  py::module::import("dune.xt.grid._grid_functors_interfaces");
  py::module::import("dune.xt.grid._grid_gridprovider_provider");
  py::module::import("dune.xt.grid._grid_intersection");

  Walker_for_all_grids<>::bind(m);
//  Walker_for_all_coupling_grids<>::bind(m);
}
