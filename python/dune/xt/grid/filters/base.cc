// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)
//   René Fritze     (2020)

#include "config.h"

#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/grid/view/coupling.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

#include "base.hh"


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct ElementFilter_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using LGV = typename G::LeafGridView;
  static const size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::XT::Grid::bindings::grid_name;
    Dune::XT::Grid::bindings::ElementFilter<LGV>::bind(m, grid_name<G>::value(), "leaf");
    if constexpr (d == 2) {
      using GridGlueType = Dune::XT::Grid::DD::Glued<G, G, Dune::XT::Grid::Layers::leaf>;
      using CGV = Dune::XT::Grid::CouplingGridView<GridGlueType>;
      Dune::XT::Grid::bindings::ElementFilter<CGV>::bind(m, grid_name<G>::value(), "coupling");
    }
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
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using LGV = typename G::LeafGridView;
  static const size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::XT::Grid::bindings::grid_name;
    Dune::XT::Grid::bindings::IntersectionFilter<LGV>::bind(m, grid_name<G>::value(), "leaf");
    if constexpr (d == 2) {
      using GridGlueType = Dune::XT::Grid::DD::Glued<G, G, Dune::XT::Grid::Layers::leaf>;
      using CGV = Dune::XT::Grid::CouplingGridView<GridGlueType>;
      Dune::XT::Grid::bindings::IntersectionFilter<CGV>::bind(m, grid_name<G>::value(), "coupling");
    }
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
