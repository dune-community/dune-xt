// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)
//   René Fritze     (2020)
//   Tim Keil        (2021)
//   Tobias Leibner  (2021)

#include "config.h"

#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/view/coupling.hh>

#include <python/dune/xt/grid/filters/element.hh>
#include <python/dune/xt/grid/grids.bindings.hh>


template <template <class> class Filter, class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct InitlessElementFilter_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using LGV = typename G::LeafGridView;
  static const size_t d = G::dimension;

  static void bind(pybind11::module& m, const std::string& class_id)
  {
    using Dune::XT::Grid::bindings::grid_name;
    Dune::XT::Grid::bindings::InitlessElementFilter<Filter, LGV>::bind(m, class_id, "leaf");
    Dune::XT::Grid::bindings::InitlessElementFilter<Filter, LGV>::bind_leaf_factory(m, class_id);
#if HAVE_DUNE_GRID_GLUE
    if constexpr (d < 3) {
      using GridGlueType = Dune::XT::Grid::DD::Glued<G, G, Dune::XT::Grid::Layers::leaf>;
      using CGV = Dune::XT::Grid::CouplingGridView<GridGlueType>;
      Dune::XT::Grid::bindings::InitlessElementFilter<Filter, CGV>::bind(m, class_id, "coupling");
      Dune::XT::Grid::bindings::InitlessElementFilter<Filter, CGV>::bind_coupling_factory(m, class_id);
    }
#endif
    InitlessElementFilter_for_all_grids<Filter, Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m, class_id);
  }
};

template <template <class> class Filter>
struct InitlessElementFilter_for_all_grids<Filter, Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/, const std::string& /*class_id*/) {}
};


PYBIND11_MODULE(_grid_filters_element, m)
{
  namespace py = pybind11;
  using namespace Dune::XT::Grid;

  py::module::import("dune.xt.grid._grid_gridprovider_provider");
  py::module::import("dune.xt.grid._grid_filters_base");

#define BIND_(NAME) InitlessElementFilter_for_all_grids<ApplyOn::NAME>::bind(m, std::string("ApplyOn") + #NAME)

  BIND_(AllElements);
  BIND_(NoElements);
  BIND_(BoundaryElements);
  //  BIND_(GenericFilteredElements); <- not initless
  //  BIND_(PartitionSetElements); <- requires partitionset template parameter

#undef BIND_
}
