// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019 - 2020)
//   René Fritze     (2019 - 2020)
//   Tobias Leibner  (2019 - 2020)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/python.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

#include "indicator.hh"


template <class G>
void addbind_for_Grid(pybind11::module& m)
{
  using namespace Dune::XT::Functions;
  const auto grid_id = Dune::XT::Grid::bindings::grid_name<G>::value();
  const auto g_dim = G::dimension;

  bind_IndicatorGridFunction<G, g_dim, 1, 1>(m, grid_id);
  bind_IndicatorGridFunction<G, g_dim, 2, 1>(m, grid_id);
  bind_IndicatorGridFunction<G, g_dim, 2, 2>(m, grid_id);
  bind_IndicatorGridFunction<G, g_dim, 3, 1>(m, grid_id);
  bind_IndicatorGridFunction<G, g_dim, 3, 3>(m, grid_id);
} // ... addbind_for_Grid(...)


template <class Tuple = Dune::XT::Grid::bindings::AvailableGridTypes>
void all_grids(pybind11::module& m)
{
  Dune::XT::Common::bindings::guarded_bind([&]() { //  different grids but same entity
    addbind_for_Grid<Dune::XT::Common::tuple_head_t<Tuple>>(m);
  });
  all_grids<Dune::XT::Common::tuple_tail_t<Tuple>>(m);
} // ... addbind_for_Grid(...)


template <>
void all_grids<Dune::XT::Common::tuple_null_type>(pybind11::module&)
{}


PYBIND11_MODULE(_functions_indicator, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_1d");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_2d");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_3d");

  all_grids(m);

  using namespace Dune::XT::Functions;

  bind_IndicatorFunction<1, 1, 1>(m);
  bind_IndicatorFunction<1, 2, 1>(m);
  bind_IndicatorFunction<1, 2, 2>(m);
  bind_IndicatorFunction<1, 3, 1>(m);
  bind_IndicatorFunction<1, 3, 3>(m);

  bind_IndicatorFunction<2, 1, 1>(m);
  bind_IndicatorFunction<2, 2, 1>(m);
  bind_IndicatorFunction<2, 2, 2>(m);
  bind_IndicatorFunction<2, 3, 1>(m);
  bind_IndicatorFunction<2, 3, 3>(m);

  bind_IndicatorFunction<3, 1, 1>(m);
  bind_IndicatorFunction<3, 2, 1>(m);
  bind_IndicatorFunction<3, 2, 2>(m);
  bind_IndicatorFunction<3, 3, 1>(m);
  bind_IndicatorFunction<3, 3, 3>(m);
}
