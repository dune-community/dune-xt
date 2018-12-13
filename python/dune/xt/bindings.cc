// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2018)
//   René Fritze     (2018)
//   Tim Keil        (2018)

#include "config.h"

#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

#include <python/dune/xt/functions/constant.hh>
#include <python/dune/xt/functions/checkerboard.hh>
#include <python/dune/xt/functions/ESV2007.bindings.hh>
#include <python/dune/xt/functions/expression.hh>
#include <python/dune/xt/functions/spe10.hh>
#include <python/dune/xt/functions/indicator.hh>

#include <python/dune/xt/common/exceptions.bindings.hh>

template <class G>
void addbind_for_Grid(pybind11::module& m)
{
  using namespace Dune::XT::Functions;
  const auto grid_id = Dune::XT::Grid::bindings::grid_name<G>::value();
  const auto g_dim = G::dimension;

  bind_CheckerboardFunction<G, g_dim, 1, 1>(m, grid_id);
  bind_CheckerboardFunction<G, g_dim, 2, 1>(m, grid_id);
  bind_CheckerboardFunction<G, g_dim, 3, 1>(m, grid_id);
  bind_CheckerboardFunction<G, g_dim, 4, 1>(m, grid_id);
  bind_CheckerboardFunction<G, g_dim, 1, 3>(m, grid_id);
  bind_CheckerboardFunction<G, g_dim, 3, 3>(m, grid_id);

  //  bindings::ESV2007::CutoffFunction<G>::bind(m);

  bind_Spe10Model1Function<G, g_dim, 1, 1>(m, grid_id);

  bind_IndicatorGridFunction<G, g_dim, 1, 1>(m, grid_id);
  bind_IndicatorGridFunction<G, g_dim, 2, 1>(m, grid_id);
  bind_IndicatorGridFunction<G, g_dim, 3, 1>(m, grid_id);
  bind_IndicatorGridFunction<G, g_dim, 4, 1>(m, grid_id);
  bind_IndicatorGridFunction<G, g_dim, 1, 3>(m, grid_id);
  bind_IndicatorGridFunction<G, g_dim, 3, 3>(m, grid_id);
} // ... addbind_for_Grid(...)


PYBIND11_MODULE(_functions, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;
  using namespace Dune::XT::Functions;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt._interfaces");

  bind_ConstantFunction<1, 1, 1>(m);
  bind_ConstantFunction<1, 2, 1>(m);
  bind_ConstantFunction<1, 3, 1>(m);
  bind_ConstantFunction<1, 4, 1>(m);
  bind_ConstantFunction<1, 1, 3>(m);
  bind_ConstantFunction<1, 3, 3>(m);

  bind_ConstantFunction<2, 1, 1>(m);
  bind_ConstantFunction<2, 2, 1>(m);
  bind_ConstantFunction<2, 3, 1>(m);
  bind_ConstantFunction<2, 4, 1>(m);
  bind_ConstantFunction<2, 1, 3>(m);
  bind_ConstantFunction<2, 3, 3>(m);

  bind_ConstantFunction<3, 1, 1>(m);
  bind_ConstantFunction<3, 2, 1>(m);
  bind_ConstantFunction<3, 3, 1>(m);
  bind_ConstantFunction<3, 4, 1>(m);
  bind_ConstantFunction<3, 1, 3>(m);
  bind_ConstantFunction<3, 3, 3>(m);

  bind_IndicatorFunction<1, 1, 1>(m);
  bind_IndicatorFunction<1, 2, 1>(m);
  bind_IndicatorFunction<1, 3, 1>(m);
  bind_IndicatorFunction<1, 4, 1>(m);
  bind_IndicatorFunction<1, 1, 3>(m);
  bind_IndicatorFunction<1, 3, 3>(m);

  bind_IndicatorFunction<2, 1, 1>(m);
  bind_IndicatorFunction<2, 2, 1>(m);
  bind_IndicatorFunction<2, 3, 1>(m);
  bind_IndicatorFunction<2, 4, 1>(m);
  bind_IndicatorFunction<2, 1, 3>(m);
  bind_IndicatorFunction<2, 3, 3>(m);

  bind_IndicatorFunction<3, 1, 1>(m);
  bind_IndicatorFunction<3, 2, 1>(m);
  bind_IndicatorFunction<3, 3, 1>(m);
  bind_IndicatorFunction<3, 4, 1>(m);
  bind_IndicatorFunction<3, 1, 3>(m);
  bind_IndicatorFunction<3, 3, 3>(m);

  bind_ExpressionFunction<1, 1, 1>(m);
  bind_ExpressionFunction<1, 2, 1>(m);
  bind_ExpressionFunction<1, 3, 1>(m);
  bind_ExpressionFunction<1, 4, 1>(m);
  bind_ExpressionFunction<1, 1, 3>(m);
  bind_ExpressionFunction<1, 2, 2>(m);
  bind_ExpressionFunction<1, 3, 3>(m);

  bind_ExpressionFunction<2, 1, 1>(m);
  bind_ExpressionFunction<2, 2, 1>(m);
  bind_ExpressionFunction<2, 3, 1>(m);
  bind_ExpressionFunction<2, 4, 1>(m);
  bind_ExpressionFunction<2, 1, 3>(m);
  bind_ExpressionFunction<2, 2, 2>(m);
  bind_ExpressionFunction<2, 3, 3>(m);

  bind_ExpressionFunction<3, 1, 1>(m);
  bind_ExpressionFunction<3, 2, 1>(m);
  bind_ExpressionFunction<3, 3, 1>(m);
  bind_ExpressionFunction<3, 4, 1>(m);
  bind_ExpressionFunction<3, 1, 3>(m);
  bind_ExpressionFunction<3, 2, 2>(m);
  bind_ExpressionFunction<3, 3, 3>(m);

  addbind_for_Grid<Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>>(m);
  addbind_for_Grid<Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>>(m);
#if HAVE_DUNE_ALUGRID
  addbind_for_Grid<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>(m);
#endif
#if HAVE_UG
  addbind_for_Grid<Dune::UGGrid<2>>(m);
#endif
  //#if HAVE_ALBERTA
  //  addbind_for_Grid<Dune::AlbertaGrid<2, 2>>(m, "2d_simplex_albertagrid");
  //#endif

  Dune::XT::Common::bindings::add_initialization(m, "dune.xt.functions");
}
