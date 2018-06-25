// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
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

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

#include <python/dune/xt/functions/interfaces.hh>
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
  const auto diff = internal::Combination::difference;
  const auto sum = internal::Combination::sum;
  const auto prod = internal::Combination::product;
  const auto g_dim = G::dimension;

  auto i_1_1 = bind_GridFunctionInterface<G, 1, 1>(m, grid_id);
  auto i_2_1 = bind_GridFunctionInterface<G, 2, 1>(m, grid_id);
  auto i_3_1 = bind_GridFunctionInterface<G, 3, 1>(m, grid_id);
  auto i_4_1 = bind_GridFunctionInterface<G, 4, 1>(m, grid_id);
  //! this generates multiple binds for the same type
  //! auto i_d_d = bind_GridFunctionInterface<G, g_dim, g_dim>(m, grid_id);
  auto i_d_d = bind_GridFunctionInterface<G, 2, 2>(m, grid_id);

  bind_combined_GridFunction<G, g_dim, diff, 1, 1, 1, 1>(m, grid_id);
  addbind_GridFunctionInterface_combined_op<G, g_dim, diff, 1, 1, 1, 1>(i_1_1);

  bind_combined_GridFunction<G, g_dim, diff, 2, 1, 2, 1>(m, grid_id);
  addbind_GridFunctionInterface_combined_op<G, g_dim, diff, 2, 1, 2, 1>(i_2_1);

  bind_combined_GridFunction<G, g_dim, diff, 3, 1, 3, 1>(m, grid_id);
  addbind_GridFunctionInterface_combined_op<G, g_dim, diff, 3, 1, 3, 1>(i_3_1);

  bind_combined_GridFunction<G, g_dim, diff, 4, 1, 4, 1>(m, grid_id);
  addbind_GridFunctionInterface_combined_op<G, g_dim, diff, 4, 1, 4, 1>(i_4_1);

  bind_combined_GridFunction<G, g_dim, diff, 2, 2, 2, 2>(m, grid_id);
  addbind_GridFunctionInterface_combined_op<G, g_dim, diff, g_dim, g_dim, g_dim, g_dim>(i_d_d);


  bind_combined_GridFunction<G, g_dim, sum, 1, 1, 1, 1>(m, grid_id);
  addbind_GridFunctionInterface_combined_op<G, g_dim, sum, 1, 1, 1, 1>(i_1_1);

  bind_combined_GridFunction<G, g_dim, sum, 2, 1, 2, 1>(m, grid_id);
  addbind_GridFunctionInterface_combined_op<G, g_dim, sum, 2, 1, 2, 1>(i_2_1);

  bind_combined_GridFunction<G, g_dim, sum, 3, 1, 3, 1>(m, grid_id);
  addbind_GridFunctionInterface_combined_op<G, g_dim, sum, 3, 1, 3, 1>(i_3_1);

  bind_combined_GridFunction<G, g_dim, sum, 4, 1, 4, 1>(m, grid_id);
  addbind_GridFunctionInterface_combined_op<G, g_dim, sum, 4, 1, 4, 1>(i_4_1);

  bind_combined_GridFunction<G, g_dim, sum, 2, 2, 2, 2>(m, grid_id);
  addbind_GridFunctionInterface_combined_op<G, g_dim, sum, g_dim, g_dim, g_dim, g_dim>(i_d_d);


  bind_combined_GridFunction<G, g_dim, prod, 1, 1, 1, 1>(m, grid_id);
  addbind_GridFunctionInterface_combined_op<G, g_dim, prod, 1, 1, 1, 1>(i_1_1);

  bind_combined_GridFunction<G, g_dim, prod, 1, 1, 2, 1>(m, grid_id);
  addbind_GridFunctionInterface_combined_op<G, g_dim, prod, 1, 1, 2, 1>(i_1_1);

  bind_combined_GridFunction<G, g_dim, prod, 1, 1, 3, 1>(m, grid_id);
  addbind_GridFunctionInterface_combined_op<G, g_dim, prod, 1, 1, 3, 1>(i_1_1);

  bind_combined_GridFunction<G, g_dim, prod, 1, 1, 4, 1>(m, grid_id);
  addbind_GridFunctionInterface_combined_op<G, g_dim, prod, 1, 1, 4, 1>(i_1_1);

  bind_combined_GridFunction<G, g_dim, prod, 1, 1, 2, 2>(m, grid_id);
  addbind_GridFunctionInterface_combined_op<G, g_dim, prod, 1, 1, 2, 2>(i_1_1);


  bind_ConstantFunction<G, g_dim, 1, 1>(m, grid_id);
  bind_ConstantFunction<G, g_dim, 2, 1>(m, grid_id);
  bind_ConstantFunction<G, g_dim, 3, 1>(m, grid_id);
  bind_ConstantFunction<G, g_dim, 4, 1>(m, grid_id);
  bind_ConstantFunction<G, g_dim, 2, 2>(m, grid_id);

  bind_CheckerboardFunction<G, g_dim, 1, 1>(m, grid_id);
  bind_CheckerboardFunction<G, g_dim, 2, 1>(m, grid_id);
  bind_CheckerboardFunction<G, g_dim, 3, 1>(m, grid_id);
  bind_CheckerboardFunction<G, g_dim, 4, 1>(m, grid_id);
  bind_CheckerboardFunction<G, g_dim, 1, 2>(m, grid_id);
  bind_CheckerboardFunction<G, g_dim, 1, 3>(m, grid_id);

  //  bindings::ESV2007::CutoffFunction<G>::bind(m);

  bind_ExpressionFunction<G, g_dim, 1, 1>(m, grid_id);
  bind_ExpressionFunction<G, g_dim, 2, 1>(m, grid_id);
  bind_ExpressionFunction<G, g_dim, 3, 1>(m, grid_id);
  bind_ExpressionFunction<G, g_dim, 4, 1>(m, grid_id);
  bind_ExpressionFunction<G, g_dim, 1, 2>(m, grid_id);
  bind_ExpressionFunction<G, g_dim, 1, 3>(m, grid_id);

  bind_Spe10Model1Function<G, g_dim, 1, 1>(m, grid_id);

  bind_IndicatorFunction<G, g_dim, 1, 1>(m, grid_id);
} // ... addbind_for_Grid(...)


PYBIND11_PLUGIN(_functions)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module m("_functions", "dune-xt-functions");

  Dune::XT::Common::bindings::addbind_exceptions(m);

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");

  addbind_for_Grid<Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>>(m);
  addbind_for_Grid<Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>>(m);
#if HAVE_DUNE_ALUGRID
  addbind_for_Grid<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>(m);
#endif
  //#if HAVE_UG
  //  addbind_for_Grid<Dune::UGGrid<2>>(m, "2d_simplex_uggrid");
  //#endif
  //#if HAVE_ALBERTA
  //  addbind_for_Grid<Dune::AlbertaGrid<2, 2>>(m, "2d_simplex_albertagrid");
  //#endif

  Dune::XT::Common::bindings::add_initialization(m, "dune.xt.functions");
  return m.ptr();
}
