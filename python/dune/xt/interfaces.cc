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

#include <python/dune/xt/common/exceptions.bindings.hh>

template <class G>
void add_bind_for_Grid_interface(pybind11::module& m)
{
  using namespace Dune::XT::Functions;
  const auto grid_id = Dune::XT::Grid::bindings::grid_name<G>::value();
  const auto diff = internal::Combinations::difference;
  const auto sum = internal::Combinations::sum;
  const auto prod = internal::Combinations::product;
  const auto g_dim = G::dimension;

  auto i_1_1 = bind_GridFunctionInterface<G, 1, 1>(m, grid_id);
  auto i_2_1 = bind_GridFunctionInterface<G, 2, 1>(m, grid_id);
  auto i_3_1 = bind_GridFunctionInterface<G, 3, 1>(m, grid_id);
  auto i_4_1 = bind_GridFunctionInterface<G, 4, 1>(m, grid_id);
  auto i_1_3 = bind_GridFunctionInterface<G, 1, 3>(m, grid_id);
  auto i_3_3 = bind_GridFunctionInterface<G, 3, 3>(m, grid_id);

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
} // ... addbind_for_Grid(...)


PYBIND11_MODULE(_interfaces, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;
  using namespace Dune::XT::Functions;

  const auto diff = internal::Combination::difference;
  const auto sum = internal::Combination::sum;
  const auto prod = internal::Combination::product;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");

  auto i_1_1_1 = bind_FunctionInterface<1, 1, 1>(m);
  auto i_1_2_1 = bind_FunctionInterface<1, 2, 1>(m);
  auto i_1_3_1 = bind_FunctionInterface<1, 3, 1>(m);
  auto i_1_4_1 = bind_FunctionInterface<1, 4, 1>(m);
  auto i_1_2_2 = bind_FunctionInterface<1, 2, 2>(m);
//  auto i_1_1_2 = bind_FunctionInterface<1, 1, 2>(m);
//  auto i_1_1_3 = bind_FunctionInterface<1, 1, 3>(m);
//  auto i_1_3_3 = bind_FunctionInterface<1, 3, 3>(m);

  auto i_2_1_1 = bind_FunctionInterface<2, 1, 1>(m);
  auto i_2_2_1 = bind_FunctionInterface<2, 2, 1>(m);
  auto i_2_3_1 = bind_FunctionInterface<2, 3, 1>(m);
  auto i_2_4_1 = bind_FunctionInterface<2, 4, 1>(m);
  auto i_2_2_2 = bind_FunctionInterface<2, 2, 2>(m);
//  auto i_2_1_2 = bind_FunctionInterface<2, 1, 2>(m);
//  auto i_2_1_3 = bind_FunctionInterface<2, 1, 3>(m);
//  auto i_2_3_3 = bind_FunctionInterface<2, 3, 3>(m);

  auto i_3_1_1 = bind_FunctionInterface<3, 1, 1>(m);
  auto i_3_2_1 = bind_FunctionInterface<3, 2, 1>(m);
  auto i_3_3_1 = bind_FunctionInterface<3, 3, 1>(m);
  auto i_3_4_1 = bind_FunctionInterface<3, 4, 1>(m);
  auto i_3_2_2 = bind_FunctionInterface<3, 2, 2>(m);
//  auto i_3_1_2 = bind_FunctionInterface<3, 1, 2>(m);
//  auto i_3_1_3 = bind_FunctionInterface<3, 1, 3>(m);
//  auto i_3_3_3 = bind_FunctionInterface<3, 3, 3>(m);

  bind_combined_Function<1, diff, 1, 1, 1, 1>(m);
  addbind_FunctionInterface_combined_op<1, diff, 1, 1, 1, 1>(i_1_1_1);

  bind_combined_Function<1, diff, 2, 1, 2, 1>(m);
  addbind_FunctionInterface_combined_op<1, diff, 2, 1, 2, 1>(i_1_2_1);

  bind_combined_Function<1, diff, 3, 1, 3, 1>(m);
  addbind_FunctionInterface_combined_op<1, diff, 3, 1, 3, 1>(i_1_3_1);

  bind_combined_Function<1, diff, 4, 1, 4, 1>(m);
  addbind_FunctionInterface_combined_op<1, diff, 4, 1, 4, 1>(i_1_4_1);

  bind_combined_Function<1, diff, 2, 2, 2, 2>(m);
  addbind_FunctionInterface_combined_op<1, diff, 2, 2, 2, 2>(i_1_2_2);


  bind_combined_Function<1, sum, 1, 1, 1, 1>(m);
  addbind_FunctionInterface_combined_op<1, sum, 1, 1, 1, 1>(i_1_1_1);

  bind_combined_Function<1, sum, 2, 1, 2, 1>(m);
  addbind_FunctionInterface_combined_op<1, sum, 2, 1, 2, 1>(i_1_2_1);

  bind_combined_Function<1, sum, 3, 1, 3, 1>(m);
  addbind_FunctionInterface_combined_op<1, sum, 3, 1, 3, 1>(i_1_3_1);

  bind_combined_Function<1, sum, 4, 1, 4, 1>(m);
  addbind_FunctionInterface_combined_op<1, sum, 4, 1, 4, 1>(i_1_4_1);

  bind_combined_Function<1, sum, 2, 2, 2, 2>(m);
  addbind_FunctionInterface_combined_op<1, sum, 2, 2, 2, 2>(i_1_2_2);


  bind_combined_Function<1, prod, 1, 1, 1, 1>(m);
  addbind_FunctionInterface_combined_op<1, prod, 1, 1, 1, 1>(i_1_1_1);

  bind_combined_Function<1, prod, 1, 1, 2, 1>(m);
  addbind_FunctionInterface_combined_op<1, prod, 1, 1, 2, 1>(i_1_1_1);

  bind_combined_Function<1, prod, 1, 1, 3, 1>(m);
  addbind_FunctionInterface_combined_op<1, prod, 1, 1, 3, 1>(i_1_1_1);

  bind_combined_Function<1, prod, 1, 1, 4, 1>(m);
  addbind_FunctionInterface_combined_op<1, prod, 1, 1, 4, 1>(i_1_1_1);

  bind_combined_Function<1, prod, 1, 1, 2, 2>(m);
  addbind_FunctionInterface_combined_op<1, prod, 1, 1, 2, 2>(i_1_1_1);


  bind_combined_Function<2, diff, 1, 1, 1, 1>(m);
  addbind_FunctionInterface_combined_op<2, diff, 1, 1, 1, 1>(i_2_1_1);

  bind_combined_Function<2, diff, 2, 1, 2, 1>(m);
  addbind_FunctionInterface_combined_op<2, diff, 2, 1, 2, 1>(i_2_2_1);

  bind_combined_Function<2, diff, 3, 1, 3, 1>(m);
  addbind_FunctionInterface_combined_op<2, diff, 3, 1, 3, 1>(i_2_3_1);

  bind_combined_Function<2, diff, 4, 1, 4, 1>(m);
  addbind_FunctionInterface_combined_op<2, diff, 4, 1, 4, 1>(i_2_4_1);

  bind_combined_Function<2, diff, 2, 2, 2, 2>(m);
  addbind_FunctionInterface_combined_op<2, diff, 2, 2, 2, 2>(i_2_2_2);


  bind_combined_Function<2, sum, 1, 1, 1, 1>(m);
  addbind_FunctionInterface_combined_op<2, sum, 1, 1, 1, 1>(i_2_1_1);

  bind_combined_Function<2, sum, 2, 1, 2, 1>(m);
  addbind_FunctionInterface_combined_op<2, sum, 2, 1, 2, 1>(i_2_2_1);

  bind_combined_Function<2, sum, 3, 1, 3, 1>(m);
  addbind_FunctionInterface_combined_op<2, sum, 3, 1, 3, 1>(i_2_3_1);

  bind_combined_Function<2, sum, 4, 1, 4, 1>(m);
  addbind_FunctionInterface_combined_op<2, sum, 4, 1, 4, 1>(i_2_4_1);

  bind_combined_Function<2, sum, 2, 2, 2, 2>(m);
  addbind_FunctionInterface_combined_op<2, sum, 2, 2, 2, 2>(i_2_2_2);


  bind_combined_Function<2, prod, 1, 1, 1, 1>(m);
  addbind_FunctionInterface_combined_op<2, prod, 1, 1, 1, 1>(i_2_1_1);

  bind_combined_Function<2, prod, 1, 1, 2, 1>(m);
  addbind_FunctionInterface_combined_op<2, prod, 1, 1, 2, 1>(i_2_1_1);

  bind_combined_Function<2, prod, 1, 1, 3, 1>(m);
  addbind_FunctionInterface_combined_op<2, prod, 1, 1, 3, 1>(i_2_1_1);

  bind_combined_Function<2, prod, 1, 1, 4, 1>(m);
  addbind_FunctionInterface_combined_op<2, prod, 1, 1, 4, 1>(i_2_1_1);

  bind_combined_Function<2, prod, 1, 1, 2, 2>(m);
  addbind_FunctionInterface_combined_op<2, prod, 1, 1, 2, 2>(i_2_1_1);


  bind_combined_Function<3, diff, 1, 1, 1, 1>(m);
  addbind_FunctionInterface_combined_op<3, diff, 1, 1, 1, 1>(i_3_1_1);

  bind_combined_Function<3, diff, 2, 1, 2, 1>(m);
  addbind_FunctionInterface_combined_op<3, diff, 2, 1, 2, 1>(i_3_2_1);

  bind_combined_Function<3, diff, 3, 1, 3, 1>(m);
  addbind_FunctionInterface_combined_op<3, diff, 3, 1, 3, 1>(i_3_3_1);

  bind_combined_Function<3, diff, 4, 1, 4, 1>(m);
  addbind_FunctionInterface_combined_op<3, diff, 4, 1, 4, 1>(i_3_4_1);

  bind_combined_Function<3, diff, 2, 2, 2, 2>(m);
  addbind_FunctionInterface_combined_op<3, diff, 2, 2, 2, 2>(i_3_2_2);


  bind_combined_Function<3, sum, 1, 1, 1, 1>(m);
  addbind_FunctionInterface_combined_op<3, sum, 1, 1, 1, 1>(i_3_1_1);

  bind_combined_Function<3, sum, 2, 1, 2, 1>(m);
  addbind_FunctionInterface_combined_op<3, sum, 2, 1, 2, 1>(i_3_2_1);

  bind_combined_Function<3, sum, 3, 1, 3, 1>(m);
  addbind_FunctionInterface_combined_op<3, sum, 3, 1, 3, 1>(i_3_3_1);

  bind_combined_Function<3, sum, 4, 1, 4, 1>(m);
  addbind_FunctionInterface_combined_op<3, sum, 4, 1, 4, 1>(i_3_4_1);

  bind_combined_Function<3, sum, 2, 2, 2, 2>(m);
  addbind_FunctionInterface_combined_op<3, sum, 2, 2, 2, 2>(i_3_2_2);


  bind_combined_Function<3, prod, 1, 1, 1, 1>(m);
  addbind_FunctionInterface_combined_op<3, prod, 1, 1, 1, 1>(i_3_1_1);

  bind_combined_Function<3, prod, 1, 1, 2, 1>(m);
  addbind_FunctionInterface_combined_op<3, prod, 1, 1, 2, 1>(i_3_1_1);

  bind_combined_Function<3, prod, 1, 1, 3, 1>(m);
  addbind_FunctionInterface_combined_op<3, prod, 1, 1, 3, 1>(i_3_1_1);

  bind_combined_Function<3, prod, 1, 1, 4, 1>(m);
  addbind_FunctionInterface_combined_op<3, prod, 1, 1, 4, 1>(i_3_1_1);

  bind_combined_Function<3, prod, 1, 1, 2, 2>(m);
  addbind_FunctionInterface_combined_op<3, prod, 1, 1, 2, 2>(i_3_1_1);

  add_bind_for_Grid_interface<Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>>(m);
  add_bind_for_Grid_interface<Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>>(m);
#if HAVE_DUNE_ALUGRID
  add_bind_for_Grid_interface<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>(m);
#endif
#if HAVE_UG
  add_bind_for_Grid_interface<Dune::UGGrid<2>>(m);
#endif
  //#if HAVE_ALBERTA
  //  addbind_for_Grid<Dune::AlbertaGrid<2, 2>>(m, "2d_simplex_albertagrid");
  //#endif

  Dune::XT::Common::bindings::add_initialization(m, "dune.xt.functions");
}
