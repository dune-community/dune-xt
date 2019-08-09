// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#include <dune/xt/common/test/main.hxx> // <- Has to come first, includes the config.h!

#include <dune/xt/common/unused.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/functions/generic/grid-function.hh>
#include <dune/xt/functions/grid-function.hh>


GTEST_TEST(ScalarGridFunction, constructible_from_grid_function_rvalueref)
{
  using namespace Dune;
  using namespace Dune::XT::Functions;

  using G = CUBEGRID_2D;
  using E = XT::Grid::extract_entity_t<G>;

  const GenericGridFunction<E> scalar_grid_function{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"};

  GridFunction<E> func{scalar_grid_function};
  EXPECT_EQ(std::string("THE_NAME"), func.name());
} // GTEST_TEST(ScalarGridFunction, constructible_from_grid_function_rvalueref)


GTEST_TEST(MatrixGridFunction, constructible_from_scalar_grid_function_rvalueref)
{
  using namespace Dune;
  using namespace Dune::XT::Functions;

  using G = CUBEGRID_2D;
  using E = XT::Grid::extract_entity_t<G>;

  const GenericGridFunction<E> scalar_grid_function{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"};

  GridFunction<E, 2, 2> DXTC_UNUSED(func){scalar_grid_function};
} // GTEST_TEST(MatrixGridFunction, constructible_from_scalar_grid_function_rvalueref)
