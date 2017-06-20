// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include <cmath>

#include <dune/xt/common/test/main.hxx>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/functions/lambda/global-function.hh>
#include <dune/xt/functions/lambda/global-flux-function.hh>

#include "functions.hh"

using namespace Dune;
using namespace Dune::XT;


struct GlobalLambdaFluxFunctionTest : public ::testing::Test
{
  void check() const
  {
    typedef TESTFUNCTIONTYPE U;
    U u([](typename U::DomainType xx) { return xx[0]; }, 1);

    typedef XT::Functions::GlobalLambdaFluxFunction<U> FluxType;
    FluxType F([](const typename FluxType::DomainType& /*xx*/,
                  const typename FluxType::StateRangeType& uu,
                  const XT::Common::Parameter& mu) { return std::pow(uu[0], mu.get("power").at(0)); },
               XT::Common::ParameterType("power", 1),
               "burgers_flux",
               [](const XT::Common::Parameter& mu) { return mu.get("power").at(0); });

    auto grid = XT::Grid::make_cube_grid<GRIDTYPE>();

    for (auto&& entity : elements(grid.leaf_view())) {
      auto xx_global = entity.geometry().center();
      auto xx_local = entity.geometry().local(xx_global);
      auto u_value = u.local_function(entity)->evaluate(xx_local)[0];
      ASSERT_EQ(std::pow(u_value, 2.), F.local_function(entity)->evaluate(xx_local, u_value, {"power", 2.}));
      ASSERT_EQ(std::pow(u_value, 2.), F.evaluate(xx_global, u_value, {"power", 2.}));
      ASSERT_EQ(2, F.local_function(entity)->order({"power", 2.}));
      ASSERT_EQ(2, F.order({"power", 2.}));
    }
  }
};


TEST_F(GlobalLambdaFluxFunctionTest, creation_and_evalution)
{
  this->check();
}
