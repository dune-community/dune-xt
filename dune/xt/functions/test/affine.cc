// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2013 - 2018)
//   Tobias Leibner  (2014 - 2017)

#include <dune/xt/common/test/main.hxx>

#include <memory>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/functions/affine.hh>

#include <dune/xt/functions/test/functions.hh>

using namespace Dune;
using namespace Dune::XT;
using namespace Dune::XT::Functions;

struct FunctionsTest : public FunctionTest<TESTFUNCTIONTYPE>
{
  void check()
  {
    auto grid = XT::Grid::make_cube_grid<GRIDTYPE>();
    const auto testfunction = TESTFUNCTIONTYPE::create();
    for (auto&& entity : elements(grid.leaf_view())) {
      auto xx_global = entity.geometry().center();
      auto xx_local = entity.geometry().local(xx_global);
      TESTFUNCTIONTYPE::JacobianRangeType jacobian;
      testfunction->local_function(entity)->jacobian(xx_local);
    }
  }
};
