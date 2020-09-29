// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2016, 2020)

/**
 * This file is intended as a starting point for quick testing.
 */

#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/xt/test/main.hxx>

#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/boundaryinfo.hh>


GTEST_TEST(BoundaryInfoFactory, all)
{
  using Factory = Dune::XT::Grid::BoundaryInfoFactory<INTERSECTIONTYPE>;

  for (auto type : Factory::available()) {
    auto opts = Factory::default_config(type);
    [[maybe_unused]] auto aa = Factory::create();
    [[maybe_unused]] auto bb = Factory::create(type);
    [[maybe_unused]] auto cc = Factory::create(opts);
  }
}
