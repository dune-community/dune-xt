// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2018)
//   Tobias Leibner  (2016)

/**
 * This file is intended as a starting point for quick testing.
 */

#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/xt/common/test/main.hxx>

#include <dune/common/unused.hh>

#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/boundaryinfo.hh>


GTEST_TEST(BoundaryInfoFactory, all)
{
  typedef Dune::XT::Grid::BoundaryInfoFactory<INTERSECTIONTYPE> Factory;

  for (auto type : Factory::available()) {
    auto opts = Factory::default_config(type);
    auto aa DUNE_UNUSED = Factory::create();
    auto bb DUNE_UNUSED = Factory::create(type);
    auto cc DUNE_UNUSED = Factory::create(opts);
  }
}
