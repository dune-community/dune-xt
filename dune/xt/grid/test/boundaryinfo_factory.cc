// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:

/**
  * This file is intended as a starting point for quick testing.
  */

#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/xt/common/test/main.hxx>

#include <dune/common/unused.hh>

#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/boundaryinfo.hh>


TEST(BoundaryInfoFactory, all)
{
  typedef Dune::XT::Grid::BoundaryInfoFactory<INTERSECTIONTYPE> Factory;

  for (auto type : Factory::available()) {
    auto opts = Factory::default_config(type);
    auto DUNE_UNUSED(aa) = Factory::create();
    auto DUNE_UNUSED(bb) = Factory::create(type);
    auto DUNE_UNUSED(cc) = Factory::create(opts);
  }
}
