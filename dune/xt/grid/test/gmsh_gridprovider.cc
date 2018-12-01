// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2017)
//   René Fritze     (2012 - 2016, 2018)
//   Tobias Leibner  (2014 - 2016)

#include <dune/xt/common/test/main.hxx>

#if HAVE_DUNE_ALUGRID

#include "provider.hh"


struct GmshGridProvider : public GridProviderBase<TESTGRIDTYPE>
{
  typedef GridProviderBase<TESTGRIDTYPE> BaseType;

  void check_layers()
  {
    BaseType::check_layers(Dune::XT::Grid::make_cube_grid<TESTGRIDTYPE>());
  }

  void check_visualize()
  {
    BaseType::check_visualize(Dune::XT::Grid::make_cube_grid<TESTGRIDTYPE>());
  }
};


TEST_F(GmshGridProvider, layers)
{
  this->check_layers();
}
TEST_F(GmshGridProvider, visualize)
{
  this->check_visualize();
}

#else // HAVE_DUNE_ALUGRID


TEST(DISABLED_GmshGridProvider, layers) {}
TEST(DISABLED_GmshGridProvider, visualize) {}


#endif // HAVE_DUNE_ALUGRID
