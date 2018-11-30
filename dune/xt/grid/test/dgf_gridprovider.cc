// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2018)

#include <dune/xt/common/test/main.hxx>

#include "provider.hh"


struct DgfGridProvider : public GridProviderBase<TESTGRIDTYPE>
{
  typedef GridProviderBase<TESTGRIDTYPE> BaseType;

  void check_layers()
  {
    BaseType::check_layers(Dune::XT::Grid::make_dgf_grid<TESTGRIDTYPE>());
  }

  void check_visualize()
  {
    BaseType::check_visualize(Dune::XT::Grid::make_dgf_grid<TESTGRIDTYPE>());
  }
};

TEST_F(DgfGridProvider, layers)
{
  this->check_layers();
}
TEST_F(DgfGridProvider, visualize)
{
  this->check_visualize();
}
