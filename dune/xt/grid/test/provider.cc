// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2016)
//   Rene Milk       (2012 - 2016)
//   Tobias Leibner  (2014 - 2016)

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/grid/provider/cube.hh>
#include <dune/xt/grid/provider/dgf.hh>
#include <dune/xt/grid/provider/gmsh.hh>

#include <dune/xt/grid/test/provider.hh>


struct TestGridProvider : public GridProviderBase<Dune::XT::Grid::Providers::TESTGRIDPROVIDERTYPE<TESTGRIDTYPE>>
{
};

TEST_F(TestGridProvider, is_default_creatable)
{
  this->is_default_creatable();
}
TEST_F(TestGridProvider, fulfills_const_interface)
{
  this->const_interface();
}
TEST_F(TestGridProvider, is_visualizable)
{
  this->visualize();
}
TEST_F(TestGridProvider, fulfills_non_const_interface)
{
  this->non_const_interface();
}
