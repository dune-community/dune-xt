// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/xt/test/main.hxx>

#include <dune/xt/grid/provider/cube.hh>
#include <dune/xt/grid/provider/dgf.hh>
#include <dune/xt/grid/provider/gmsh.hh>

#include "grid_provider.hh"


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
