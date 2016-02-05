// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2012 - 2016)
//   Rene Milk       (2012 - 2016)
//   Tobias Leibner  (2014 - 2016)

#include <dune/xt/common/test/main.hxx>

#include "provider.hh"


struct DgfGridProvider : public GridProviderBase<TESTGRIDTYPE>
{
  typedef GridProviderBase<TESTGRIDTYPE> BaseType;

  void const_interface()
  {
    BaseType::const_interface(Dune::XT::Grid::make_dgf_grid<TESTGRIDTYPE>());
  }
};


TEST_F(DgfGridProvider, const_interface)
{
  this->const_interface();
}
// TEST_F(TestGridProvider, fulfills_const_interface)
//{
//  this->const_interface();
//}
// TEST_F(TestGridProvider, is_visualizable)
//{
//  this->visualize();
//}
// TEST_F(TestGridProvider, fulfills_non_const_interface)
//{
//  this->non_const_interface();
//}
