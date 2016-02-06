// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:

#include <dune/xt/common/test/main.hxx>

#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID

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

#else // HAVE_ALUGRID || HAVE_DUNE_ALUGRID


TEST(DISABLED_GmshGridProvider, layers)
{
}
TEST(DISABLED_GmshGridProvider, visualize)
{
}


#endif // HAVE_ALUGRID || HAVE_DUNE_ALUGRID
