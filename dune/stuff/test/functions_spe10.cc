// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "main.hxx"
#include "functions.hh"

#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/spe10.hh>
#include <dune/stuff/functions/spe10model2.hh>

#if HAVE_DUNE_GRID
#include <dune/grid/yaspgrid.hh>
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif // HAVE_ALUGRID

struct Spe10Model2FunctionTest : public DS::FunctionTest<TESTFUNCTIONTYPE>
{
  void check() const
  {
    const FunctionType zero;
  }
};


TEST_F(Spe10Model2FunctionTest, provides_required_methods)
{
  this->check();
}

#else // HAVE_DUNE_GRID

TEST(DISABLED_Spe10Model2FunctionTest, provides_required_methods)
{
}

#endif // HAVE_DUNE_GRID
