// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2015)
//   Rene Milk       (2014 - 2015)
//   Tobias Leibner  (2014)

#include "main.hxx"
#include "functions.hh"

#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/global.hh>

#if HAVE_DUNE_GRID
#include <dune/grid/yaspgrid.hh>
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif // HAVE_ALUGRID

struct GlobalLambdaFunctionTest : public DS::FunctionTest<TESTFUNCTIONTYPE>
{
  virtual void check() const
  {
    const FunctionType zero([](DomainType /*x*/) { return RangeType(0); }, 0);
    const FunctionType one([](DomainType /*x*/) { return RangeType(1); }, 0);
    const auto diff = zero - one;
    const auto xx = DomainType(666);
    EXPECT_EQ(zero.evaluate(xx), RangeType(0));
    EXPECT_EQ(one.evaluate(xx), RangeType(1));
  }
};

TEST_F(GlobalLambdaFunctionTest, provides_required_methods)
{
  this->check();
}

#else // HAVE_DUNE_GRID

TEST(DISABLED_GlobalLambdaFunctionTest, provides_required_methods)
{
}

#endif // HAVE_DUNE_GRID
