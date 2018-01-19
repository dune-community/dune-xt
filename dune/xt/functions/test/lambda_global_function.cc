// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014 - 2016, 2018)
//   Tobias Leibner  (2014 - 2017)

#include <dune/xt/common/test/main.hxx>

#include <memory>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/functions/lambda/global-function.hh>
#include <dune/xt/functions/interfaces.hh>

#include "functions.hh"

using namespace Dune;
using namespace Dune::XT;

struct GlobalLambdaFunctionTest : public FunctionTest<TESTFUNCTIONTYPE>
{
  virtual void check() const
  {
    const FunctionType zero([](DomainType /*x*/, const Common::Parameter& /*mu*/ = {}) { return RangeType(0); }, 0);
    const FunctionType one([](DomainType /*x*/, const Common::Parameter& /*mu*/ = {}) { return RangeType(1); }, 0);
    const auto diff = zero - one;
    const auto xx = DomainType(666);
    EXPECT_EQ(zero.evaluate(xx), RangeType(0));
    EXPECT_EQ(one.evaluate(xx), RangeType(1));
  }
};
