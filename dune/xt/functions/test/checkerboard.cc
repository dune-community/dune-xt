// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014 - 2016)
//   Tobias Leibner  (2014 - 2016)

#include <dune/xt/common/test/main.hxx>

#include <memory>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/functions/checkerboard.hh>

#include "functions.hh"

using namespace Dune;
using namespace Dune::XT;

struct CheckerboardFunctionTest : public FunctionTest<TESTFUNCTIONTYPE>
{
  void check() const
  {
    typedef typename FunctionType::DomainType DomainType;
    typedef typename FunctionType::RangeType RangeType;
    static const size_t dimDomain = FunctionType::dimDomain;
    const std::unique_ptr<const FunctionType> function(FunctionType::create(FunctionType::default_config()));
    const FunctionType function2(DomainType(0),
                                 DomainType(1),
                                 FieldVector<size_t, dimDomain>(2),
                                 std::vector<RangeType>{RangeType(1),
                                                        RangeType(2),
                                                        RangeType(3),
                                                        RangeType(4),
                                                        RangeType(5),
                                                        RangeType(6),
                                                        RangeType(7),
                                                        RangeType(8)});
  }
};

TEST_F(CheckerboardFunctionTest, provides_required_methods)
{
  this->check();
}
