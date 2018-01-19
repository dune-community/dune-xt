// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2015 - 2016, 2018)
//   Tobias Leibner  (2015 - 2016)

#include <dune/xt/common/test/main.hxx>

#include <memory>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/functions/spe10/model1.hh>
#include <dune/xt/functions/spe10/model2.hh>

#include "functions.hh"

using namespace Dune;
using namespace Dune::XT;

struct Spe10Model2FunctionTest : public FunctionTest<TESTFUNCTIONTYPE>
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
