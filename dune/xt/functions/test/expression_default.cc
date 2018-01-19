// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2017)
//   Rene Milk       (2010, 2012 - 2016, 2018)
//   Tobias Leibner  (2014 - 2016)

#include <dune/xt/common/test/main.hxx>

#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/xt/grid/grids.hh>
#include <dune/xt/functions/expression.hh>
#include <dune/xt/functions/interfaces.hh>

#include "functions.hh"

using namespace Dune;
using namespace Dune::XT;

struct ExpressionFunctionTest : public FunctionTest<TESTFUNCTIONTYPE>
{
  virtual void check() const
  {
    Common::Configuration config = FunctionType::default_config();
    const std::unique_ptr<const FunctionType> function(FunctionType::create(config));
    config["expression"] = "[2*x[0] 3*x[0] 4*x[0]; 1 sin(x[0]) 0; cos(x[0]) x[0] 0]";
    if (dimRangeCols == 1)
      config["gradient"] = "[2 0 0; 0 0 0; -sin(x[0]) 0 0]";
    else {
      config["gradient.0"] = "[2 0 0; 0 0 0; -sin(x[0]) 0 0]";
      config["gradient.1"] = "[3 0 0; cos(x[0]) 0 0; 1 0 0]";
      config["gradient.2"] = "[4 0 0; 0 0 0; 0 0 0]";
    }
    const std::unique_ptr<const FunctionType> function2(FunctionType::create(config));
    const std::unique_ptr<const FunctionType> function3(
        new FunctionType("x", "sin(x[0])", 3, FunctionType::static_id(), {"cos(x[0])", "0", "0"}));
  }
};

TEST_F(ExpressionFunctionTest, provides_required_methods)
{
  this->check();
}
