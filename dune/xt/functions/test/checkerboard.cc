// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2016)
//   Rene Milk       (2014 - 2015)
//   Tobias Leibner  (2014 - 2015)

#include <dune/xt/common/test/main.hxx>

#include <memory>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dune/xt/common/exceptions.hh>

#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/functions/checkerboard.hh>

#include "functions.hh"

using namespace Dune;
using namespace Dune::XT;

struct CheckerboardFunctionTest : public FunctionTest<TESTFUNCTIONTYPE>
{
  void check() const
  {
    const std::unique_ptr<const FunctionType> function(FunctionType::create(FunctionType::default_config()));
  }
};

TEST_F(CheckerboardFunctionTest, provides_required_methods)
{
  this->check();
}
