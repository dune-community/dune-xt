// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2013 - 2015)
//   Rene Milk       (2013 - 2015)
//   Tobias Leibner  (2014)

#include <dune/xt/test/main.hxx>


#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/functions/constant.hh>

#include <dune/grid/yaspgrid.hh>
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif // HAVE_ALUGRID

#include "functions.hh"

using namespace Dune;
using namespace Dune::XT;

struct ConstantFunctionTest : public FunctionTest<TESTFUNCTIONTYPE>
{
  virtual void check() const
  {
    const std::unique_ptr<const FunctionType> function(FunctionType::create(FunctionType::default_config()));
  }
};


TEST_F(ConstantFunctionTest, provides_required_methods)
{
  this->check();
}
