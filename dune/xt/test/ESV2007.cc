// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015)
//   Rene Milk       (2015)

#include <dune/xt/test/main.hxx>
#include "functions.hh"

#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/functions/ESV2007.hh>

#include <dune/grid/yaspgrid.hh>
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif // HAVE_ALUGRID

using namespace Dune;
using namespace Dune::XT;

struct ESV2007FunctionTest : public FunctionTest<TESTFUNCTIONTYPE>
{
  virtual void check() const
  {
    const std::unique_ptr<const FunctionType> function(FunctionType::create(FunctionType::default_config()));
  }
};

TEST_F(ESV2007FunctionTest, provides_required_methods)
{
  this->check();
}
