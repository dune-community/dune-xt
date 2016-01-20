// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2013 - 2015)
//   Rene Milk       (2013 - 2015)
//   Tobias Leibner  (2014)

#include <dune/xt/test/main.hxx>
#include "functions.hh"

#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/xt/functions.hh>

#if HAVE_DUNE_GRID
#include <dune/grid/yaspgrid.hh>
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif // HAVE_ALUGRID

using namespace Dune;
using namespace Dune::XT;

/* we just take the constant function as a container for the types we need */
/* since this one always exists for all combinations */
struct FunctionsTest : public FunctionTest<TESTFUNCTIONTYPE>
{
  typedef FunctionsProvider<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols>
      FunctionsProvider;
  typedef typename FunctionsProvider::InterfaceType InterfaceType;

  virtual void check() const
  {
    for (const auto& type : FunctionsProvider::available()) {
      const Common::Configuration config = FunctionsProvider::default_config(type);
      try {
        const std::unique_ptr<InterfaceType> function = FunctionsProvider::create(type, config);
      } catch (Exceptions::spe10_data_file_missing&) {
      }
    }
  }
};


TEST_F(FunctionsTest, provides_required_methods)
{
  this->check();
}

#else // HAVE_DUNE_GRID

TEST(DISABLED_FunctionsTest, provides_required_methods)
{
}

#endif // HAVE_DUNE_GRID
