// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "main.hxx"

#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/stuff/functions.hh>

#if HAVE_DUNE_GRID
#include <dune/grid/yaspgrid.hh>
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif // HAVE_ALUGRID


/* we just take the constant function as a container for the types we need */
/* since this one always exists for all combinations */
struct FunctionsTest : public ::testing::Test
{
  typedef typename TESTFUNCTIONTYPE ConstantFunctionType;
  typedef typename ConstantFunctionType::EntityType E;
  typedef typename ConstantFunctionType::DomainFieldType D;
  static const size_t d = ConstantFunctionType::dimDomain;
  typedef typename ConstantFunctionType::RangeFieldType R;
  static const size_t r  = ConstantFunctionType::dimRange;
  static const size_t rC = ConstantFunctionType::dimRangeCols;
  typedef Dune::Stuff::FunctionsProvider<E, D, d, R, r, rC> FunctionsProvider;
  typedef typename FunctionsProvider::InterfaceType InterfaceType;

  void check() const
  {
    for (const auto& type : FunctionsProvider::available()) {
      const Dune::Stuff::Common::Configuration config = FunctionsProvider::default_config(type);
      try {
        const std::unique_ptr<InterfaceType> function = FunctionsProvider::create(type, config);
      } catch (Dune::Stuff::Exceptions::spe10_data_file_missing&) {
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
