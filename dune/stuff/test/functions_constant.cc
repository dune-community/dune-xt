// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "main.hxx"

#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>

#if HAVE_DUNE_GRID
#include <dune/grid/yaspgrid.hh>
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif // HAVE_ALUGRID

// we need this nasty code generation because the testing::Types< ... > only accepts 50 arguments
// and all combinations of functions and entities and dimensions and fieldtypes would be way too much
struct ConstantFunctionTest : public ::testing::Test
{
  typedef typename TESTFUNCTIONTYPE LocalizableFunctionType;
  typedef typename LocalizableFunctionType::EntityType EntityType;
  typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  typedef typename LocalizableFunctionType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = LocalizableFunctionType::dimDomain;
  typedef typename LocalizableFunctionType::DomainType DomainType;
  typedef typename LocalizableFunctionType::RangeFieldType RangeFieldType;
  static const size_t dimRange     = LocalizableFunctionType::dimRange;
  static const size_t dimRangeCols = LocalizableFunctionType::dimRangeCols;
  typedef typename LocalizableFunctionType::RangeType RangeType;
  typedef typename LocalizableFunctionType::JacobianRangeType JacobianRangeType;

  void check() const
  {
    const std::unique_ptr<const LocalizableFunctionType> function(
        LocalizableFunctionType::create(LocalizableFunctionType::default_config()));
  }
};


TEST_F(ConstantFunctionTest, provides_required_methods)
{
  this->check();
}

#else // HAVE_DUNE_GRID

TEST(DISABLED_ConstantFunctionTest, provides_required_methods)
{
}

#endif // HAVE_DUNE_GRID
