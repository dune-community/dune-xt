// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2015)
//   Rene Milk       (2015)

#ifndef DUNE_XT_TEST_FUNCTIONS_HH
#define DUNE_XT_TEST_FUNCTIONS_HH

#include <string>
#include <type_traits>
#include <memory>

#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/test/gtest/gtest.h>
#include <dune/xt/functions/interfaces.hh>

namespace Dune {
namespace XT {

template <class FunctionImp>
class FunctionTest : public ::testing::Test
{
protected:
  // required types and static members
  typedef FunctionImp FunctionType;
  typedef typename FunctionImp::EntityType EntityType;
  typedef typename FunctionImp::DomainFieldType DomainFieldType;
  static const size_t dimDomain = FunctionImp::dimDomain;
  typedef typename FunctionImp::RangeFieldType RangeFieldType;
  static const size_t dimRange     = FunctionImp::dimRange;
  static const size_t dimRangeCols = FunctionImp::dimRangeCols;
  typedef typename FunctionImp::LocalfunctionType LocalfunctionType;
  typedef typename FunctionImp::DomainType DomainType;
  typedef typename FunctionImp::RangeType RangeType;
  typedef typename FunctionImp::JacobianRangeType JacobianRangeType;
  typedef LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols>
      InterfaceType;

  static void static_interface_check()
  {
    // compare static information to interface
    static_assert(std::is_same<EntityType, typename InterfaceType::EntityType>::value, "");
    static_assert(std::is_same<DomainFieldType, typename InterfaceType::DomainFieldType>::value, "");
    static_assert(dimDomain == InterfaceType::dimDomain, "");
    static_assert(std::is_same<RangeFieldType, typename InterfaceType::RangeFieldType>::value, "");
    static_assert(dimRange == InterfaceType::dimRange, "");
    static_assert(dimRangeCols == InterfaceType::dimRangeCols, "");
    static_assert(std::is_same<LocalfunctionType, typename InterfaceType::LocalfunctionType>::value, "");
    static_assert(std::is_same<DomainType, typename InterfaceType::DomainType>::value, "");
    static_assert(std::is_same<RangeType, typename InterfaceType::RangeType>::value, "");
    static_assert(std::is_same<JacobianRangeType, typename InterfaceType::JacobianRangeType>::value, "");
    // required static methods
    std::string static_id = FunctionImp::static_id();
  } // ... static_interface_check(...)

  static void static_create_check()
  {
    Common::Configuration default_cfg = FunctionImp::default_config();
    std::unique_ptr<FunctionImp> func = FunctionImp::create(default_cfg);
  } // ... static_create_check(...)

  template <class GridType>
  void dynamic_interface_check(const FunctionImp& func, GridType& grid) const
  {
#if HAVE_DUNE_GRID
    for (const auto& entity : Common::entityRange(grid.leafGridView()))
      std::unique_ptr<LocalfunctionType> local_func = func.local_function(entity);
#endif
    std::string tp = func.type();
    std::string nm = func.name();
  } // ... dynamic_interface_check(...)

}; // class FunctionTest

} // namespace XT
} // namespace Dune

#endif // DUNE_XT_TEST_FUNCTIONS_HH
