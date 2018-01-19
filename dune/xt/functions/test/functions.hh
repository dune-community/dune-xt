// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2015 - 2016, 2018)
//   Tobias Leibner  (2015 - 2016)

#ifndef DUNE_XT_TEST_FUNCTIONS_HH
#define DUNE_XT_TEST_FUNCTIONS_HH

#include <memory>
#include <string>
#include <type_traits>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/test/gtest/gtest.h>

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
  static const size_t dimRange = FunctionImp::dimRange;
  static const size_t dimRangeCols = FunctionImp::dimRangeCols;
  typedef typename FunctionImp::LocalfunctionType LocalfunctionType;
  typedef typename FunctionImp::DomainType DomainType;
  typedef typename FunctionImp::RangeType RangeType;
  typedef typename FunctionImp::JacobianRangeType JacobianRangeType;
  typedef Functions::
      LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols>
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
    for (auto&& entity : elements(grid.leafGridView()))
      std::unique_ptr<LocalfunctionType> local_func = func.local_function(entity);
    std::string tp = func.type();
    std::string nm = func.name();
  } // ... dynamic_interface_check(...)

}; // class FunctionTest

} // namespace XT
} // namespace Dune

#endif // DUNE_XT_TEST_FUNCTIONS_HH
