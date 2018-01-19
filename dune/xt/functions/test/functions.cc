// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2013 - 2018)
//   Tobias Leibner  (2014 - 2016)

#include <dune/xt/common/test/main.hxx>

#include <memory>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/functions.hh>
#include <dune/xt/common/test/gtest/gtest.h>

#include "functions.hh"

using namespace Dune;
using namespace Dune::XT;
using namespace Dune::XT::Functions;

/* we just take the constant function as a container for the types we need */
/* since this one always exists for all combinations */
struct FunctionsTest : public FunctionTest<TESTFUNCTIONTYPE>
{
  typedef FunctionsFactory<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols>
      FunctionsProviderType;
  typedef typename FunctionsProviderType::InterfaceType InterfaceType;

  virtual void check() const
  {
    for (const auto& type : FunctionsProviderType::available()) {
      const Common::Configuration config = FunctionsProviderType::default_config(type);
      try {
        const std::unique_ptr<InterfaceType> function = FunctionsProviderType::create(type, config);
      } catch (Exceptions::spe10_data_file_missing&) {
      }
    }
  }
};

TEST_F(FunctionsTest, creation_and_evalution)
{
  this->check();
}
