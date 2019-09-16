// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014, 2016 - 2017)
//   René Fritze     (2012 - 2013, 2015 - 2016, 2018)
//   Tobias Leibner  (2014, 2016, 2018)

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/common/debug.hh>

using namespace Dune::XT::Common;

struct Dummy
{
  static std::string name()
  {
    return std::string(__CLASS__);
  }
};

GTEST_TEST(debug, main)
{
#ifndef NDEBUG
  EXPECT_THROW(DXT_ASSERT(false), Dune::XT::Common::Exceptions::debug_assertion);
#else
  EXPECT_NO_THROW(DXT_ASSERT(false));
#endif

  EXPECT_THROW(Dummy::name(), Dune::NotImplemented);
  EXPECT_THROW(__CLASS__, Dune::NotImplemented);
}