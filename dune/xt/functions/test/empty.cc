// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018)

/**
  * This file is intended as a starting point for quick testing.
  */

#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_DEBUG_LOGGING 1
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING 1

#include <dune/xt/common/test/main.hxx> // <- Has to come first, includes the config.h!


GTEST_TEST(empty, main)
{
}
