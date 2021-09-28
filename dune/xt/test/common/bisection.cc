// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2019 - 2020)

#include <dune/xt/test/main.hxx> // <- Has to come first, includes the config.h!
#include <dune/xt/common/bisect.hh>
#include <dune/xt/test/common/float_cmp.hh>

using namespace Dune::XT::Common;

GTEST_TEST(find_largest_by_bisection_test, works_for_left_smaller_right)
{
  const auto result =
      find_largest_by_bisection(/*left=*/0, /*right=*/1, /*condition=*/[&](const auto& x) { return x < 0.5; });
  DXTC_EXPECT_FLOAT_EQ(0.5, result);
}

GTEST_TEST(find_largest_by_bisection_test, works_for_right_smaller_left)
{
  const auto result =
      find_largest_by_bisection(/*left=*/1, /*right=*/0, /*condition=*/[&](const auto& x) { return x < 0.5; });
  DXTC_EXPECT_FLOAT_EQ(0.5, result);
}

GTEST_TEST(find_largest_by_bisection_test, breaks_for_broken_condition)
{
  EXPECT_THROW(
      find_largest_by_bisection(/*left=*/0, /*right=*/1, /*condition=*/[&](const auto& /*x*/) { return false; }),
      Exceptions::bisection_error);
}
