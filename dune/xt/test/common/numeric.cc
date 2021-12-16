// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019 - 2020)

#include <dune/xt/test/main.hxx>
#include <vector>
#include <dune/xt/common/numeric.hh>

GTEST_TEST(MoveIfTest, All)
{
  std::vector<double> vals1(10, 1.);
  std::vector<double> vals2(10);
  // fill vals2 with values [1, 2, ..., 10]
  std::iota(vals2.begin(), vals2.end(), 1);
  double expected1 = (10. * 11.) / 2.; // sum of numbers up to 10
  double expected2 = 3628800.; // 10!
  double expected3 = 39916800.; // 11!
  EXPECT_EQ(Dune::XT::Common::reduce(vals2.begin(), vals2.end(), 0.), expected1);
  EXPECT_EQ(Dune::XT::Common::transform_reduce(vals1.begin(), vals1.end(), vals2.begin(), 0.), expected1);
  EXPECT_EQ(Dune::XT::Common::reduce(vals2.begin(), vals2.end(), 0., std::plus<>()), expected1);
  EXPECT_EQ(Dune::XT::Common::transform_reduce(
                vals1.begin(), vals1.end(), vals2.begin(), 0., std::plus<>(), std::multiplies<>()),
            expected1);
  EXPECT_EQ(Dune::XT::Common::reduce(vals2.begin(), vals2.end(), 1., std::multiplies<>()), expected2);
  EXPECT_EQ(Dune::XT::Common::transform_reduce(
                vals1.begin(), vals1.end(), vals2.begin(), 1., std::multiplies<>(), std::plus<>()),
            expected3);
}
