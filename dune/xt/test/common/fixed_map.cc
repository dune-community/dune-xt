// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2014, 2016 - 2017)
//   René Fritze     (2012 - 2018)
//   Tobias Leibner  (2014, 2016 - 2017)

#include <dune/xt/common/test/main.hxx>

#include <string>
#include <array>
#include <initializer_list>
#include <vector>
#include <dune/xt/common/fixed_map.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/ranges.hh>

using namespace Dune::XT::Common;

enum class TestEnum
{
  one,
  two,
  three
};


GTEST_TEST(FixedMapTest, All)
{
  const std::initializer_list<std::pair<std::string, int>> values{{"0", 0}, {"1", 1}, {"2", 2}};

  const FixedMap<std::string, int, 1> too_small(values);
  FixedMap<std::string, int, 3> fits(values);
  const FixedMap<std::string, int, 6> too_big(values);
  EXPECT_EQ(0, too_small[to_string(0)]);
  EXPECT_THROW(too_small["1"], Dune::RangeError);
  for (int i : {0, 1, 2}) {
    EXPECT_EQ(i, too_big[to_string(i)]);
    EXPECT_NE(too_big.find(to_string(i)), too_big.end());
  }
  for (int i : {3, 4, 5}) {
    EXPECT_EQ(int(), too_big[std::string()]);
    EXPECT_EQ(too_big.find(to_string(i)), too_big.end());
  }
  auto size = fits.size();
  const auto end = fits.end();
  for (auto& pt : fits) {
    const int value = pt.second;
    const std::string key = pt.first;
    EXPECT_EQ(key, to_string(value));
    EXPECT_NE(fits.find(key), end);
    size--;
  }
  EXPECT_EQ(0, size);
  for (int i : {0, 1, 2}) {
    fits[to_string(i)] += 1;
    EXPECT_EQ(i + 1, fits[to_string(i)]);
  }

  EXPECT_EQ(std::make_pair(std::string("0"), 0), *too_big.begin());
  // EXPECT_DEATH(*too_big.end(), ".*");

  FixedMap<TestEnum, std::string, 1> enum_names = {{TestEnum::one, "one"}};
  EXPECT_EQ(std::string("one"), enum_names[TestEnum::one]);
}
