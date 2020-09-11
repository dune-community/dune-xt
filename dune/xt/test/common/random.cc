// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2016 - 2019)
//   Tobias Leibner  (2018, 2020)

#include <dune/xt/test/main.hxx>

#include <vector>

#include <dune/common/dynvector.hh>
#include <dune/common/densevector.hh>

#include <dune/xt/common/random.hh>
#include <dune/xt/common/float_cmp.hh>

using namespace Dune::XT::Common;
using namespace Dune::XT::Test;
using namespace std;

GTEST_TEST(Init, Random)
{
  typedef RNG_FIELD_TYPE T;
  typedef DefaultRNG<T> RNG;
  //  typedef typename VectorAbstraction<T>::S Scalar;
  const auto lower_bound = init_bound<T>(1);
  // for the int check we need at least on int in between lower and upper bound
  const auto upper_bound = init_bound<T>(3);
  const size_t seed = 0;
  RNG rng1;
  RNG rng2(lower_bound);
  RNG rng3(lower_bound, upper_bound);
  RNG rng4(lower_bound, upper_bound, seed);
  for (auto& r : std::array<RNG, 3>{{rng2, rng3, rng4}}) {
    const auto val = r();
    EXPECT_TRUE(FloatCmp::ge(val, lower_bound));
  }
  for (auto& r : std::array<RNG, 2>{{rng3, rng4}}) {
    const auto val = r();
    EXPECT_TRUE(FloatCmp::le(val, upper_bound));
  }
  {
    RNG rng_a(lower_bound, upper_bound, seed);
    RNG rng_b(lower_bound, upper_bound, seed);
    // same seed -> same values
    for ([[maybe_unused]] auto i : value_range(1000))
      EXPECT_EQ(rng_a(), rng_b());
  }
  DefaultRNG<std::string> str_rng(2);
  std::string rstr = str_rng();
  EXPECT_EQ(rstr.size(), 2);
}
