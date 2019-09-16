// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014, 2016 - 2017)
//   René Fritze     (2012 - 2013, 2015 - 2018)
//   Tobias Leibner  (2014, 2016)

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/vector.hh>
#include <dune/xt/common/vector_statistics.hh>

using namespace Dune::XT::Common;
using namespace std;

GTEST_TEST(STDDEV, Vector)
{
  const auto vec = VectorAbstraction<VECTOR_TYPE>::create(VECTOR_SIZE, 0);
  EXPECT_DOUBLE_EQ(standard_deviation(vec), 0);
}
