// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2015 - 2016, 2018)
//   Tobias Leibner  (2016)

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/common/print.hh>

using namespace Dune::XT::Common;
using namespace std;

GTEST_TEST(OutputIterator, All)
{
  const vector<int> ints{0, 1, 2};
  const string s_ints("0,1,2");
  const vector<string> strings{"a", "b", "c"};
  const string s_strings("a,b,c");
  stringstream stream;
  std::copy(ints.begin(), ints.end(), PrefixOutputIterator<int>(stream, ","));
  EXPECT_EQ(stream.str(), s_ints);
  stream.str("");
  std::copy(strings.begin(), strings.end(), PrefixOutputIterator<string>(stream, ","));
  EXPECT_EQ(stream.str(), s_strings);
}
