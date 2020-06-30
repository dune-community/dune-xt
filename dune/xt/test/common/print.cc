// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2015 - 2016, 2018 - 2019)
//   Tobias Leibner  (2016, 2020)

#include <dune/xt/test/main.hxx>

#include <dune/xt/common/print.hh>

using namespace Dune::XT::Common;
using namespace std;


GTEST_TEST(print, some_types)
{
  std::cout << print(int(1)) << std::endl;
  std::cout << print(double(2)) << std::endl;
  std::cout << print(std::string("foo")) << std::endl;
  std::cout << print(std::vector<double>({1, 2})) << std::endl;
}


GTEST_TEST(repr, some_types)
{
  std::cout << repr(int(1)) << std::endl;
  std::cout << repr(double(2)) << std::endl;
  std::cout << repr(std::string("foo")) << std::endl;
  std::cout << repr(std::vector<double>({1, 2})) << std::endl;
}


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
