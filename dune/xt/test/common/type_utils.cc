// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   René Fritze     (2012 - 2016, 2018)
//   Tobias Leibner  (2014, 2016)

#include <dune/xt/common/test/main.hxx>

#include <iostream>

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/validation.hh>

DUNE_XT_COMMON_TYPENAME(Dune::MPIHelper)

using namespace Dune::XT::Common;
using namespace std;

GTEST_TEST(Typename, Knowns)
{
  EXPECT_EQ(Typename<unsigned long>::value(), "unsigned long");
  EXPECT_EQ(Typename<unsigned long>::value(), "unsigned long");
  EXPECT_EQ(Typename<int>::value(), "int");
  EXPECT_EQ(Typename<double>::value(), "double");
  EXPECT_EQ(Typename<float>::value(), "float");
  EXPECT_EQ(get_typename(1.f), "float");
  EXPECT_EQ(Typename<long>::value(), "long");
  EXPECT_EQ(Typename<unsigned int>::value(), "unsigned int");
  EXPECT_EQ(Typename<unsigned long>::value(), "unsigned long");
  EXPECT_EQ(Typename<char>::value(), "char");
}

GTEST_TEST(Typename, Unknowns)
{
  EXPECT_NE(Typename<Dune::Exception>::value(), string());
}

GTEST_TEST(Typename, Extended)
{
  EXPECT_EQ(Typename<Dune::MPIHelper>::value(), string("Dune::MPIHelper"));
  {
    using namespace Dune;
    EXPECT_EQ(Typename<MPIHelper>::value(), string("Dune::MPIHelper"));
  }
  typedef ValidateAny<int> Valid;
  Valid v;
  EXPECT_EQ(get_typename(v), string("Parameter::ValidateAny<T>"));
}
