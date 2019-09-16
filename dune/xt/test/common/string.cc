// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2018)
//   René Fritze     (2012 - 2016, 2018)
//   Tobias Leibner  (2014, 2016, 2018)

#include <dune/xt/common/test/main.hxx>

#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/densematrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/densevector.hh>

#include <dune/xt/common/string.hh>
#include <dune/xt/common/logstreams.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/fmatrix.hh>

using namespace Dune::XT::Common;
using namespace std;

// define types and expected values for the typed tests
template <class MatrixType>
struct MatrixStringTestDouble : public ::testing::Test
{
  void check() const
  {
    EXPECT_EQ("[1 2; 3 4]", to_string(from_string<MatrixType>("[1.0 2; 3.0 4]")));
    EXPECT_THROW(from_string<MatrixType>("[1 2; 3 4]", 3, 3), Exceptions::conversion_error);
  }
};

template <class MatrixType>
struct MatrixStringTestChar : public ::testing::Test
{
  void check() const
  {
    EXPECT_EQ("[1 2; 3 4]", to_string(from_string<MatrixType>("[1 2; 3 4]")));
    EXPECT_THROW(from_string<MatrixType>("[1 2; 3 4]", 3, 3), Exceptions::conversion_error);
  }
};

template <class VectorType>
struct VectorStringTestDouble : public ::testing::Test
{
  void check() const
  {
    EXPECT_EQ("[1 2 3]", to_string(from_string<VectorType>("[1.0 2 3.0]")));
    EXPECT_THROW(from_string<VectorType>("[1.0 2 3.0]", 4), Exceptions::conversion_error);
  }
};

template <class VectorType>
struct VectorStringTestInt : public ::testing::Test
{
  void check() const
  {
    EXPECT_EQ("[1 2 3]", to_string(from_string<VectorType>("[1 2 3]")));
    EXPECT_THROW(from_string<VectorType>("[1 2 3]", 4), Exceptions::conversion_error);
  }
};

typedef testing::
    Types<std::vector<double>, Dune::DynamicVector<double>, Dune::FieldVector<double, 3>, FieldVector<double, 3>>
        VectorTypesDouble;

typedef testing::Types<std::vector<int>, Dune::DynamicVector<int>, Dune::FieldVector<int, 3>, FieldVector<int, 3>>
    VectorTypesInt;

typedef testing::Types<Dune::DynamicMatrix<double>, FieldMatrix<double, 2, 2>, Dune::FieldMatrix<double, 2, 2>>
    MatrixTypesDouble;

typedef testing::Types<Dune::DynamicMatrix<char>, FieldMatrix<char, 2, 2>, Dune::FieldMatrix<char, 2, 2>>
    MatrixTypesChar;

// from_string/to_string tests for vector and matrix types
TYPED_TEST_CASE(MatrixStringTestDouble, MatrixTypesDouble);
TYPED_TEST(MatrixStringTestDouble, CheckDouble)
{
  this->check();
}

TYPED_TEST_CASE(MatrixStringTestChar, MatrixTypesChar);
TYPED_TEST(MatrixStringTestChar, CheckChar)
{
  this->check();
}

TYPED_TEST_CASE(VectorStringTestDouble, VectorTypesDouble);
TYPED_TEST(VectorStringTestDouble, CheckDouble)
{
  this->check();
}

TYPED_TEST_CASE(VectorStringTestInt, VectorTypesInt);
TYPED_TEST(VectorStringTestInt, CheckInt)
{
  this->check();
}

// Additional from_string/to_string tests
GTEST_TEST(StringTest, ConvertTo)
{
  EXPECT_EQ("9", to_string(from_string<int>("9")));
  EXPECT_EQ("P", to_string(from_string<char>("P")));
  EXPECT_EQ(double(0.1), from_string<double>(to_string<double>(0.1)));
  EXPECT_EQ("0.100006", to_string(from_string<double>("0.1000055511151231257827021181583404541015625")));
  EXPECT_EQ(1e-14, from_string<double>("1e-14"));
  EXPECT_EQ(1e-14, from_string<double>("1E-14"));
  EXPECT_EQ("1", to_string(from_string<bool>("1")));
  EXPECT_EQ("0", to_string(from_string<bool>("0")));

  EXPECT_EQ("-1", to_string(from_string<long>("-1")));
  typedef std::complex<double> Complex;
  EXPECT_EQ("-1.3+0i", to_string(from_string<Complex>("-1.3")));
  EXPECT_EQ("-1.3+1i", to_string(from_string<Complex>("-1.3+1i")));
  EXPECT_EQ("1.3-1i", to_string(from_string<Complex>("1.3-1i")));
  EXPECT_EQ("0-1i", to_string(from_string<Complex>("0-1i")));
}

GTEST_TEST(StringTest, ConvertFrom)
{
  EXPECT_EQ(9, from_string<int>("9"));
  EXPECT_EQ(0, from_string<int>("0"));
  EXPECT_EQ('p', from_string<char>(to_string('p')));
  EXPECT_EQ(-1, from_string<char>(to_string(char(-1))));
  EXPECT_THROW(from_string<char>("sd"), Exceptions::conversion_error);
  EXPECT_EQ(true, from_string<bool>("1"));
  EXPECT_EQ(true, from_string<bool>("true"));
  EXPECT_EQ(true, from_string<bool>("True"));
  EXPECT_EQ(false, from_string<bool>("0"));
  EXPECT_EQ(false, from_string<bool>("false"));
  EXPECT_EQ(false, from_string<bool>("False"));
  EXPECT_THROW(from_string<int>(""), Exceptions::conversion_error);
  typedef std::complex<double> Complex;
  EXPECT_EQ(Complex(1, -1), from_string<Complex>("+1-1i"));
  EXPECT_EQ(Complex(1, 0), from_string<Complex>("+1"));
  EXPECT_EQ(Complex(0, -1), from_string<Complex>("0-1i"));
}

// Hex, whitespacify, tokenize, stringFromTime tests
GTEST_TEST(StringTest, Hex)
{
  EXPECT_GT(boost::lexical_cast<HexToString<unsigned long>>(&cout), 0u);
  EXPECT_EQ(boost::lexical_cast<HexToString<unsigned long>>("0x00000F"), 15u);
}

GTEST_TEST(StringTest, Whitespace)
{
  EXPECT_EQ("---------", whitespaceify("\t\t\t\t\t\t\t\t\t", '-'));
}

GTEST_TEST(StringTest, Trim)
{
  std::vector<std::string> v{"  ---- ", "----  ", "\t----"};
  trim(v);
  for (const auto& str : v)
    EXPECT_EQ(str, "----");
}

GTEST_TEST(StringTest, Tokenizer)
{
  const string seps(" \t;");
  const string msg("a t\tkk;;g");
  const vector<string> tokens_default = {"a", "t", "kk", "", "g"};
  const vector<string> tokens_compressed = {"a", "t", "kk", "g"};
  EXPECT_EQ(tokens_default, tokenize(msg, seps, boost::algorithm::token_compress_off));
  EXPECT_EQ(tokens_compressed, tokenize(msg, seps, boost::algorithm::token_compress_on));
  const string num_msg("-1 2;;4");
  vector<int> numbers_default = {-1, 2, 0, 4};
  vector<int> numbers_compressed = {-1, 2, 4};
  EXPECT_EQ(numbers_default, tokenize<int>(num_msg, seps, boost::algorithm::token_compress_off));
  EXPECT_EQ(numbers_compressed, tokenize<int>(num_msg, seps, boost::algorithm::token_compress_on));
}

GTEST_TEST(StringTest, TimeString)
{
  string ts = stringFromTime(-1);
}
