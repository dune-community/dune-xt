// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   René Fritze     (2015 - 2018)
//   Tim Keil        (2018)
//   Tobias Leibner  (2016 - 2018)

#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/xt/common/test/main.hxx>

#include <vector>

#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/test/float_cmp.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/fmatrix.hh>

using namespace Dune;
using XT::Common::create;
using XT::Common::FloatCmp::Style;

// add operator+= for std::array and std::vector
template <typename T, size_t N>
std::array<T, N>& operator+=(std::array<T, N>& arr, const std::array<T, N>& other)
{
  std::transform(arr.begin(), arr.end(), other.begin(), arr.begin(), std::plus<T>());
  return arr;
}

template <typename T>
std::vector<T>& operator+=(std::vector<T>& vec, const std::vector<T>& other)
{
  DXT_ASSERT(vec.size() == other.size());
  std::transform(vec.begin(), vec.end(), other.begin(), vec.begin(), std::plus<T>());
  return vec;
}


template <class T, class S>
typename std::enable_if<XT::Common::is_arithmetic<T>::value, T>::type make_number(const S& number)
{
  return T(number);
}

template <class T, class S>
typename std::enable_if<XT::Common::is_complex<T>::value && !XT::Common::is_complex<S>::value, T>::type
make_number(const S& number)
{
  return T(number, number); // <- This has to match the DEFAULT_EPSILON below!
}

template <class T>
typename std::enable_if<XT::Common::is_complex<T>::value, T>::type make_number(const T& number)
{
  return number;
}


template <class T>
typename std::enable_if<XT::Common::is_arithmetic<T>::value, T>::type make_type(const T& number)
{
  return number;
}

template <class T>
typename std::enable_if<XT::Common::is_complex<T>::value, T>::type make_type(const T& number)
{
  return number; // <- No special treatment here, already done above!
}

template <class T, class S>
typename std::enable_if<XT::Common::is_vector<T>::value, T>::type make_type(const S& number)
{
  using Vec = XT::Common::VectorAbstraction<T>;
  if (Vec::has_static_size)
    return Vec::create(Vec::static_size, number);
  else
    return Vec::create(VECSIZE, number);
}

template <class T, class S>
typename std::enable_if<XT::Common::is_matrix<T>::value, T>::type make_type(const S& number)
{
  using Mat = XT::Common::MatrixAbstraction<T>;
  if (Mat::has_static_size)
    return Mat::create(Mat::static_rows, Mat::static_cols, number);
  else
    return Mat::create(VECSIZE, NUMCOLS, number);
}

template <class T,
          bool is_number = (XT::Common::is_arithmetic<T>::value || XT::Common::is_complex<T>::value),
          bool is_vector = XT::Common::is_vector<T>::value>
struct suitable_entry
{
  static T change_an_entry(const T& input, const int value)
  {
    DUNE_UNUSED_PARAMETER(value);
    return input;
  }
};

template <class T>
struct suitable_entry<T, false, true>
{
  static T change_an_entry(const T& input, const int value)
  {
    T ret = input;
    ret[0] = value;
    return ret;
  }
};

template <class T>
struct suitable_entry<T, false, false>
{
  static T change_an_entry(const T& input, const int value)
  {
    T ret = input;
    ret[0][0] = value;
    return ret;
  }
};

template <class T>
struct FloatCmpTest : public testing::Test
{
  typedef typename std::conditional<XT::Common::is_matrix<T>::value,
                                    typename XT::Common::MatrixAbstraction<T>::ScalarType,
                                    typename std::conditional<XT::Common::is_vector<T>::value,
                                                              typename XT::Common::VectorAbstraction<T>::ScalarType,
                                                              T>::type>::type S;
  typedef typename std::conditional<XT::Common::is_matrix<T>::value,
                                    typename XT::Common::MatrixAbstraction<T>::RealType,
                                    typename std::conditional<XT::Common::is_vector<T>::value,
                                                              typename XT::Common::VectorAbstraction<T>::RealType,
                                                              typename Dune::FieldTraits<T>::real_type>::type>::type R;
  static constexpr bool fieldtype_is_float = std::is_floating_point<R>::value;

  FloatCmpTest()
    : zero(make_type<T>(make_number<S>(0)))
    , one_with_one_zero_entry(suitable_entry<T>::change_an_entry(make_type<T>(make_number<S>(1)), 0))
    , one(make_type<T>(make_number<S>(1)))
    , neg_one(make_type<T>(make_number<S>(-1)))
    , epsilon(make_type<T>(make_number<S>(XT::Common::FloatCmp::DEFAULT_EPSILON::value())))
    , eps_plus(make_type<T>(make_number<S>(XT::Common::FloatCmp::DEFAULT_EPSILON::value() * 1.1)))
    , eps_minus(make_type<T>(make_number<S>(XT::Common::FloatCmp::DEFAULT_EPSILON::value() * 0.9)))
    , two(make_type<T>(make_number<S>(2)))
    , test_config(DXTC_CONFIG.sub("test_common_float_cmp"))
  {
    one_plus_eps_minus_ = one;
    one_plus_eps_minus_ += eps_minus;
    one_plus_eps_plus_ = one;
    one_plus_eps_plus_ += eps_plus;
  }

  const T zero;
  const T one_with_one_zero_entry;
  const T one;
  const T neg_one;
  const T epsilon;
  const T eps_plus;
  const T eps_minus;
  const T two;
  T one_plus_eps_minus_;
  T one_plus_eps_plus_;
  const typename XT::Common::Configuration test_config;

  void check_eq()
  {
    TEST_DXTC_EXPECT_FLOAT_EQ(zero, zero);
    EXPECT_FALSE(FLOATCMP_EQ(zero, one));
    TEST_DXTC_EXPECT_FLOAT_EQ(one, one_plus_eps_minus_);
    EXPECT_FALSE(FLOATCMP_EQ(one, two));

    if (test_config["cmpstyle_is_relative"] == "true" && fieldtype_is_float)
      EXPECT_FALSE(FLOATCMP_EQ(zero, eps_minus));
    else
      TEST_DXTC_EXPECT_FLOAT_EQ(zero, eps_minus);

    if (test_config["cmpstyle_is_numpy"] == "true" || !fieldtype_is_float)
      TEST_DXTC_EXPECT_FLOAT_EQ(one, one_plus_eps_plus_);
    else
      EXPECT_FALSE(FLOATCMP_EQ(one, one_plus_eps_plus_));

    if (fieldtype_is_float)
      EXPECT_FALSE(FLOATCMP_EQ(zero, eps_plus));
    else
      TEST_DXTC_EXPECT_FLOAT_EQ(zero, eps_plus);
  }

  void check_ne()
  {
    EXPECT_FALSE(FLOATCMP_NE(zero, zero));
    TEST_DXTC_EXPECT_FLOAT_NE(zero, one);
    EXPECT_FALSE(FLOATCMP_NE(one, one_plus_eps_minus_));
    TEST_DXTC_EXPECT_FLOAT_NE(one, two);
    if (test_config["cmpstyle_is_relative"] == "true" && fieldtype_is_float)
      TEST_DXTC_EXPECT_FLOAT_NE(zero, eps_minus);
    else
      EXPECT_FALSE(FLOATCMP_NE(zero, eps_minus));

    if (test_config["cmpstyle_is_numpy"] == "true" || !fieldtype_is_float)
      EXPECT_FALSE(FLOATCMP_NE(one, one_plus_eps_plus_));
    else
      TEST_DXTC_EXPECT_FLOAT_NE(one, one_plus_eps_plus_);

    if (fieldtype_is_float)
      TEST_DXTC_EXPECT_FLOAT_NE(zero, eps_plus);
    else
      EXPECT_FALSE(FLOATCMP_NE(zero, eps_plus));
  }

  void check_gt()
  {
    EXPECT_FALSE(FLOATCMP_GT(zero, zero));
    TEST_DXTC_EXPECT_FLOAT_GT(one, zero);
    EXPECT_FALSE(FLOATCMP_GT(one_plus_eps_minus_, one));
    TEST_DXTC_EXPECT_FLOAT_GT(two, one);

    if (test_config["cmpstyle_is_relative"] == "true" && fieldtype_is_float)
      TEST_DXTC_EXPECT_FLOAT_GT(eps_minus, zero);
    else
      EXPECT_FALSE(FLOATCMP_GT(eps_minus, zero));

    if (test_config["cmpstyle_is_numpy"] == "true" || !fieldtype_is_float)
      EXPECT_FALSE(FLOATCMP_GT(one_plus_eps_plus_, one));
    else
      TEST_DXTC_EXPECT_FLOAT_GT(one_plus_eps_plus_, one);

    if (fieldtype_is_float)
      TEST_DXTC_EXPECT_FLOAT_GT(eps_plus, zero);
    else
      EXPECT_FALSE(FLOATCMP_GT(eps_plus, zero));
  }

  void check_lt()
  {
    EXPECT_FALSE(FLOATCMP_LT(zero, zero));
    TEST_DXTC_EXPECT_FLOAT_LT(zero, one);
    EXPECT_FALSE(FLOATCMP_LT(one, one_plus_eps_minus_));
    TEST_DXTC_EXPECT_FLOAT_LT(one, two);

    if (test_config["cmpstyle_is_relative"] == "true" && fieldtype_is_float)
      TEST_DXTC_EXPECT_FLOAT_LT(zero, eps_minus);
    else
      EXPECT_FALSE(FLOATCMP_LT(zero, eps_minus));

    if (test_config["cmpstyle_is_numpy"] == "true" || !fieldtype_is_float)
      EXPECT_FALSE(FLOATCMP_LT(one, one_plus_eps_plus_));
    else
      TEST_DXTC_EXPECT_FLOAT_LT(one, one_plus_eps_plus_);

    if (fieldtype_is_float)
      TEST_DXTC_EXPECT_FLOAT_LT(zero, eps_plus);
    else
      EXPECT_FALSE(FLOATCMP_LT(zero, eps_plus));
  }

  void check_ge()
  {
    TEST_DXTC_EXPECT_FLOAT_GE(zero, zero);
    TEST_DXTC_EXPECT_FLOAT_GE(eps_minus, zero);
    TEST_DXTC_EXPECT_FLOAT_GE(eps_plus, zero);
    TEST_DXTC_EXPECT_FLOAT_GE(one, zero);
    TEST_DXTC_EXPECT_FLOAT_GE(one_plus_eps_minus_, one);
    TEST_DXTC_EXPECT_FLOAT_GE(one_plus_eps_plus_, one);
    TEST_DXTC_EXPECT_FLOAT_GE(two, one);
    TEST_DXTC_EXPECT_FLOAT_GE(one, one_with_one_zero_entry);
  }

  void check_le()
  {
    TEST_DXTC_EXPECT_FLOAT_LE(zero, zero);
    TEST_DXTC_EXPECT_FLOAT_LE(zero, eps_minus);
    TEST_DXTC_EXPECT_FLOAT_LE(zero, eps_plus);
    TEST_DXTC_EXPECT_FLOAT_LE(zero, one);
    TEST_DXTC_EXPECT_FLOAT_LE(one, one_plus_eps_minus_);
    TEST_DXTC_EXPECT_FLOAT_LE(one, one_plus_eps_plus_);
    TEST_DXTC_EXPECT_FLOAT_LE(one, two);
    TEST_DXTC_EXPECT_FLOAT_LE(one_with_one_zero_entry, one);
  }

  void negative_numbers()
  {
    TEST_DXTC_EXPECT_FLOAT_LT(neg_one, one);
    TEST_DXTC_EXPECT_FLOAT_LT(neg_one, zero);
    EXPECT_FALSE(FLOATCMP_LT(neg_one, neg_one));
    TEST_DXTC_EXPECT_FLOAT_GT(one, neg_one);
    EXPECT_FALSE(FLOATCMP_GT(neg_one, neg_one));
    EXPECT_FALSE(FLOATCMP_LT(neg_one, neg_one));
    TEST_DXTC_EXPECT_FLOAT_NE(neg_one, two);
    EXPECT_FALSE(FLOATCMP_NE(neg_one, neg_one));
    EXPECT_FALSE(FLOATCMP_EQ(zero, neg_one));
    TEST_DXTC_EXPECT_FLOAT_EQ(neg_one, neg_one);
  }
}; // struct FloatCmpBase


TYPED_TEST_CASE(FloatCmpTest, TESTTYPES);
TYPED_TEST(FloatCmpTest, eq)
{
  this->check_eq();
}

TYPED_TEST(FloatCmpTest, ne)
{
  this->check_ne();
}

TYPED_TEST(FloatCmpTest, gt)
{
  this->check_gt();
}

TYPED_TEST(FloatCmpTest, lt)
{
  this->check_lt();
}

TYPED_TEST(FloatCmpTest, ge)
{
  this->check_ge();
}

TYPED_TEST(FloatCmpTest, le)
{
  this->check_le();
}

TYPED_TEST(FloatCmpTest, neg)
{
  if (std::is_signed<typename XT::Common::VectorAbstraction<TypeParam>::ScalarType>::value)
    this->negative_numbers();
}
