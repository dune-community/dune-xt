// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2018, 2020)

#include <dune/xt/test/main.hxx>

#include <dune/xt/common/vector.hh>
#include <dune/xt/common/fvector.hh>

using VectorTypes = ::testing::Types<std::tuple<Dune::FieldVector<int, 1>, Int<1>>,
                                     std::tuple<Dune::FieldVector<int, 2>, Int<2>>,
                                     std::tuple<Dune::FieldVector<int, 3>, Int<3>>,
                                     std::tuple<Dune::XT::Common::FieldVector<int, 1>, Int<1>>,
                                     std::tuple<Dune::XT::Common::FieldVector<int, 2>, Int<2>>,
                                     std::tuple<Dune::XT::Common::FieldVector<int, 3>, Int<3>>,
                                     std::tuple<Dune::DynamicVector<int>, Int<1>>,
                                     std::tuple<Dune::DynamicVector<int>, Int<2>>,
                                     std::tuple<Dune::DynamicVector<int>, Int<3>>>;


template <class Tuple>
struct SerializeTest : public ::testing::Test
{
  using VectorType = typename std::tuple_element<0, Tuple>::type;
  using V = Dune::XT::Common::VectorAbstraction<VectorType>;
  using ScalarType = typename V::ScalarType;
  static constexpr size_t size = std::tuple_element<1, Tuple>::type::value;

  SerializeTest()
    : vector_(V::create(size))
  {
    for (size_t ii = 0; ii < size; ++ii)
      V::set_entry(vector_, ii, static_cast<ScalarType>(100 * ii));
  }

  template <class T>
  void check_type()
  {
    auto data = Dune::XT::Common::serialize<T>(vector_);
    for (size_t ii = 0; ii < size; ++ii)
      EXPECT_EQ(V::get_entry(vector_, ii), data[ii]);
  }

  void check_default()
  {
    auto data = Dune::XT::Common::serialize(vector_);
    for (size_t ii = 0; ii < size; ++ii)
      EXPECT_EQ(V::get_entry(vector_, ii), data[ii]);
  }

  VectorType vector_;
}; // struct SerializeTest


TYPED_TEST_SUITE(SerializeTest, VectorTypes);
TYPED_TEST(SerializeTest, default)
{
  this->check_default();
}
TYPED_TEST(SerializeTest, int)
{
  this->template check_type<int>();
}
TYPED_TEST(SerializeTest, long_int)
{
  this->template check_type<long int>();
}
TYPED_TEST(SerializeTest, long_long_int)
{
  this->template check_type<long long int>();
}
TYPED_TEST(SerializeTest, ssize_t)
{
  this->template check_type<ssize_t>();
}
TYPED_TEST(SerializeTest, unsigned_int)
{
  this->template check_type<unsigned int>();
}
TYPED_TEST(SerializeTest, unsigned_long_int)
{
  this->template check_type<unsigned long int>();
}
TYPED_TEST(SerializeTest, unsigned_long_long_int)
{
  this->template check_type<unsigned long long int>();
}
TYPED_TEST(SerializeTest, size_t)
{
  this->template check_type<size_t>();
}
TYPED_TEST(SerializeTest, float)
{
  this->template check_type<float>();
}
TYPED_TEST(SerializeTest, double)
{
  this->template check_type<double>();
}
