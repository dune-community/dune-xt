// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018)
//   Tobias Leibner  (2018)

#include <dune/xt/common/test/main.hxx> // <- Needs to come first, include the config.h.

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/float_cmp.hh>

// using this function to allow for {...} syntax below (would not work in the EXPECT_EQ macro directly).
template <class A, int a, class B, int b>
void check(const Dune::XT::Common::FieldVector<A, a>& expected_result, const Dune::FieldVector<B, b>& actual_result)
{
  EXPECT_EQ(expected_result, actual_result);
}

GTEST_TEST(dune_xt_common_field_vector, creation_and_calculations)
{
  using VectorType = Dune::XT::Common::FieldVector<double, 4>;
  using BaseVectorType = Dune::FieldVector<double, 4>;
  using RowMatrixType = Dune::FieldMatrix<double, 1, 4>;
  using ColMatrixType = Dune::FieldMatrix<double, 4, 1>;
  using ArrayType = std::array<double, 4>;
  VectorType zeros;
  for (const auto& val : zeros)
    EXPECT_EQ(val, 0.);
  VectorType test_vec{0, 1, 2, 3};
  for (size_t ii = 0; ii < 4; ++ii)
    EXPECT_EQ(test_vec[ii], ii);
  EXPECT_ANY_THROW(VectorType({1, 2, 3}));
  VectorType copy_assigned = test_vec;
  check(copy_assigned, test_vec);
  VectorType copy_constructed(test_vec);
  check(copy_constructed, test_vec);
  std::vector<double> std_vec{0, 1, 2, 3};
  VectorType constructed_from_std_vec(std_vec);
  check(constructed_from_std_vec, test_vec);
  std::vector<double> std_vec2 = test_vec;
  for (size_t ii = 0; ii < 4; ++ii)
    EXPECT_EQ(std_vec2[ii], ii);
  std_vec.resize(2);
  EXPECT_ANY_THROW(VectorType{std_vec});
  BaseVectorType base_vec{0, 1, 2, 3};
  VectorType constructed_from_base(base_vec);
  check(constructed_from_base, test_vec);
  RowMatrixType row_mat = test_vec;
  for (size_t ii = 0; ii < 4; ++ii)
    EXPECT_EQ(row_mat[0][ii], ii);
  // TODO: reenable this test
  //  ColMatrixType col_mat = test_vec;
  //  for (size_t ii = 0; ii < 4; ++ii)
  //    EXPECT_EQ(col_mat[ii][0], ii);
  ArrayType array = test_vec;
  for (size_t ii = 0; ii < 4; ++ii)
    EXPECT_EQ(array[ii], ii);
  // TODO: remove this when the construction above is reenabled
  ColMatrixType col_mat;
  for (size_t ii = 0; ii < 4; ++ii)
    col_mat[ii][0] = ii;
  EXPECT_DOUBLE_EQ(test_vec * test_vec, 14.);
  EXPECT_DOUBLE_EQ(test_vec * col_mat, 14.);
  for (size_t ii = 0; ii < 4; ++ii)
    EXPECT_DOUBLE_EQ((test_vec * 2.)[ii], 2. * ii);
  for (size_t ii = 0; ii < 4; ++ii)
    EXPECT_DOUBLE_EQ((test_vec / 2.)[ii], ii / 2.);
}


GTEST_TEST(blockedfieldvector, creation_and_calculations)
{
  using BlockedVectorType = Dune::XT::Common::BlockedFieldVector<double, 2, 2>;
  using VectorType = Dune::FieldVector<double, 4>;
  using BlockType = typename BlockedVectorType::BlockType;
  //  using RowMatrixType = Dune::FieldMatrix<double, 1, 4>;
  //  using ColMatrixType = Dune::FieldMatrix<double, 4, 1>;
  //  using ArrayType = std::array<double, 4>;
  BlockedVectorType zeros;
  const BlockedVectorType const_zeros;
  for (const auto& val : zeros)
    EXPECT_EQ(val, 0.);
  for (const auto& val : const_zeros)
    EXPECT_EQ(val, 0.);
  VectorType vec{0, 1, 2, 3};
  BlockedVectorType vec_blocked(vec);
  BlockedVectorType vec_blocked2 = vec;
  const BlockedVectorType const_vec_blocked(vec);
  BlockedVectorType vec_blocked3(BlockType{0., 1.});
  for (size_t ii = 0; ii < 4; ++ii) {
    EXPECT_EQ(vec_blocked[ii], ii);
    EXPECT_EQ(vec_blocked2[ii], ii);
    EXPECT_EQ(const_vec_blocked[ii], ii);
  }
  for (size_t jj = 0; jj < 2; ++jj) {
    for (size_t ii = 0; ii < 2; ++ii) {
      EXPECT_EQ(vec_blocked.block(jj)[ii], ii + 2 * jj);
      EXPECT_EQ(vec_blocked.get_entry(jj, ii), ii + 2 * jj);
      EXPECT_EQ(vec_blocked3.get_entry(jj, ii), ii);
      EXPECT_EQ(const_vec_blocked.block(jj)[ii], ii + 2 * jj);
      EXPECT_EQ(const_vec_blocked.get_entry(jj, ii), ii + 2 * jj);
    }
  }
  EXPECT_DOUBLE_EQ(vec_blocked.one_norm(), 6);
  EXPECT_DOUBLE_EQ(vec_blocked.two_norm(), std::sqrt(14));
  EXPECT_DOUBLE_EQ(vec_blocked.two_norm2(), 14);
  EXPECT_DOUBLE_EQ(vec_blocked * vec_blocked, 14);
  VectorType assigned_vec = vec_blocked;
  for (size_t ii = 0; ii < 4; ++ii)
    EXPECT_EQ(assigned_vec[ii], vec_blocked[ii]);
  EXPECT_EQ(*vec_blocked.data(), 0.);
  EXPECT_TRUE(Dune::XT::Common::FloatCmp::eq(vec_blocked - vec_blocked, zeros));
  EXPECT_TRUE(Dune::XT::Common::FloatCmp::eq(vec_blocked + (vec_blocked * (-1)), zeros));
  EXPECT_TRUE(Dune::XT::Common::FloatCmp::eq(vec_blocked * 0., zeros));
  auto vec_blocked_copy = vec_blocked;
  vec_blocked_copy *= 0.;
  EXPECT_TRUE(Dune::XT::Common::FloatCmp::eq(vec_blocked_copy, zeros));
  vec_blocked_copy = vec_blocked;
  vec_blocked_copy += vec_blocked;
  for (size_t ii = 0; ii < 4; ++ii)
    EXPECT_DOUBLE_EQ(vec_blocked_copy[ii], 2. * ii);
  vec_blocked_copy -= vec_blocked;
  EXPECT_TRUE(Dune::XT::Common::FloatCmp::eq(vec_blocked_copy, vec_blocked));
}


GTEST_TEST(hstack, dune_field_vector)
{
  const Dune::FieldVector<double, 1> one(1);
  const Dune::FieldVector<double, 2> two(2);
  // hstack is specialized for one and two arguments, so testing up to three should be enough
  using Dune::XT::Common::hstack;
  check(Dune::XT::Common::FieldVector<double, 2>({1, 3}), hstack(one, 3));
  check(Dune::XT::Common::FieldVector<double, 2>({3, 1}), hstack(3, one));
  check(Dune::XT::Common::FieldVector<double, 3>({1, 2, 2}), hstack(one, two));
  check(Dune::XT::Common::FieldVector<double, 3>({2, 2, 1}), hstack(two, one));
  check(Dune::XT::Common::FieldVector<double, 4>({3, 1, 2, 2}), hstack(3, one, two));
  check(Dune::XT::Common::FieldVector<double, 4>({1, 3, 2, 2}), hstack(one, 3, two));
  check(Dune::XT::Common::FieldVector<double, 4>({1, 2, 2, 3}), hstack(one, two, 3));
}


GTEST_TEST(hstack, dune_xt_common_field_vector)
{
  using namespace Dune;
  const Dune::XT::Common::FieldVector<double, 1> one(1);
  const Dune::XT::Common::FieldVector<double, 2> two(2);
  // hstack is specialized for one and two arguments, so testing up to three should be enough
  using Dune::XT::Common::hstack;
  check(Dune::XT::Common::FieldVector<double, 2>({1, 3}), hstack(one, 3));
  check(Dune::XT::Common::FieldVector<double, 2>({3, 1}), hstack(3, one));
  check(Dune::XT::Common::FieldVector<double, 3>({1, 2, 2}), hstack(one, two));
  check(Dune::XT::Common::FieldVector<double, 3>({2, 2, 1}), hstack(two, one));
  check(Dune::XT::Common::FieldVector<double, 4>({3, 1, 2, 2}), hstack(3, one, two));
  check(Dune::XT::Common::FieldVector<double, 4>({1, 3, 2, 2}), hstack(one, 3, two));
  check(Dune::XT::Common::FieldVector<double, 4>({1, 2, 2, 3}), hstack(one, two, 3));
}


GTEST_TEST(hstack, both_field_vectors)
{
  using namespace Dune;
  const Dune::XT::Common::FieldVector<double, 1> one(1);
  const Dune::FieldVector<double, 2> two(2);
  // hstack is specialized for one and two arguments, so testing up to three should be enough
  using Dune::XT::Common::hstack;
  check(Dune::XT::Common::FieldVector<double, 3>({1, 2, 2}), hstack(one, two));
  check(Dune::XT::Common::FieldVector<double, 3>({2, 2, 1}), hstack(two, one));
  check(Dune::XT::Common::FieldVector<double, 4>({3, 1, 2, 2}), hstack(3, one, two));
  check(Dune::XT::Common::FieldVector<double, 4>({1, 3, 2, 2}), hstack(one, 3, two));
  check(Dune::XT::Common::FieldVector<double, 4>({1, 2, 2, 3}), hstack(one, two, 3));
}
