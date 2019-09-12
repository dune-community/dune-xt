// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner       (2019)

#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_DEBUG_LOGGING 1
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING 1

#include <dune/xt/common/test/main.hxx> // <- This one has to come first, includes config.h!
#include <dune/xt/common/test/gtest/gtest.h>

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/la/container/pattern.hh>

GTEST_TEST(SparsityPatternDefaultTest, test_interface)
{
  using namespace Dune;
  using Pattern = XT::LA::SparsityPatternDefault;
  constexpr size_t ROWS = 5, COLS = 4;

  Pattern pattern(ROWS);
  const Pattern const_pattern(ROWS);
  EXPECT_EQ(pattern.size(), ROWS);
  for (size_t ii = 0; ii < ROWS; ++ii)
    pattern.insert(ii, std::min(ii, COLS - 1));
  for (size_t ii = 0; ii < ROWS; ++ii) {
    EXPECT_EQ(pattern.inner(ii).size(), 1);
    EXPECT_EQ(const_pattern.inner(ii).size(), 0);
  }
  {
    size_t ii = 0;
    for (const auto& inner : pattern) {
      EXPECT_TRUE(pattern.inner(ii) == inner);
      ++ii;
    }
    ii = 0;
    auto pattern_copy = pattern;
    for (auto& inner : pattern_copy) {
      EXPECT_TRUE(pattern_copy.inner(ii) == inner);
      pattern_copy.inner(ii).push_back(0);
      ++ii;
    }
  }
  auto pattern2 = pattern;
  EXPECT_TRUE(pattern2.contains(pattern));
  EXPECT_TRUE(pattern.contains(pattern2));
  EXPECT_TRUE(pattern.contains(const_pattern));
  pattern.sort();
  for (size_t ii = 0; ii < ROWS; ++ii)
    pattern2.sort(ii);
  EXPECT_TRUE(pattern == pattern2);
  EXPECT_TRUE(pattern != const_pattern);
  auto test_pattern = pattern + pattern2;
  EXPECT_TRUE(test_pattern.contains(pattern));
  EXPECT_TRUE(test_pattern.contains(pattern2));
  EXPECT_TRUE(pattern + pattern2 == pattern);
  Pattern pattern3(ROWS), pattern4(ROWS), pattern5(ROWS);
  for (size_t ii = 0; ii < ROWS; ++ii) {
    pattern3.insert(ii, 0);
    pattern4.insert(ii, COLS - 1);
    pattern5.insert(ii, 0);
    pattern5.insert(ii, COLS - 1);
  }
  const auto pattern5_T = pattern5.transposed(COLS);
  EXPECT_TRUE(pattern5_T.size() == COLS);
  for (size_t ii = 0; ii < ROWS; ++ii) {
    EXPECT_TRUE(pattern3.contains(ii, 0));
    EXPECT_FALSE(pattern3.contains(ii, 1));
    EXPECT_TRUE(pattern5_T.contains(0, ii));
    EXPECT_TRUE(pattern5_T.contains(COLS - 1, ii));
    EXPECT_FALSE(pattern5_T.contains(1, ii));
    EXPECT_FALSE(pattern5_T.contains(2, ii));
  }
  EXPECT_FALSE(pattern3 == pattern4);
  EXPECT_TRUE(pattern3 + pattern4 == pattern5);
}

GTEST_TEST(SparsityPatternDefaultTest, test_creation_functions)
{
  using namespace Dune;
  constexpr size_t ROWS = 5, COLS = 4;
  auto dense_patt = XT::LA::dense_pattern(ROWS, COLS);
  EXPECT_TRUE(dense_patt.size() == ROWS);
  for (size_t ii = 0; ii < ROWS; ++ii) {
    EXPECT_TRUE(dense_patt.inner(ii).size() == COLS);
    for (size_t jj = 0; jj < COLS; ++jj)
      EXPECT_TRUE(dense_patt.contains(ii, jj));
  }
  auto lower_triangular_patt = XT::LA::triangular_pattern(ROWS, COLS);
  auto lower_triangular_patt2 = XT::LA::triangular_pattern(ROWS, COLS, XT::Common::MatrixPattern::lower_triangular);
  EXPECT_TRUE(lower_triangular_patt == lower_triangular_patt2);
  auto upper_triangular_patt = XT::LA::triangular_pattern(ROWS, COLS, XT::Common::MatrixPattern::upper_triangular);
  EXPECT_TRUE(lower_triangular_patt.size() == ROWS);
  EXPECT_TRUE(upper_triangular_patt.size() == ROWS);
  for (size_t ii = 0; ii < ROWS; ++ii) {
    if (ii < COLS) {
      EXPECT_TRUE(lower_triangular_patt.contains(ii, ii));
      EXPECT_TRUE(upper_triangular_patt.contains(ii, ii));
    }
    for (size_t jj = 0; jj < ii; ++jj) {
      EXPECT_TRUE(lower_triangular_patt.contains(ii, jj));
      EXPECT_FALSE(upper_triangular_patt.contains(ii, jj));
    }
    for (size_t jj = 0; jj > ii; ++jj) {
      EXPECT_FALSE(lower_triangular_patt.contains(ii, jj));
      EXPECT_TRUE(upper_triangular_patt.contains(ii, jj));
    }
  } // ii
  auto diagonal_patt = XT::LA::diagonal_pattern(ROWS, COLS);
  auto subdiagonal_patt = XT::LA::diagonal_pattern(ROWS, COLS, -1);
  EXPECT_TRUE(diagonal_patt.size() == ROWS);
  EXPECT_TRUE(subdiagonal_patt.size() == ROWS);
  for (size_t ii = 0; ii < ROWS; ++ii) {
    if (ii < COLS) {
      EXPECT_TRUE(diagonal_patt.contains(ii, ii));
    }
    if (ii > 0) {
      EXPECT_TRUE(subdiagonal_patt.contains(ii, ii - 1));
    }
    for (size_t jj = 0; jj < COLS; ++jj) {
      if (jj != ii) {
        EXPECT_FALSE(diagonal_patt.contains(ii, jj));
      }
      if (jj != ii - 1) {
        EXPECT_FALSE(subdiagonal_patt.contains(ii, jj));
      }
    } // jj
  } // ii
  const auto diagonal_subdiagonal_patt = diagonal_patt + subdiagonal_patt;
  EXPECT_TRUE(diagonal_subdiagonal_patt.contains(diagonal_patt));
  EXPECT_TRUE(diagonal_subdiagonal_patt.contains(subdiagonal_patt));
  EXPECT_FALSE(diagonal_patt.contains(diagonal_subdiagonal_patt));
  EXPECT_FALSE(subdiagonal_patt.contains(diagonal_subdiagonal_patt));

  XT::LA::SparsityPatternDefault pattern(ROWS), pattern2(COLS);
  pattern.insert(0, 0);
  pattern.insert(0, COLS - 1);
  pattern.insert(ROWS - 1, 0);
  pattern.insert(ROWS - 1, COLS - 1);
  for (size_t ii = 0; ii < COLS; ++ii) {
    pattern2.insert(ii, 0);
    pattern2.insert(ii, ROWS - 1);
  }
  auto mult_pattern = XT::LA::multiplication_pattern(pattern, pattern2, ROWS);
  EXPECT_TRUE(mult_pattern.size() == ROWS);
  for (size_t ii = 0; ii < ROWS; ++ii) {
    for (size_t jj = 0; jj < ROWS; ++jj) {
      if ((ii == 0 || ii == ROWS - 1) && (jj == 0 || jj == ROWS - 1))
        EXPECT_TRUE(mult_pattern.contains(ii, jj));
      else
        EXPECT_FALSE(mult_pattern.contains(ii, jj));
    }
  }
}
