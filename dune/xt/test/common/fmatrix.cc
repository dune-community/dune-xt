// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   René Fritze    (2018)
//   Tobias Leibner (2018)

#include <dune/xt/common/test/main.hxx> // <- Needs to come first, include the config.h.

#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/float_cmp.hh>

GTEST_TEST(dune_xt_common_field_matrix, creation_and_calculations)
{
  static constexpr size_t rows = 2;
  static constexpr size_t cols = 3;
  using MatrixType = Dune::XT::Common::FieldMatrix<double, rows, cols>;
  using SquareMatrixType = Dune::XT::Common::FieldMatrix<double, rows, rows>;
  using BaseMatrixType = Dune::FieldMatrix<double, rows, cols>;
  using ScalarMatrixType = Dune::XT::Common::FieldMatrix<double, 1, 1>;
  using ScalarBaseMatrixType = Dune::FieldMatrix<double, 1, 1>;
  using RowVectorType = Dune::XT::Common::FieldVector<double, cols>;
  using ColVectorType = Dune::XT::Common::FieldVector<double, rows>;
  MatrixType zeros;
  ScalarMatrixType zero;
  MatrixType ones(1.);
  MatrixType ones2(rows, cols, 1.);
  ScalarMatrixType one(1.);
  ScalarMatrixType one2{1.};
  MatrixType mat({{0., 1., 2.}, {3., 4., 5.}});
  EXPECT_ANY_THROW(MatrixType({{0., 2.}, {3., 4., 5.}}));
  EXPECT_ANY_THROW(MatrixType({{0., 1., 2.}, {4., 5.}}));
  EXPECT_ANY_THROW(MatrixType({{0., 2.}, {4., 5.}}));
  EXPECT_ANY_THROW(ScalarMatrixType({{0., 2.}, {4., 5.}}));
#ifdef NDEBUG
  EXPECT_ANY_THROW(MatrixType(rows + 1, cols - 1, 1.));
  EXPECT_ANY_THROW(ScalarMatrixType(2, 3, 1.));
#endif
  EXPECT_EQ(zero[0][0], 0.);
  EXPECT_EQ(one[0][0], 1.);
  EXPECT_EQ(one2[0][0], 1.);
  EXPECT_EQ(one, 1.);
  EXPECT_EQ(ones, ones2);
  for (size_t ii = 0; ii < rows; ++ii) {
    for (size_t jj = 0; jj < cols; ++jj) {
      EXPECT_DOUBLE_EQ(mat[ii][jj], ii * cols + jj);
      EXPECT_EQ(zeros[ii][jj], 0.);
      EXPECT_EQ(ones[ii][jj], 1.);
    }
  }
  MatrixType copy_assigned = mat;
  MatrixType copy_constructed(mat);
  ScalarMatrixType scalar_copy_assigned = one;
  ScalarMatrixType scalar_copy_constructed(one);
  BaseMatrixType base_mat = mat;
  MatrixType base_constructed(base_mat);
  ScalarBaseMatrixType scalar_base_mat = one;
  ScalarMatrixType scalar_base_constructed(scalar_base_mat);
  EXPECT_EQ(copy_assigned, mat);
  EXPECT_EQ(copy_constructed, mat);
  EXPECT_EQ(base_constructed, mat);
  EXPECT_EQ(scalar_copy_assigned, one);
  EXPECT_EQ(scalar_copy_constructed, one);
  EXPECT_EQ(scalar_base_constructed, one);
  auto transposed_mat = mat.transpose();
  auto scalar_transposed_mat = one.transpose();
  EXPECT_EQ(scalar_transposed_mat, one);
  for (size_t ii = 0; ii < rows; ++ii)
    for (size_t jj = 0; jj < cols; ++jj)
      EXPECT_EQ(transposed_mat[jj][ii], mat[ii][jj]);
  RowVectorType ones_vec(1.);
  ColVectorType mat_times_vec2;
  auto mat_times_vec = mat * ones_vec;
  mat.mv(ones_vec, mat_times_vec2);
  for (size_t ii = 0; ii < rows; ++ii) {
    EXPECT_DOUBLE_EQ(mat_times_vec[ii], std::accumulate(mat[ii].begin(), mat[ii].end(), 0.));
    EXPECT_DOUBLE_EQ(mat_times_vec[ii], mat_times_vec2[ii]);
  }
  auto mat_scaled = mat;
  mat_scaled *= 2.;
  auto mat_scaled2 = mat * 2.;
  for (size_t ii = 0; ii < rows; ++ii)
    for (size_t jj = 0; jj < cols; ++jj)
      EXPECT_DOUBLE_EQ(mat_scaled[ii][jj], mat_scaled2[ii][jj]);
  SquareMatrixType square_mat{{0, 1}, {2, 3}};
  SquareMatrixType expected_inv{{-1.5, 0.5}, {1, 0}};
  auto inv_square_mat = square_mat;
  inv_square_mat.invert();
  for (size_t ii = 0; ii < rows; ++ii)
    for (size_t jj = 0; jj < rows; ++jj)
      EXPECT_DOUBLE_EQ(expected_inv[ii][jj], inv_square_mat[ii][jj]);
  ScalarMatrixType fvec_constructed(Dune::XT::Common::FieldVector<double, 1>{1.});
  ScalarMatrixType fvec_constructed2(Dune::FieldVector<double, 1>{1.});
  ScalarMatrixType fvec_assigned = Dune::FieldVector<double, 1>{1.};
  EXPECT_EQ(fvec_constructed, one);
  EXPECT_EQ(fvec_constructed2, one);
  EXPECT_EQ(fvec_assigned, one);
  auto two = one * 2.;
  EXPECT_EQ(two[0][0], 2.);


#if 0
template <class L, int ROWS, int COLS, class R>
Dune::XT::Common::FieldMatrix<typename PromotionTraits<L, R>::PromotedType, ROWS, COLS>
operator-(const Dune::FieldMatrix<L, ROWS, COLS>& left, const Dune::FieldMatrix<R, ROWS, COLS>& right)
{
  Dune::XT::Common::FieldMatrix<typename PromotionTraits<L, R>::PromotedType, ROWS, COLS> ret = left;
  ret -= right;
  return ret;
}


template <class L, int ROWS, int COLS, class R>
Dune::XT::Common::FieldMatrix<typename PromotionTraits<L, R>::PromotedType, ROWS, COLS>
operator+(const Dune::FieldMatrix<L, ROWS, COLS>& left, const Dune::FieldMatrix<R, ROWS, COLS>& right)
{
  Dune::XT::Common::FieldMatrix<typename PromotionTraits<L, R>::PromotedType, ROWS, COLS> ret = left;
  ret += right;
  return ret;
}


template <class K, int L_ROWS, int L_COLS, int R_COLS>
Dune::XT::Common::FieldMatrix<K, L_ROWS, R_COLS> operator*(const Dune::FieldMatrix<K, L_ROWS, L_COLS>& left,
                                                           const Dune::FieldMatrix<K, L_COLS, R_COLS>& right)
{
  return left.rightmultiplyany(right);
}

// we need this explicit overload to fix an ambiguous operator* error due to the automatic conversion from
// FieldMatrix<K, 1, 1> to K
template <class K, int L_ROWS>
Dune::XT::Common::FieldMatrix<K, L_ROWS, 1> operator*(const Dune::XT::Common::FieldMatrix<K, L_ROWS, 1>& left,
                                                      const Dune::FieldMatrix<K, 1, 1>& right)
{
  return left.rightmultiplyany(right);
}


template <class K, int L_ROWS, int L_COLS, int R_COLS>
void rightmultiply(Dune::FieldMatrix<K, L_ROWS, R_COLS>& ret,
                   const Dune::FieldMatrix<K, L_ROWS, L_COLS>& left,
                   const Dune::FieldMatrix<K, L_COLS, R_COLS>& right)
{
  for (size_t ii = 0; ii < L_ROWS; ++ii) {
    for (size_t jj = 0; jj < R_COLS; ++jj) {
      ret[ii][jj] = 0.;
      for (size_t kk = 0; kk < L_COLS; ++kk)
        ret[ii][jj] += left[ii][kk] * right[kk][jj];
    }
  }
}

template <class L, int L_ROWS, int L_COLS, class R, int R_COLS>
typename std::enable_if<!std::is_same<L, R>::value,
                        Dune::XT::Common::FieldMatrix<typename PromotionTraits<L, R>::PromotedType, L_ROWS, R_COLS>>::
    type
    operator*(const Dune::FieldMatrix<L, L_ROWS, L_COLS>& left, const Dune::FieldMatrix<R, L_COLS, R_COLS>& right)
{
  using Promoted = Dune::XT::Common::FieldMatrix<typename PromotionTraits<L, R>::PromotedType, L_ROWS, R_COLS>;
  using Dune::XT::Common::convert_to;
  return convert_to<Promoted>(left).rightmultiplyany(convert_to<Promoted>(right));
}

// versions that do not allocate matrices on the stack (for large matrices)
template <class K, int L_ROWS, int L_COLS, int R_COLS>
std::unique_ptr<Dune::XT::Common::FieldMatrix<K, L_ROWS, R_COLS>>
operator*(const std::unique_ptr<Dune::FieldMatrix<K, L_ROWS, L_COLS>>& left,
          const Dune::FieldMatrix<K, L_COLS, R_COLS>& right)
{
  auto ret = std::make_unique<Dune::XT::Common::FieldMatrix<K, L_ROWS, R_COLS>>();
  rightmultiply(*ret, *left, right);
  return ret;
}

template <class K, int L_ROWS, int L_COLS, int R_COLS>
std::unique_ptr<Dune::XT::Common::FieldMatrix<K, L_ROWS, R_COLS>>
operator*(const Dune::FieldMatrix<K, L_ROWS, L_COLS>& left,
          const std::unique_ptr<Dune::FieldMatrix<K, L_COLS, R_COLS>>& right)
{
  auto ret = std::make_unique<Dune::XT::Common::FieldMatrix<K, L_ROWS, R_COLS>>();
  rightmultiply(*ret, left, *right);
  return ret;
}

template <class K, int L_ROWS, int L_COLS, int R_COLS>
std::unique_ptr<Dune::XT::Common::FieldMatrix<K, L_ROWS, R_COLS>>
operator*(const std::unique_ptr<Dune::FieldMatrix<K, L_ROWS, L_COLS>>& left,
          const std::unique_ptr<Dune::FieldMatrix<K, L_COLS, R_COLS>>& right)
{
  return left * *right;
}
#endif
}

GTEST_TEST(blockedfieldmatrix, creation_and_calculations)
{
  static constexpr size_t num_blocks = 2;
  static constexpr size_t block_rows = 2;
  static constexpr size_t block_cols = 3;
  static constexpr size_t rows = num_blocks * block_rows;
  static constexpr size_t cols = num_blocks * block_cols;
  using BlockedMatrixType = Dune::XT::Common::BlockedFieldMatrix<double, num_blocks, block_rows, block_cols>;
  using MatrixType = Dune::FieldMatrix<double, rows, cols>;
  using XtMatrixType = Dune::XT::Common::FieldMatrix<double, rows, cols>;
  using RowVectorType = Dune::XT::Common::FieldVector<double, cols>;
  using ColVectorType = Dune::XT::Common::FieldVector<double, rows>;
  using BlockedRowVectorType = Dune::XT::Common::BlockedFieldVector<double, num_blocks, block_cols>;
  using BlockedColVectorType = Dune::XT::Common::BlockedFieldVector<double, num_blocks, block_rows>;
  using BlockType = BlockedMatrixType::BlockType;
  BlockedMatrixType zeros;
  BlockedMatrixType ones(1.);
  MatrixType mat({{0, 1, 2, 0, 0, 0}, {3, 4, 5, 0, 0, 0}, {0, 0, 0, 0, 1, 2}, {0, 0, 0, 3, 4, 5}});
  BlockedMatrixType block_mat(mat);
  BlockedMatrixType block_mat2 = mat;
  BlockType block{{0, 1, 2}, {3, 4, 5}};
  BlockedMatrixType block_mat3(block);
  const BlockedMatrixType const_block_mat = block_mat;
  EXPECT_EQ(block_mat, block_mat2);
  EXPECT_EQ(block_mat, block_mat3);
  for (size_t jj = 0; jj < num_blocks; ++jj) {
    EXPECT_EQ(block_mat.block(jj), block);
    EXPECT_EQ(const_block_mat.block(jj), block);
  }
  for (size_t ii = 0; ii < rows; ++ii) {
    for (size_t jj = 0; jj < cols; ++jj) {
      EXPECT_EQ(zeros.get_entry(ii, jj), 0.);
      EXPECT_EQ(block_mat.get_entry(ii, jj), mat[ii][jj]);
      EXPECT_EQ(const_block_mat.get_entry(ii, jj), mat[ii][jj]);
    }
  }
  for (size_t jj = 0; jj < num_blocks; ++jj) {
    for (size_t ll = 0; ll < block_rows; ++ll) {
      for (size_t mm = 0; mm < block_cols; ++mm) {
        EXPECT_EQ(ones.block(jj)[ll][mm], 1.);
        EXPECT_EQ(block_mat.get_entry(jj, ll, mm), mat[jj * block_rows + ll][jj * block_cols + mm]);
        EXPECT_EQ(const_block_mat.get_entry(jj, ll, mm), mat[jj * block_rows + ll][jj * block_cols + mm]);
        EXPECT_EQ(block_mat.block(jj)[ll][mm], mat[jj * block_rows + ll][jj * block_cols + mm]);
        EXPECT_EQ(const_block_mat.block(jj)[ll][mm], mat[jj * block_rows + ll][jj * block_cols + mm]);
      }
    }
  }
  RowVectorType ones_vec(1.);
  ColVectorType mat_times_vec;
  BlockedRowVectorType block_ones_vec(1.);
  BlockedColVectorType block_mat_times_vec;
  block_mat.mv(ones_vec, mat_times_vec);
  block_mat.mv(block_ones_vec, block_mat_times_vec);
  for (size_t ii = 0; ii < rows; ++ii) {
    EXPECT_DOUBLE_EQ(mat_times_vec[ii], std::accumulate(mat[ii].begin(), mat[ii].end(), 0.));
    EXPECT_DOUBLE_EQ(mat_times_vec[ii], block_mat_times_vec[ii]);
  }
  ColVectorType ones_vec2(1.);
  RowVectorType mtv_vec;
  BlockedColVectorType block_ones_vec2(1.);
  BlockedRowVectorType block_mtv_vec;
  block_mat.mtv(ones_vec2, mtv_vec);
  block_mat.mtv(block_ones_vec2, block_mtv_vec);
  auto transposed_mat = XtMatrixType(mat).transpose();
  auto transposed_block_mat = block_mat.transpose();
  for (size_t ii = 0; ii < cols; ++ii) {
    EXPECT_DOUBLE_EQ(mtv_vec[ii], std::accumulate(transposed_mat[ii].begin(), transposed_mat[ii].end(), 0.));
    EXPECT_DOUBLE_EQ(mtv_vec[ii], block_mtv_vec[ii]);
  }
  for (size_t ii = 0; ii < rows; ++ii)
    for (size_t jj = 0; jj < cols; ++jj)
      EXPECT_EQ(transposed_block_mat.get_entry(jj, ii), transposed_mat[jj][ii]);
  Dune::XT::Common::BlockedFieldMatrix<double, 2, 2, 2> block_square_mat;
  size_t val = 0;
  for (size_t jj = 0; jj < 2; ++jj)
    for (size_t ll = 0; ll < 2; ++ll)
      for (size_t mm = 0; mm < 2; ++mm)
        block_square_mat.block(jj)[ll][mm] = val++;
  auto block_square_mat_copy = block_square_mat;
  block_square_mat.rightmultiply(block_square_mat_copy);
  Dune::FieldMatrix<double, 4, 4> expected_mat{{2, 3, 0, 0}, {6, 11, 0, 0}, {0, 0, 46, 55}, {0, 0, 66, 79}};
  for (size_t ii = 0; ii < 4; ++ii)
    for (size_t jj = 0; jj < 4; ++jj)
      EXPECT_EQ(block_square_mat.get_entry(ii, jj), expected_mat[ii][jj]);
}
