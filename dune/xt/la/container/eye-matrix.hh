// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_XT_LA_CONTAINER_EYE_MATRIX_HH
#define DUNE_XT_LA_CONTAINER_EYE_MATRIX_HH

#include <dune/xt/common/matrix.hh>

#include <dune/xt/la/type_traits.hh>

#include "pattern.hh"

namespace Dune {
namespace XT {
namespace LA {
namespace internal {


template <class MatrixType>
typename std::enable_if<Common::is_matrix<MatrixType>::value, void>::type set_diagonal_to_one(MatrixType& mat)
{
  using M = Common::MatrixAbstraction<MatrixType>;
  for (size_t ii = 0; ii < std::min(M::rows(mat), M::cols(mat)); ++ii)
    M::set_entry(mat, ii, ii, 1.);
}


} // namespace internal


template <class MatrixType>
typename std::enable_if<is_matrix<MatrixType>::value, MatrixType>::type
eye_matrix(const size_t rows, const size_t cols, const SparsityPatternDefault& pattern = SparsityPatternDefault())
{
  MatrixType mat = MatrixType(rows, cols, pattern.size() == 0 ? diagonal_pattern(rows, cols) : pattern);
  internal::set_diagonal_to_one(mat);
  return mat;
}

template <class MatrixType>
typename std::enable_if<Common::is_matrix<MatrixType>::value, void>::type eye_matrix(MatrixType& matrix)
{
  matrix *= 0.;
  internal::set_diagonal_to_one(matrix);
}

template <class MatrixType>
typename std::enable_if<Common::is_matrix<MatrixType>::value && !is_matrix<MatrixType>::value, MatrixType>::type
eye_matrix(const size_t rows, const size_t cols, const SparsityPatternDefault& /*pattern*/ = SparsityPatternDefault())
{
  using M = Common::MatrixAbstraction<MatrixType>;
  auto mat = M::create(rows, cols, typename M::ScalarType(0.));
  for (size_t ii = 0; ii < std::min(rows, cols); ++ii)
    M::set_entry(mat, ii, ii, 1);
  return mat;
}


template <class MatrixType>
typename std::enable_if<Common::is_matrix<MatrixType>::value, MatrixType>::type
eye_matrix(const size_t size, const SparsityPatternDefault& pattern = SparsityPatternDefault())
{
  return eye_matrix<MatrixType>(size, size, pattern);
}


} // namespace LA
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_LA_CONTAINER_EYE_MATRIX_HH
