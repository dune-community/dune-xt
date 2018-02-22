// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2012 - 2017)
//   Rene Milk       (2010 - 2016, 2018)
//   Sven Kaulmann   (2010 - 2011)
//   Tobias Leibner  (2014, 2017)

#ifndef DUNE_XT_COMMON_MATRIX_HH
#define DUNE_XT_COMMON_MATRIX_HH

#include <memory>

#include <dune/common/dynmatrix.hh>
#include <dune/common/fmatrix.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/common/pattern.hh>

namespace Dune {
namespace XT {
namespace Common {


template <class FirstMatrixType, class SecondMatrixType, class ReturnMatrixType>
void multiply_helper(const FirstMatrixType& first, const SecondMatrixType& second, ReturnMatrixType& ret);

/**
 * \brief Traits to uniformly handle dense matrices.
 *
 *        If you want your matrix class to benefit from the functionality defined in this header you have to manually
 *        specify a specialization of this class in your code with is_matrix defined to true and the appropriate
 *        static methods implemented and members defined (see the specializations below).
 */
template <class MatType>
struct MatrixAbstraction
{
  typedef MatType MatrixType;
  typedef MatType ScalarType;
  typedef MatType RealType;
  typedef MatType S;
  typedef MatType R;

  template <size_t rows = 0, size_t cols = 0, class Field = ScalarType>
  using MatrixTypeTemplate = MatrixType;

  static const bool is_matrix = false;

  static const bool has_static_size = false;

  static const size_t static_rows = std::numeric_limits<size_t>::max();

  static const size_t static_cols = std::numeric_limits<size_t>::max();

  static const constexpr StorageLayout storage_layout = XT::Common::StorageLayout::other;

  template <size_t ROWS = static_rows, size_t COLS = static_cols, class Field = ScalarType>
  static inline /*MatrixType*/ void create(const size_t /*rows*/,
                                           const size_t /*cols*/,
                                           const SparsityPatternDefault& /*pattern*/ = SparsityPatternDefault())
  {
    static_assert(AlwaysFalse<MatType>::value, "Do not call me if is_matrix is false!");
  }

  template <size_t ROWS = static_rows, size_t COLS = static_cols, class Field = ScalarType>
  static inline /*MatrixType*/ void create(const size_t /*rows*/,
                                           const size_t /*cols*/,
                                           const ScalarType& /*val*/,
                                           const SparsityPatternDefault& /*pattern*/ = SparsityPatternDefault())
  {
    static_assert(AlwaysFalse<MatType>::value, "Do not call me if is_matrix is false!");
  }

  static inline /*std::unique_ptr<MatrixType>*/ void
  create_unique_ptr(const size_t /*rows*/, const size_t /*cols*/, const ScalarType& /*val*/)
  {
    static_assert(AlwaysFalse<MatType>::value, "Do not call me if is_matrix is false!");
  }

  static inline /*size_t*/ void rows(const MatrixType& /*mat*/)
  {
    static_assert(AlwaysFalse<MatType>::value, "Do not call me if is_matrix is false!");
  }

  static inline /*size_t*/ void cols(const MatrixType& /*mat*/)
  {
    static_assert(AlwaysFalse<MatType>::value, "Do not call me if is_matrix is false!");
  }

  static inline void
  set_entry(MatrixType& /*mat*/, const size_t /*row*/, const size_t /*col*/, const ScalarType& /*val*/)
  {
    static_assert(AlwaysFalse<MatType>::value, "Do not call me if is_matrix is false!");
  }

  static inline /*ScalarType*/ void get_entry(const MatrixType& /*mat*/, const size_t /*row*/, const size_t /*col*/)
  {
    static_assert(AlwaysFalse<MatType>::value, "Do not call me if is_matrix is false!");
  }

  static inline void
  add_to_entry(MatrixType& /*mat*/, const size_t /*row*/, const size_t /*col*/, const ScalarType& /*val*/)
  {
    static_assert(AlwaysFalse<MatType>::value, "Do not call me if is_matrix is false!");
  }

  static inline ScalarType* data(MatrixType& /*mat*/)
  {
    static_assert(AlwaysFalse<MatType>::value, "Do not call me if storage_layout is not dense!");
  }

  static inline const ScalarType* data(const MatrixType& /*mat*/)
  {
    static_assert(AlwaysFalse<MatType>::value, "Do not call me if storage_layout is not dense!");
  }
};


//! logically and structurally this belongs in type_utils.hh, but the dependent implementation prohibits that
template <class MatrixType>
struct is_matrix
{
  static const bool value = MatrixAbstraction<MatrixType>::is_matrix;
};


template <class K>
struct MatrixAbstraction<Dune::DynamicMatrix<K>>
{
  typedef Dune::DynamicMatrix<K> MatrixType;
  typedef typename Dune::FieldTraits<K>::field_type ScalarType;
  typedef typename Dune::FieldTraits<K>::real_type RealType;
  typedef ScalarType S;
  typedef RealType R;
  template <size_t rows = 0, size_t cols = 0, class Field = K>
  using MatrixTypeTemplate = DynamicMatrix<K>;

  static const bool is_matrix = true;

  static const bool has_static_size = false;

  static const size_t static_rows = std::numeric_limits<size_t>::max();

  static const size_t static_cols = std::numeric_limits<size_t>::max();

  static const constexpr StorageLayout storage_layout = StorageLayout::other;

  template <size_t ROWS = static_rows, size_t COLS = static_cols, class Field = ScalarType>
  static inline MatrixTypeTemplate<ROWS, COLS, Field>
  create(const size_t rows, const size_t cols, const SparsityPatternDefault& /*pattern*/ = SparsityPatternDefault())
  {
    return MatrixType(rows, cols);
  }

  template <size_t ROWS = static_rows, size_t COLS = static_cols, class Field = ScalarType>
  static inline MatrixTypeTemplate<ROWS, COLS, Field>
  create(const size_t rows,
         const size_t cols,
         const ScalarType& val,
         const SparsityPatternDefault& /*pattern*/ = SparsityPatternDefault())
  {
    return MatrixType(rows, cols, val);
  }

  static inline std::unique_ptr<MatrixType> create_dynamic(const size_t rows, const size_t cols)
  {
    return std::make_unique<MatrixType>(rows, cols);
  }

  static inline std::unique_ptr<MatrixType> create_dynamic(const size_t rows, const size_t cols, const ScalarType& val)
  {
    return std::make_unique<MatrixType>(rows, cols, val);
  }

  static inline size_t rows(const MatrixType& mat)
  {
    return mat.rows();
  }

  static inline size_t cols(const MatrixType& mat)
  {
    return mat.cols();
  }

  static inline void set_entry(MatrixType& mat, const size_t row, const size_t col, const ScalarType& val)
  {
    mat[row][col] = val;
  }

  static inline ScalarType get_entry(const MatrixType& mat, const size_t row, const size_t col)
  {
    return mat[row][col];
  }

  static inline void add_to_entry(MatrixType& mat, const size_t row, const size_t col, const ScalarType& val)
  {
    mat[row][col] += val;
  }

  static inline ScalarType* data(MatrixType& mat)
  {
    return &(mat[0][0]);
  }

  static inline const ScalarType* data(const MatrixType& mat)
  {
    return &(mat[0][0]);
  }
};

template <class K, int N, int M>
struct MatrixAbstraction<Dune::FieldMatrix<K, N, M>>
{
  typedef Dune::FieldMatrix<K, N, M> MatrixType;
  typedef typename Dune::FieldTraits<K>::field_type ScalarType;
  typedef typename Dune::FieldTraits<K>::real_type RealType;
  typedef ScalarType S;
  typedef RealType R;
  template <size_t rows = N, size_t cols = M, class Field = K>
  using MatrixTypeTemplate = Dune::FieldMatrix<Field, rows, cols>;

  static const bool is_matrix = true;

  static const bool has_static_size = true;

  static const size_t static_rows = N;

  static const size_t static_cols = M;

  static const constexpr StorageLayout storage_layout = StorageLayout::dense_row_major;

  template <size_t ROWS = static_rows, size_t COLS = static_cols, class Field = ScalarType>
  static inline MatrixTypeTemplate<ROWS, COLS, Field>
  create(const size_t rows, const size_t cols, const SparsityPatternDefault& /*pattern*/ = SparsityPatternDefault())
  {
    if (rows != ROWS)
      DUNE_THROW(Exceptions::shapes_do_not_match, "rows = " << rows << "\nN = " << int(N));
    if (cols != COLS)
      DUNE_THROW(Exceptions::shapes_do_not_match, "cols = " << cols << "\nM = " << int(M));
    return MatrixTypeTemplate<ROWS, COLS, Field>();
  }

  template <size_t ROWS = static_rows, size_t COLS = static_cols, class Field = ScalarType>
  static inline MatrixTypeTemplate<ROWS, COLS, Field>
  create(const size_t rows,
         const size_t cols,
         const ScalarType& val,
         const SparsityPatternDefault& /*pattern*/ = SparsityPatternDefault())
  {
    if (rows != ROWS)
      DUNE_THROW(Exceptions::shapes_do_not_match, "rows = " << rows << "\nN = " << int(N));
    if (cols != COLS)
      DUNE_THROW(Exceptions::shapes_do_not_match, "cols = " << cols << "\nM = " << int(M));
    return MatrixTypeTemplate<ROWS, COLS, Field>(val);
  }

  static inline std::unique_ptr<MatrixType> create_dynamic(const size_t rows, const size_t cols)
  {
    if (rows != N)
      DUNE_THROW(Exceptions::shapes_do_not_match, "rows = " << rows << "\nN = " << int(N));
    if (cols != M)
      DUNE_THROW(Exceptions::shapes_do_not_match, "cols = " << cols << "\nM = " << int(M));
    return std::make_unique<MatrixType>();
  }

  static inline std::unique_ptr<MatrixType> create_dynamic(const size_t rows, const size_t cols, const ScalarType& val)
  {
    if (rows != N)
      DUNE_THROW(Exceptions::shapes_do_not_match, "rows = " << rows << "\nN = " << int(N));
    if (cols != M)
      DUNE_THROW(Exceptions::shapes_do_not_match, "cols = " << cols << "\nM = " << int(M));
    return std::make_unique<MatrixType>(val);
  }

  static inline size_t rows(const MatrixType& /*mat*/)
  {
    return numeric_cast<size_t>(N);
  }

  static inline size_t cols(const MatrixType& /*mat*/)
  {
    return numeric_cast<size_t>(M);
  }

  static inline void set_entry(MatrixType& mat, const size_t row, const size_t col, const ScalarType& val)
  {
    mat[row][col] = val;
  }

  static inline ScalarType get_entry(const MatrixType& mat, const size_t row, const size_t col)
  {
    return mat[row][col];
  }

  static inline void add_to_entry(MatrixType& mat, const size_t row, const size_t col, const ScalarType& val)
  {
    mat[row][col] += val;
  }

  static inline ScalarType* data(MatrixType& mat)
  {
    return &(mat[0][0]);
  }

  static inline const ScalarType* data(const MatrixType& mat)
  {
    return &(mat[0][0]);
  }
};


template <class MatrixType>
typename std::enable_if<is_matrix<MatrixType>::value, size_t>::type get_matrix_rows(const MatrixType& matrix)
{
  return MatrixAbstraction<MatrixType>::rows(matrix);
}


template <class MatrixType>
typename std::enable_if<is_matrix<MatrixType>::value, size_t>::type get_matrix_cols(const MatrixType& matrix)
{
  return MatrixAbstraction<MatrixType>::cols(matrix);
}


template <class MatrixType>
typename std::enable_if<is_matrix<MatrixType>::value, typename MatrixAbstraction<MatrixType>::S>::type
get_matrix_entry(const MatrixType& matrix, const size_t ii, const size_t jj)
{
  return MatrixAbstraction<MatrixType>::get_entry(matrix, ii, jj);
}


template <class MatrixType, class S>
typename std::enable_if<is_matrix<MatrixType>::value, void>::type
set_matrix_entry(MatrixType& matrix, const size_t ii, const size_t jj, const S& value)
{
  MatrixAbstraction<MatrixType>::set_entry(matrix, ii, jj, value);
}


template <class MatrixType,
          size_t ROWS = MatrixAbstraction<MatrixType>::static_rows,
          size_t COLS = MatrixAbstraction<MatrixType>::static_cols,
          class Field = typename MatrixAbstraction<MatrixType>::ScalarType>
typename std::enable_if<is_matrix<MatrixType>::value,
                        typename MatrixAbstraction<MatrixType>::template MatrixTypeTemplate<ROWS, COLS, Field>>::type
create(const size_t rows,
       const size_t cols,
       const typename MatrixAbstraction<MatrixType>::S& val = 0,
       const SparsityPatternDefault& pattern = SparsityPatternDefault())
{
  return MatrixAbstraction<MatrixType>::template create<ROWS, COLS, Field>(rows, cols, val, pattern);
}


template <class MatrixType>
typename std::enable_if<is_matrix<MatrixType>::value, MatrixType>::type
create(const size_t rows, const size_t cols, const typename MatrixAbstraction<MatrixType>::S& val = 0)
{
  return MatrixAbstraction<MatrixType>::create(rows, cols, val);
}


template <class TargetMatrixType, class SourceMatrixType>
typename std::enable_if<is_matrix<TargetMatrixType>::value && is_matrix<SourceMatrixType>::value,
                        TargetMatrixType>::type
zeros_like(const SourceMatrixType& source)
{
  return create<TargetMatrixType>(get_matrix_rows(source), get_matrix_cols(source), 0);
}


template <class MatrixType>
typename std::enable_if<is_matrix<MatrixType>::value, MatrixType>::type zeros_like(const MatrixType& source)
{
  return zeros_like<MatrixType, MatrixType>(source);
}


template <class MatrixType>
typename std::enable_if<is_matrix<MatrixType>::value, typename MatrixAbstraction<MatrixType>::ScalarType*>::type
data(MatrixType& source)
{
  return MatrixAbstraction<MatrixType>::data(source);
}


template <class MatrixType>
typename std::enable_if<is_matrix<MatrixType>::value, typename MatrixAbstraction<MatrixType>::ScalarType*>::type
data(const MatrixType& source)
{
  return MatrixAbstraction<MatrixType>::data(source);
}


template <class T, class M>
typename std::enable_if<is_matrix<M>::value && is_arithmetic<T>::value, std::unique_ptr<T[]>>::type
serialize_rowwise(const M& mat)
{
  using Mat = MatrixAbstraction<M>;
  const size_t rows = Mat::rows(mat);
  const size_t cols = Mat::cols(mat);
  auto data = std::make_unique<T[]>(rows * cols);
  size_t ii = 0;
  for (size_t rr = 0; rr < rows; ++rr)
    for (size_t cc = 0; cc < cols; ++cc)
#ifdef DXT_DISABLE_CHECKS
      data[ii++] = Mat::get_entry(mat, rr, cc);
#else
      data[ii++] = numeric_cast<T>(Mat::get_entry(mat, rr, cc));
#endif
  return data;
} // ... serialize_rowwise(...)


template <class M>
typename std::enable_if<is_matrix<M>::value, std::unique_ptr<typename MatrixAbstraction<M>::S[]>>::type
serialize_rowwise(const M& mat)
{
  return serialize_rowwise<typename MatrixAbstraction<M>::ScalarType>(mat);
}


template <class T, class M>
typename std::enable_if<is_matrix<M>::value && is_arithmetic<T>::value, std::unique_ptr<T[]>>::type
serialize_colwise(const M& mat)
{
  using Mat = MatrixAbstraction<M>;
  const size_t rows = Mat::rows(mat);
  const size_t cols = Mat::cols(mat);
  auto data = std::make_unique<T[]>(rows * cols);
  size_t ii = 0;
  for (size_t cc = 0; cc < cols; ++cc)
    for (size_t rr = 0; rr < rows; ++rr)
#ifdef DXT_DISABLE_CHECKS
      data[ii++] = numeric_cast<T>(Mat::get_entry(mat, rr, cc));
#else
      data[ii++] = Mat::get_entry(mat, rr, cc);
#endif
  return data;
} // ... serialize_colwise(...)


template <class M>
typename std::enable_if<is_matrix<M>::value, std::unique_ptr<typename MatrixAbstraction<M>::S[]>>::type
serialize_colwise(const M& mat)
{
  return serialize_colwise<typename MatrixAbstraction<M>::ScalarType>(mat);
}


template <class RangeType, class SourceType>
typename std::enable_if<is_matrix<SourceType>::value && is_matrix<RangeType>::value, RangeType>::type
convert_to(const SourceType& source)
{
  const size_t rows = get_matrix_rows(source);
  const size_t cols = get_matrix_cols(source);
  auto ret = create<RangeType>(rows, cols);
  for (size_t ii = 0; ii < rows; ++ii)
    for (size_t jj = 0; jj < cols; ++jj)
      set_matrix_entry(ret,
                       ii,
                       jj,
#ifndef DXT_DISABLE_CHECKS
                       numeric_cast<typename MatrixAbstraction<RangeType>::S>(
#endif
                           get_matrix_entry(source, ii, jj)
#ifndef DXT_DISABLE_CHECKS
                               )
#endif
                           );
  return ret;
} // ... convert_to(...)


template <class MatrixType, class M = MatrixAbstraction<MatrixType>>
typename std::enable_if<is_matrix<MatrixType>::value,
                        typename M::template MatrixTypeTemplate<M::static_cols, M::static_rows>>::type
transposed(const MatrixType& mat)
{
  auto ret = M::template create<M::static_cols, M::static_rows>(M::cols(mat), M::rows(mat), 0.);
  for (size_t ii = 0; ii < M::rows(mat); ++ii)
    for (size_t jj = 0; jj < M::cols(mat); ++jj)
      set_matrix_entry(ret, jj, ii, get_matrix_entry(mat, ii, jj));
  return ret;
}


} // namespace Common
} // namespace XT


template <class K>
Dune::DynamicMatrix<K> operator*(const Dune::DynamicMatrix<K>& lhs, const Dune::DynamicMatrix<K>& rhs)
{
  Dune::DynamicMatrix<K> ret(lhs.rows(), rhs.cols(), 0.);
  XT::Common::multiply_helper(lhs, rhs, ret);
  return ret;
}

template <class K>
Dune::DynamicMatrix<K> operator+(const Dune::DynamicMatrix<K>& lhs, const Dune::DynamicMatrix<K>& rhs)
{
  Dune::DynamicMatrix<K> ret(lhs);
  ret += rhs;
  return ret;
}

template <class K>
Dune::DynamicMatrix<K> operator-(const Dune::DynamicMatrix<K>& lhs, const Dune::DynamicMatrix<K>& rhs)
{
  Dune::DynamicMatrix<K> ret(lhs);
  ret -= rhs;
  return ret;
}


} // namespace Dune

#endif // DUNE_XT_COMMON_MATRIX_HH
