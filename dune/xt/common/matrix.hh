// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2012 - 2017, 2019)
//   René Fritze     (2010 - 2016, 2018 - 2019)
//   Sven Kaulmann   (2010 - 2011)
//   Tobias Leibner  (2014, 2017 - 2018)

#ifndef DUNE_XT_COMMON_MATRIX_HH
#define DUNE_XT_COMMON_MATRIX_HH

#include <memory>
#include <ostream>

#include <dune/common/dynmatrix.hh>
#include <dune/common/fmatrix.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/common/type_traits.hh>

namespace Dune {
namespace XT {
namespace Common {


// Sparsity pattern placeholder representing a dense matrix pattern.
struct FullPattern
{};

/**
 * \brief Traits to uniformly handle matrices.
 *
 *        If you want your matrix class to benefit from the functionality defined in this header you have to manually
 *        specify a specialization of this class in your code with is_matrix defined to true and the appropriate
 *        static methods implemented and members defined (see the specializations below).
 */
template <class MatType>
struct MatrixAbstraction
{
  using M = std::conditional_t<std::is_same<MatType, void>::value, int, MatType>; // avoid reference to void
  typedef M MatrixType;
  typedef typename Dune::FieldTraits<M>::field_type ScalarType;
  typedef typename Dune::FieldTraits<M>::real_type RealType;
  typedef ScalarType S;
  typedef RealType R;

  template <size_t rows = 0, size_t cols = 0, class FieldType = ScalarType>
  using MatrixTypeTemplate = MatrixType;

  static const bool is_matrix = false;

  static const bool has_static_size = false;

  static const size_t static_rows = std::numeric_limits<size_t>::max();

  static const size_t static_cols = std::numeric_limits<size_t>::max();

  static const constexpr StorageLayout storage_layout = XT::Common::StorageLayout::other;

  static const bool has_ostream = false;

  template <class SparsityPatternType = FullPattern>
  static inline /*MatrixType*/ void create(const size_t /*rows*/,
                                           const size_t /*cols*/,
                                           const ScalarType& /*val*/ = suitable_default<ScalarType>::value(),
                                           const SparsityPatternType& /*pattern*/ = SparsityPatternType())
  {
    static_assert(AlwaysFalse<MatType>::value, "Do not call me if is_matrix is false!");
  }

  template <class SparsityPatternType = FullPattern>
  static inline /*std::unique_ptr<MatrixType>*/ void
  make_unique(const size_t /*rows*/,
              const size_t /*cols*/,
              const ScalarType& /*val*/ = suitable_default<ScalarType>::value(),
              const SparsityPatternType& /*pattern*/ = SparsityPatternType())
  {
    static_assert(AlwaysFalse<MatType>::value, "Do not call me if is_matrix is false!");
  }

  static inline size_t rows(const MatrixType& /*mat*/)
  {
    DUNE_THROW(Dune::InvalidStateException, "Do not call me if is_matrix is false!");
    return 0;
  }

  static inline size_t cols(const MatrixType& /*mat*/)
  {
    DUNE_THROW(Dune::InvalidStateException, "Do not call me if is_matrix is false!");
    return 0;
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

  static inline ScalarType* data(std::remove_const_t<MatrixType>& /*mat*/)
  {
    static_assert(AlwaysFalse<MatType>::value, "Do not call me if storage_layout is not dense!");
    return nullptr;
  }

  template <bool is_mat = is_matrix>
  static inline const std::enable_if_t<is_mat, ScalarType>* data(std::add_const_t<MatrixType>& /*vec*/)
  {
    static_assert(AlwaysFalse<MatType>::value, "Do not call me if storage_layout is not dense!");
    return nullptr;
  }
};

template <class K>
struct MatrixAbstraction<Dune::DynamicMatrix<K>>
{
  typedef Dune::DynamicMatrix<K> MatrixType;
  typedef typename Dune::FieldTraits<K>::field_type ScalarType;
  typedef typename Dune::FieldTraits<K>::real_type RealType;
  typedef ScalarType S;
  typedef RealType R;
  template <size_t rows = 0, size_t cols = 0, class FieldType = K>
  using MatrixTypeTemplate = DynamicMatrix<K>;

  static const bool is_matrix = true;

  static const bool has_static_size = false;

  static const size_t static_rows = std::numeric_limits<size_t>::max();

  static const size_t static_cols = std::numeric_limits<size_t>::max();

  static const constexpr StorageLayout storage_layout = StorageLayout::other;

  static const bool has_ostream = true;

  template <class SparsityPatternType = FullPattern>
  static inline MatrixType create(const size_t rows,
                                  const size_t cols,
                                  const ScalarType& val = suitable_default<ScalarType>::value(),
                                  const SparsityPatternType& /*pattern*/ = SparsityPatternType())
  {
    return MatrixType(rows, cols, val);
  }

  template <class SparsityPatternType = FullPattern>
  static inline std::unique_ptr<MatrixType> make_unique(const size_t rows,
                                                        const size_t cols,
                                                        const ScalarType& val = suitable_default<ScalarType>::value(),
                                                        const SparsityPatternType& /*pattern*/ = SparsityPatternType())
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
  template <size_t rows = N, size_t cols = M, class FieldType = K>
  using MatrixTypeTemplate = Dune::FieldMatrix<FieldType, rows, cols>;

  static const bool is_matrix = true;

  static const bool has_static_size = true;

  static const size_t static_rows = N;

  static const size_t static_cols = M;

  static const constexpr StorageLayout storage_layout = StorageLayout::dense_row_major;

  static const bool has_ostream = true;

  template <class SparsityPatternType = FullPattern>
  static inline MatrixType create(const size_t rows,
                                  const size_t cols,
                                  const ScalarType& val = suitable_default<ScalarType>::value(),
                                  const SparsityPatternType& /*pattern*/ = SparsityPatternType())
  {
    if (rows != N)
      DUNE_THROW(Exceptions::shapes_do_not_match, "rows = " << rows << "\nN = " << int(N));
    if (cols != M)
      DUNE_THROW(Exceptions::shapes_do_not_match, "cols = " << cols << "\nM = " << int(M));
    return MatrixType(val);
  }

  template <class SparsityPatternType = FullPattern>
  static inline std::unique_ptr<MatrixType> make_unique(const size_t rows,
                                                        const size_t cols,
                                                        const ScalarType& val = suitable_default<ScalarType>::value(),
                                                        const SparsityPatternType& /*pattern*/ = SparsityPatternType())
  {
    if (rows != N)
      DUNE_THROW(Exceptions::shapes_do_not_match, "rows = " << rows << "\nN = " << int(N));
    if (cols != M)
      DUNE_THROW(Exceptions::shapes_do_not_match, "cols = " << cols << "\nM = " << int(M));
    return std::make_unique<MatrixType>(val);
  }

  static constexpr size_t rows(const MatrixType& /*mat*/)
  {
    return static_cast<size_t>(N);
  }

  static constexpr size_t cols(const MatrixType& /*mat*/)
  {
    return static_cast<size_t>(M);
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
          class FieldType = typename MatrixAbstraction<MatrixType>::S,
          class SparsityPatternType = FullPattern>
typename std::enable_if<
    is_matrix<MatrixType>::value,
    typename MatrixAbstraction<MatrixType>::template MatrixTypeTemplate<ROWS, COLS, FieldType>>::type
create(const size_t rows,
       const size_t cols,
       const FieldType& val = 0,
       const SparsityPatternType& pattern = SparsityPatternType())
{
  return MatrixAbstraction<
      typename MatrixAbstraction<MatrixType>::template MatrixTypeTemplate<ROWS, COLS, FieldType>>::create(rows,
                                                                                                          cols,
                                                                                                          val,
                                                                                                          pattern);
}


template <class TargetMatrixType, class SourceMatrixType>
typename std::enable_if<is_matrix<TargetMatrixType>::value && is_matrix<SourceMatrixType>::value,
                        TargetMatrixType>::type
zeros_like(const SourceMatrixType& source)
{
  return create<TargetMatrixType>(
      get_matrix_rows(source), get_matrix_cols(source), typename MatrixAbstraction<TargetMatrixType>::S(0));
}


template <class MatrixType>
typename std::enable_if<is_matrix<MatrixType>::value, MatrixType>::type zeros_like(const MatrixType& source)
{
  return zeros_like<MatrixType, MatrixType>(source);
}


template <class MatrixType, class SourceVectorType, class TargetVectorType>
std::enable_if_t<is_matrix<MatrixType>::value && is_vector<SourceVectorType>::value
                     && is_vector<TargetVectorType>::value,
                 void>
mv(const MatrixType& A, const SourceVectorType& x, TargetVectorType& b)
{
  using M = MatrixAbstraction<MatrixType>;
  const size_t rows = M::rows(A);
  const size_t cols = M::rows(A);
  DUNE_THROW_IF(
      x.size() != cols, Exceptions::shapes_do_not_match, "A.cols() = " << cols << "\n x.size() = " << x.size());
  DUNE_THROW_IF(
      b.size() != rows, Exceptions::shapes_do_not_match, "A.rows() = " << rows << "\n b.size() = " << b.size());
  for (size_t ii = 0; ii < rows; ++ii)
    b[ii] = 0;
  for (size_t ii = 0; ii < rows; ++ii)
    for (size_t jj = 0; jj < cols; ++jj)
      b[ii] += M::get_entry(A, ii, jj) * x[jj];
} // ... mv(...)

template <class F>
void mv(const DynamicMatrix<F>& A, const DynamicVector<F>& x, DynamicVector<F>& b)
{
  A.mv(x, b);
}

template <class F, int ROWS, int COLS>
void mv(const FieldMatrix<F, ROWS, COLS>& A, const FieldVector<F, COLS>& x, FieldVector<F, ROWS>& b)
{
  A.mv(x, b);
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
  auto ret = create<MatrixType, M::static_cols, M::static_rows>(M::cols(mat), M::rows(mat), typename M::S(0.));
  for (size_t ii = 0; ii < M::rows(mat); ++ii)
    for (size_t jj = 0; jj < M::cols(mat); ++jj)
      set_matrix_entry(ret, jj, ii, get_matrix_entry(mat, ii, jj));
  return ret;
}


template <class M, class CharType, class CharTraits>
std::basic_ostream<CharType, CharTraits>& output_matrix(std::basic_ostream<CharType, CharTraits>& out, const M& mat)
{
  using Matrix = Dune::XT::Common::MatrixAbstraction<M>;
  out << "[";
  for (size_t ii = 0; ii < Matrix::rows(mat); ++ii) {
    for (size_t jj = 0; jj < Matrix::cols(mat); ++jj) {
      out << Matrix::get_entry(mat, ii, jj);
      if (jj == Matrix::cols(mat) - 1 && ii < Matrix::rows(mat) - 1)
        out << ";";
      else if (jj < Matrix::cols(mat) - 1)
        out << " ";
    }
  }
  out << "]";
  return out;
}


} // namespace Common
} // namespace XT


template <class K>
Dune::DynamicMatrix<K> operator*(const Dune::DynamicMatrix<K>& lhs, const Dune::DynamicMatrix<K>& rhs)
{
  Dune::DynamicMatrix<K> ret(lhs.rows(), rhs.cols(), 0.);
  for (size_t ii = 0; ii < lhs.rows(); ++ii) {
    for (size_t jj = 0; jj < rhs.cols(); ++jj) {
      ret[ii][jj] = 0.;
      for (size_t kk = 0; kk < lhs.cols(); ++kk)
        ret[ii][jj] += lhs[ii][kk] * rhs[kk][jj];
    }
  }
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
