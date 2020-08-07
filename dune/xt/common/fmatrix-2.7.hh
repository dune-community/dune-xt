// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2018)
//   René Fritze     (2015 - 2016, 2018 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2014, 2018 - 2019)

#ifndef DUNE_XT_COMMON_FMATRIX_27_HH
#define DUNE_XT_COMMON_FMATRIX_27_HH

#include <initializer_list>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/debug.hh>
#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/type_traits.hh>

namespace Dune {
namespace XT {
namespace Common {


/**
 * \todo We need to implement all operators from the base which return the base, to rather return ourselfes!
 */
template <class K, int ROWS, int COLS>
class FieldMatrix : public Dune::FieldMatrix<K, ROWS, COLS>
{
  using BaseType = Dune::FieldMatrix<K, ROWS, COLS>;
  using ThisType = FieldMatrix;

public:
  using typename BaseType::field_type;
  using typename BaseType::size_type;
  using typename BaseType::value_type;

  FieldMatrix(const K& kk = suitable_default<K>::value())
    : BaseType()
  {
    // This is required because BaseType(kk) does not work for std::string
    for (size_t i = 0; i < ROWS; ++i) {
      for (size_t j = 0; j < COLS; ++j)
        (*this)[i][j] = kk;
    }
  }

  FieldMatrix(const size_t DXTC_DEBUG_ONLY(rr),
              const size_t DXTC_DEBUG_ONLY(cc),
              const K& kk = suitable_default<K>::value())
    : BaseType()
  {
#ifndef NDEBUG
    if (rr != ROWS || cc != COLS)
      DUNE_THROW(Exceptions::wrong_input_given,
                 "You are trying to construct a FieldMatrix< ..., " << ROWS << ", " << COLS << " > (of "
                                                                    << "static size) with " << rr << " rows and " << cc
                                                                    << " columns!");
#endif // NDEBUG
    // This is required because BaseType(kk) does not work for std::string
    for (size_t i = 0; i < ROWS; ++i) {
      for (size_t j = 0; j < COLS; ++j)
        (*this)[i][j] = kk;
    }

  } // ... FieldMatrix(...)

  FieldMatrix(std::initializer_list<std::initializer_list<K>> list_of_rows)
    : BaseType()
  {
#ifndef NDEBUG
    if (list_of_rows.size() != ROWS)
      DUNE_THROW(Exceptions::wrong_input_given,
                 "You are trying to construct a FieldMatrix< ..., "
                     << ROWS << ", " << COLS << " > (of "
                     << "static size) from a list modeling a matrix with " << list_of_rows.size() << " rows!");
#endif // NDEBUG
    size_t rr = 0;
    for (auto row : list_of_rows) {
#ifndef NDEBUG
      if (row.size() != COLS)
        DUNE_THROW(Exceptions::wrong_input_given,
                   "You are trying to construct a FieldMatrix< ..., "
                       << ROWS << ", " << COLS << " > (of "
                       << "static size) from a list with a row of length " << row.size() << "!");
#endif // NDEBUG
      size_t cc = 0;
      for (auto entry : row)
        (*this)[rr][cc++] = entry;
      ++rr;
    }
  } // FieldMatrix(std::initializer_list<std::initializer_list<K>> list_of_rows)

  template <class OtherMatrixType>
  FieldMatrix(const OtherMatrixType& other, std::enable_if_t<is_matrix<OtherMatrixType>::value, int> /*dummy*/ = 0)
  {
    *this = other;
  }

  template <class OtherMatrixType>
  std::enable_if_t<is_matrix<OtherMatrixType>::value, ThisType>& operator=(const OtherMatrixType& other)
  {
    for (size_t rr = 0; rr < ROWS; ++rr)
      for (size_t cc = 0; cc < COLS; ++cc)
        (*this)[rr][cc] = MatrixAbstraction<OtherMatrixType>::get_entry(other, rr, cc);
    return *this;
  }

  Dune::XT::Common::FieldMatrix<K, COLS, ROWS> transpose() const
  {
    Dune::XT::Common::FieldMatrix<K, COLS, ROWS> ret;
    for (size_t rr = 0; rr < ROWS; ++rr)
      for (size_t cc = 0; cc < COLS; ++cc)
        ret[cc][rr] = (*this)[rr][cc];
    return ret;
  }

  Dune::XT::Common::FieldVector<K, ROWS> operator*(const Dune::FieldVector<K, COLS>& vec) const
  {
    Dune::FieldVector<K, ROWS> ret;
    this->mv(vec, ret);
    return ret;
  }

  //! vector space addition -- two-argument version
  template <class OtherScalar>
  friend auto operator+(const ThisType& matrixA, const Dune::FieldMatrix<OtherScalar, ROWS, COLS>& matrixB)
  {
    FieldMatrix<typename PromotionTraits<K, OtherScalar>::PromotedType, ROWS, COLS> result;
    for (size_type i = 0; i < ROWS; ++i)
      for (size_type j = 0; j < COLS; ++j)
        result[i][j] = matrixA[i][j] + matrixB[i][j];
    return result;
  }

  //! vector space subtraction -- two-argument version
  template <class OtherScalar>
  friend auto operator-(const ThisType& matrixA, const Dune::FieldMatrix<OtherScalar, ROWS, COLS>& matrixB)
  {
    FieldMatrix<typename PromotionTraits<K, OtherScalar>::PromotedType, ROWS, COLS> result;
    for (size_type i = 0; i < ROWS; ++i)
      for (size_type j = 0; j < COLS; ++j)
        result[i][j] = matrixA[i][j] - matrixB[i][j];
    return result;
  }

  //! vector space multiplication with scalar
  template <class Scalar, std::enable_if_t<IsNumber<Scalar>::value, int> = 0>
  friend auto operator*(const ThisType& matrix, Scalar scalar)
  {
    FieldMatrix<typename PromotionTraits<K, Scalar>::PromotedType, ROWS, COLS> result;
    for (size_type i = 0; i < ROWS; ++i)
      for (size_type j = 0; j < COLS; ++j)
        result[i][j] = matrix[i][j] * scalar;
    return result;
  }

  //! vector space multiplication with scalar
  template <class Scalar, std::enable_if_t<IsNumber<Scalar>::value, int> = 0>
  friend auto operator*(Scalar scalar, const ThisType& matrix)
  {
    FieldMatrix<typename PromotionTraits<K, Scalar>::PromotedType, ROWS, COLS> result;
    for (size_type i = 0; i < ROWS; ++i)
      for (size_type j = 0; j < COLS; ++j)
        result[i][j] = scalar * matrix[i][j];
    return result;
  }

  //! vector space division by scalar
  template <class Scalar, std::enable_if_t<IsNumber<Scalar>::value, int> = 0>
  friend auto operator/(const ThisType& matrix, Scalar scalar)
  {
    FieldMatrix<typename PromotionTraits<K, Scalar>::PromotedType, ROWS, COLS> result;
    for (size_type i = 0; i < ROWS; ++i)
      for (size_type j = 0; j < COLS; ++j)
        result[i][j] = matrix[i][j] / scalar;
    return result;
  }

  /** \brief Matrix-matrix multiplication
   */
  template <class OtherScalar, int otherCols>
  friend auto operator*(const ThisType& matrixA, const Dune::FieldMatrix<OtherScalar, COLS, otherCols>& matrixB)
  {
    FieldMatrix<typename PromotionTraits<K, OtherScalar>::PromotedType, ROWS, otherCols> result;

    for (size_type i = 0; i < matrixA.mat_rows(); ++i)
      for (size_type j = 0; j < matrixB.mat_cols(); ++j) {
        result[i][j] = 0;
        for (size_type k = 0; k < matrixA.mat_cols(); ++k)
          result[i][j] += matrixA[i][k] * matrixB[k][j];
      }

    return result;
  }

  template <class OtherScalar, int otherRows>
  friend auto operator*(const Dune::FieldMatrix<OtherScalar, otherRows, ROWS>& matrixA, const ThisType& matrixB)
  {
    FieldMatrix<typename PromotionTraits<K, OtherScalar>::PromotedType, otherRows, COLS> result;
    for (size_type i = 0; i < matrixA.mat_rows(); ++i)
      for (size_type j = 0; j < matrixB.mat_cols(); ++j) {
        result[i][j] = 0;
        for (size_type k = 0; k < matrixA.mat_cols(); ++k)
          result[i][j] += matrixA[i][k] * matrixB[k][j];
      }
    return result;
  }

  template <class OtherScalar, int otherCols>
  friend auto operator*(const ThisType& matrixA, const FieldMatrix<OtherScalar, COLS, otherCols>& matrixB)
  {
    FieldMatrix<typename PromotionTraits<K, OtherScalar>::PromotedType, ROWS, otherCols> result;
    for (size_type i = 0; i < matrixA.mat_rows(); ++i)
      for (size_type j = 0; j < matrixB.mat_cols(); ++j) {
        result[i][j] = 0;
        for (size_type k = 0; k < matrixA.mat_cols(); ++k)
          result[i][j] += matrixA[i][k] * matrixB[k][j];
      }
    return result;
  }
}; // class FieldMatrix<...>


/**
 * \todo We need to implement all operators from the base which return the base, to rather return ourselfes!
 */
template <class K>
class FieldMatrix<K, 1, 1> : public Dune::FieldMatrix<K, 1, 1>
{
  static const int ROWS = 1;
  static const int COLS = 1;
  using BaseType = Dune::FieldMatrix<K, ROWS, COLS>;
  using ThisType = FieldMatrix;

public:
  using typename BaseType::size_type;

  FieldMatrix(const K& kk = suitable_default<K>::value())
    : BaseType()
  {
    (*this)[0][0] = kk;
  }

  FieldMatrix(std::initializer_list<std::initializer_list<K>> list_of_rows)
    : BaseType()
  {
#ifndef NDEBUG
    if (list_of_rows.size() != ROWS)
      DUNE_THROW(Exceptions::wrong_input_given,
                 "You are trying to construct a FieldMatrix< ..., 1, 1 > (of "
                     << "static size) from a list modeling a matrix with " << list_of_rows.size() << " rows!");
#endif // NDEBUG
    for (auto row : list_of_rows) {
#ifndef NDEBUG
      if (row.size() != COLS)
        DUNE_THROW(Exceptions::wrong_input_given,
                   "You are trying to construct a FieldMatrix< ..., 1, 1 > (of "
                       << "static size) from a list with a row of length " << row.size() << "!");
#endif // NDEBUG
      for (auto entry : row)
        (*this)[0][0] = entry;
    }
  } // FieldMatrix(std::initializer_list<std::initializer_list<K>> list_of_rows)

  FieldMatrix(const size_t DXTC_DEBUG_ONLY(rr),
              const size_t DXTC_DEBUG_ONLY(cc),
              const K& kk = suitable_default<K>::value())
    : BaseType()
  {
#ifndef NDEBUG
    if (rr != ROWS || cc != COLS)
      DUNE_THROW(Exceptions::wrong_input_given,
                 "You are trying to construct a FieldMatrix< ..., " << ROWS << ", " << COLS << " > (of "
                                                                    << "static size) with " << rr << " rows and " << cc
                                                                    << " columns!");
#endif // NDEBUG
    (*this)[0][0] = kk;
  } // ... FieldMatrix(...)

  FieldMatrix(const BaseType& other)
    : BaseType(other)
  {}

  FieldMatrix(const Dune::XT::Common::FieldVector<K, 1>& other)
    : BaseType()
  {
    (*this)[0][0] = other[0];
  }

  FieldMatrix(const Dune::FieldVector<K, 1>& other)
    : BaseType()
  {
    (*this)[0][0] = other[0];
  }

  Dune::XT::Common::FieldMatrix<K, COLS, ROWS> transpose() const
  {
    return *this;
  }

  using BaseType::operator=;

  ThisType& operator=(const FieldVector<K, 1>& other)
  {
    BaseType::operator=(other[0]);
    return *this;
  }

  //! vector space addition -- two-argument version
  template <class OtherScalar>
  friend auto operator+(const ThisType& matrixA, const Dune::FieldMatrix<OtherScalar, 1, 1>& matrixB)
  {
    return FieldMatrix<typename PromotionTraits<K, OtherScalar>::PromotedType, 1, 1>{matrixA[0][0] + matrixB[0][0]};
  }

  //! Binary addition when treating FieldMatrix<K,1,1> like K
  template <class Scalar, std::enable_if_t<IsNumber<Scalar>::value, int> = 0>
  friend auto operator+(const ThisType& matrix, const Scalar& scalar)
  {
    return FieldMatrix<typename PromotionTraits<K, Scalar>::PromotedType, 1, 1>{matrix[0][0] + scalar};
  }

  //! Binary addition when treating FieldMatrix<K,1,1> like K
  template <class Scalar, std::enable_if_t<IsNumber<Scalar>::value, int> = 0>
  friend auto operator+(const Scalar& scalar, const ThisType& matrix)
  {
    return FieldMatrix<typename PromotionTraits<Scalar, K>::PromotedType, 1, 1>{scalar + matrix[0][0]};
  }

  //! vector space subtraction -- two-argument version
  template <class OtherScalar>
  friend auto operator-(const ThisType& matrixA, const Dune::FieldMatrix<OtherScalar, 1, 1>& matrixB)
  {
    return FieldMatrix<typename PromotionTraits<K, OtherScalar>::PromotedType, 1, 1>{matrixA[0][0] - matrixB[0][0]};
  }

  //! Binary subtraction when treating FieldMatrix<K,1,1> like K
  template <class Scalar, std::enable_if_t<IsNumber<Scalar>::value, int> = 0>
  friend auto operator-(const ThisType& matrix, const Scalar& scalar)
  {
    return FieldMatrix<typename PromotionTraits<K, Scalar>::PromotedType, 1, 1>{matrix[0][0] - scalar};
  }

  //! Binary subtraction when treating FieldMatrix<K,1,1> like K
  template <class Scalar, std::enable_if_t<IsNumber<Scalar>::value, int> = 0>
  friend auto operator-(const Scalar& scalar, const ThisType& matrix)
  {
    return FieldMatrix<typename PromotionTraits<Scalar, K>::PromotedType, 1, 1>{scalar - matrix[0][0]};
  }

  //! vector space multiplication with scalar
  template <class Scalar, std::enable_if_t<IsNumber<Scalar>::value, int> = 0>
  friend auto operator*(const ThisType& matrix, Scalar scalar)
  {
    return FieldMatrix<typename PromotionTraits<K, Scalar>::PromotedType, 1, 1>{matrix[0][0] * scalar};
  }

  //! vector space multiplication with scalar
  template <class Scalar, std::enable_if_t<IsNumber<Scalar>::value, int> = 0>
  friend auto operator*(Scalar scalar, const ThisType& matrix)
  {
    return FieldMatrix<typename PromotionTraits<K, Scalar>::PromotedType, 1, 1>{scalar * matrix[0][0]};
  }

  //! vector space division by scalar
  template <class Scalar, std::enable_if_t<IsNumber<Scalar>::value, int> = 0>
  friend auto operator/(const ThisType& matrix, Scalar scalar)
  {
    return FieldMatrix<typename PromotionTraits<K, Scalar>::PromotedType, 1, 1>{matrix[0][0] / scalar};
  }

  /** \brief Matrix-matrix multiplication
   */
  template <class OtherScalar, int otherCols>
  friend auto operator*(const ThisType& matrixA, const Dune::FieldMatrix<OtherScalar, 1, otherCols>& matrixB)
  {
    FieldMatrix<typename PromotionTraits<K, OtherScalar>::PromotedType, 1, otherCols> result;
    for (size_type j = 0; j < matrixB.mat_cols(); ++j)
      result[0][j] = matrixA[0][0] * matrixB[0][j];
    return result;
  }

  template <class OtherScalar, int otherRows>
  friend auto operator*(const Dune::FieldMatrix<OtherScalar, otherRows, 1>& matrixA, const ThisType& matrixB)
  {
    FieldMatrix<typename PromotionTraits<K, OtherScalar>::PromotedType, otherRows, 1> result;
    for (size_type i = 0; i < matrixA.mat_rows(); ++i)
      result[i][0] = matrixA[i][0] * matrixB[0][0];
    return result;
  }

  template <class OtherScalar, int otherCols>
  friend auto operator*(const ThisType& matrixA, const FieldMatrix<OtherScalar, 1, otherCols>& matrixB)
  {
    FieldMatrix<typename PromotionTraits<K, OtherScalar>::PromotedType, 1, otherCols> result;
    for (size_type j = 0; j < matrixB.mat_cols(); ++j)
      result[0][j] = matrixA[0][0] * matrixB[0][j];
    return result;
  }
}; // class FieldMatrix


template <class K, size_t num_blocks, size_t block_rows, size_t block_cols = block_rows>
class BlockedFieldMatrix
{
  using ThisType = BlockedFieldMatrix;

public:
  static constexpr size_t num_rows = num_blocks * block_rows;
  static constexpr size_t num_cols = num_blocks * block_cols;
  using MatrixType = Dune::FieldMatrix<K, num_rows, num_cols>;
  using BlockType = FieldMatrix<K, block_rows, block_cols>;

  BlockedFieldMatrix(const K& val = K(0.))
    : backend_(BlockType(val))
  {}

  BlockedFieldMatrix(const size_t DXTC_DEBUG_ONLY(rows), const size_t DXTC_DEBUG_ONLY(cols), const K& val = K(0.))
    : backend_(BlockType(val))
  {
    assert(rows == num_rows && cols == num_cols && "Requested shape has to match static shape!");
  }

  template <class OtherMatrixType>
  BlockedFieldMatrix(const OtherMatrixType& other,
                     std::enable_if_t<is_matrix<OtherMatrixType>::value, int> /*dummy*/ = 0)
  {
    *this = other;
  }

  template <class OtherMatrixType>
  std::enable_if_t<is_matrix<OtherMatrixType>::value, ThisType>& operator=(const OtherMatrixType& other)
  {
    for (size_t jj = 0; jj < num_blocks; ++jj)
      for (size_t ll = 0; ll < block_rows; ++ll)
        for (size_t mm = 0; mm < block_cols; ++mm)
          backend_[jj][ll][mm] =
              MatrixAbstraction<OtherMatrixType>::get_entry(other, jj * block_rows + ll, jj * block_cols + mm);
    return *this;
  }


  BlockedFieldMatrix(const MatrixType& other)
  {
    *this = other;
  }

  BlockedFieldMatrix(const BlockType& block)
    : backend_(block)
  {}

  bool operator==(const ThisType& other) const
  {
    for (size_t jj = 0; jj < num_blocks; ++jj)
      if (block(jj) != other.block(jj))
        return false;
    return true;
  }

  K get_entry(const size_t ii, const size_t jj) const
  {
    assert(ii < num_rows && jj < num_cols);
    if (!in_pattern(ii, jj))
      return K(0.);
    return backend_[ii / block_rows][ii % block_rows][jj % block_cols];
  }

  K& get_entry(const size_t jj, const size_t ll, const size_t mm)
  {
    assert(is_valid_entry(jj, ll, mm));
    return backend_[jj][ll][mm];
  }

  const K& get_entry(const size_t jj, const size_t ll, const size_t mm) const
  {
    assert(is_valid_entry(jj, ll, mm));
    return backend_[jj][ll][mm];
  }

  void set_entry(const size_t ii, const size_t jj, const K& val)
  {
    assert(ii < num_rows && jj < num_cols);
    if (in_pattern(ii, jj))
      backend_[ii / block_rows][ii % block_rows][jj % block_cols] = val;
    else if (XT::Common::FloatCmp::ne(val, K(0)))
      DUNE_THROW(Dune::MathError, "Tried to modify a value that is not in the pattern!");
  }

  void set_entry(const size_t jj, const size_t ll, const size_t mm, const K& val)
  {
    assert(is_valid_entry(jj, ll, mm));
    backend_[jj][ll][mm] = val;
  }

  void add_to_entry(const size_t ii, const size_t jj, const K& val)
  {
    assert(ii < num_rows && jj < num_cols);
    if (!in_pattern(ii, jj) && XT::Common::FloatCmp::ne(val, K(0)))
      DUNE_THROW(Dune::MathError, "Tried to modify a value that is not in the pattern!");
    backend_[ii / block_rows][ii % block_rows][jj % block_cols] += val;
  }

  void add_to_entry(const size_t jj, const size_t ll, const size_t mm, const K& val)
  {
    assert(is_valid_entry(jj, ll, mm));
    backend_[jj][ll][mm] += val;
  }

  BlockType& block(const size_t jj)
  {
    assert(jj < num_blocks);
    return backend_[jj];
  }

  const BlockType& block(const size_t jj) const
  {
    assert(jj < num_blocks);
    return backend_[jj];
  }

  void mv(const Dune::FieldVector<K, num_cols>& x, Dune::FieldVector<K, num_rows>& ret) const
  {
    std::fill(ret.begin(), ret.end(), 0.);
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const auto row_offset = block_rows * jj;
      const auto col_offset = block_cols * jj;
      for (size_t ll = 0; ll < block_rows; ++ll)
        for (size_t mm = 0; mm < block_cols; ++mm)
          ret[row_offset + ll] += backend_[jj][ll][mm] * x[col_offset + mm];
    } // jj
  } // void mv(...)

  void mv(const BlockedFieldVector<K, num_blocks, block_cols>& x,
          BlockedFieldVector<K, num_blocks, block_rows>& ret) const
  {
    for (size_t jj = 0; jj < num_blocks; ++jj)
      backend_[jj].mv(x.block(jj), ret.block(jj));
  } // void mv(...)

  void mtv(const Dune::FieldVector<K, num_rows>& x, Dune::FieldVector<K, num_cols>& ret) const
  {
    std::fill(ret.begin(), ret.end(), 0.);
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const auto row_offset = block_rows * jj;
      const auto col_offset = block_cols * jj;
      for (size_t mm = 0; mm < block_cols; ++mm)
        for (size_t ll = 0; ll < block_rows; ++ll)
          ret[col_offset + mm] += backend_[jj][ll][mm] * x[row_offset + ll];
    } // jj
  } // void mtv(...)

  void mtv(const BlockedFieldVector<K, num_blocks, block_rows>& x,
           BlockedFieldVector<K, num_blocks, block_cols>& ret) const
  {
    for (size_t jj = 0; jj < num_blocks; ++jj)
      backend_[jj].mtv(x.block(jj), ret.block(jj));
  } // void mv(...)

  template <size_t br, size_t bc>
  ThisType& rightmultiply(const BlockedFieldMatrix<K, num_blocks, br, bc>& other)
  {
    assert((this != &other) && "Multiplying a matrix by itself gives wrong results, please copy before!");
    static_assert(br == bc, "Cannot rightmultiply with non-square matrix");
    static_assert(br == block_cols, "Size mismatch");
    for (size_t jj = 0; jj < num_blocks; ++jj)
      backend_[jj].rightmultiply(other.backend_[jj]);
    return *this;
  }

  BlockedFieldMatrix<K, num_blocks, block_cols, block_rows> transpose()
  {
    BlockedFieldMatrix<K, num_blocks, block_cols, block_rows> ret;
    for (size_t jj = 0; jj < num_blocks; ++jj)
      ret.block(jj) = block(jj).transpose();
    return ret;
  }

  ThisType operator*(const ThisType& other) const
  {
    ThisType ret(*this);
    ret.rightmultiply(other);
    return ret;
  }

  ThisType& operator*=(const K& val)
  {
    for (size_t jj = 0; jj < num_blocks; ++jj)
      block(jj) *= val;
    return *this;
  }

  ThisType& operator+=(const ThisType& other)
  {
    for (size_t jj = 0; jj < num_blocks; ++jj)
      block(jj) += other.block(jj);
    return *this;
  }

  ThisType operator+(const ThisType& other) const
  {
    ThisType ret(*this);
    ret += other;
    return ret;
  }

  ThisType& operator-=(const ThisType& other)
  {
    for (size_t jj = 0; jj < num_blocks; ++jj)
      block(jj) -= other.block(jj);
    return *this;
  }

  ThisType operator-(const ThisType& other) const
  {
    ThisType ret(*this);
    ret -= other;
    return ret;
  }

  static bool in_pattern(const size_t ii, const size_t jj)
  {
    return (ii / block_rows == jj / block_cols);
  }

  static bool is_valid_entry(const size_t jj, const size_t ll, const size_t mm)
  {
    return (jj < num_blocks && ll < block_rows && mm < block_cols);
  }

  DynamicMatrix<K> convert_to_dynamic_matrix() const
  {
    DynamicMatrix<K> ret(num_rows, num_cols, 0.);
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const size_t row_offset = jj * block_rows;
      const size_t col_offset = jj * block_cols;
      for (size_t rr = 0; rr < block_rows; ++rr)
        for (size_t cc = 0; cc < block_cols; ++cc)
          ret[row_offset + rr][col_offset + cc] = block(jj)[rr][cc];
    } // jj
    return ret;
  }

  template <class CharType, class CharTraits>
  friend std::basic_ostream<CharType, CharTraits>& operator<<(std::basic_ostream<CharType, CharTraits>& out,
                                                              const ThisType& mat)
  {
    return output_matrix(out, mat);
  } // ... operator<<(...)

private:
  FieldVector<BlockType, num_blocks> backend_;
};


template <class K, int N, int M>
struct MatrixAbstraction<Dune::XT::Common::FieldMatrix<K, N, M>>
{
  typedef Dune::XT::Common::FieldMatrix<K, N, M> MatrixType;
  using ScalarType = typename Dune::FieldTraits<K>::field_type;
  using RealType = typename Dune::FieldTraits<K>::real_type;
  using S = ScalarType;
  using R = RealType;
  template <size_t rows = N, size_t cols = M, class FieldType = K>
  using MatrixTypeTemplate = Dune::XT::Common::FieldMatrix<FieldType, rows, cols>;

  static const bool is_matrix = true;

  static const bool has_static_size = true;

  static const size_t static_rows = N;

  static const size_t static_cols = M;

  static const constexpr StorageLayout storage_layout = StorageLayout::dense_row_major;

  static constexpr bool has_ostream = true;

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

  static constexpr size_t rows(const MatrixType& /*mat*/)
  {
    return N;
  }

  static constexpr size_t cols(const MatrixType& /*mat*/)
  {
    return M;
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

template <class K, size_t num_blocks, size_t block_rows, size_t block_cols>
struct MatrixAbstraction<Dune::XT::Common::BlockedFieldMatrix<K, num_blocks, block_rows, block_cols>>
{
  typedef Dune::XT::Common::BlockedFieldMatrix<K, num_blocks, block_rows, block_cols> MatrixType;
  using ScalarType = typename Dune::FieldTraits<K>::field_type;
  using RealType = typename Dune::FieldTraits<K>::real_type;
  using S = ScalarType;
  using R = RealType;

  static const bool is_matrix = true;

  static const bool has_static_size = true;

  static const size_t static_rows = MatrixType::num_rows;

  static const size_t static_cols = MatrixType::num_cols;

  template <size_t rows = static_rows, size_t cols = static_cols, class FieldType = K>
  using MatrixTypeTemplate = Dune::XT::Common::BlockedFieldMatrix<FieldType, rows / block_rows, block_rows, block_cols>;

  static const constexpr StorageLayout storage_layout = StorageLayout::other;

  static constexpr bool has_ostream = true;

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

  static constexpr size_t rows(const MatrixType& /*mat*/)
  {
    return static_rows;
  }

  static constexpr size_t cols(const MatrixType& /*mat*/)
  {
    return static_cols;
  }

  static inline void set_entry(MatrixType& mat, const size_t row, const size_t col, const ScalarType& val)
  {
    mat.set_entry(row, col, val);
  }

  static inline ScalarType get_entry(const MatrixType& mat, const size_t row, const size_t col)
  {
    return mat.get_entry(row, col);
  }

  static inline void add_to_entry(MatrixType& mat, const size_t row, const size_t col, const ScalarType& val)
  {
    mat.add_to_entry(row, col, val);
  }

  static inline ScalarType* data(MatrixType& /*mat*/)
  {
    DUNE_THROW(InvalidStateException, "Do not call me if storage layout is not dense!");
    return nullptr;
  }

  static inline const ScalarType* data(const MatrixType& /*mat*/)
  {
    DUNE_THROW(InvalidStateException, "Do not call me if storage_layout is not dense!");
    return nullptr;
  }
};


template <class M>
typename std::enable_if<is_matrix<M>::value && MatrixAbstraction<M>::has_static_size,
                        std::unique_ptr<FieldMatrix<typename MatrixAbstraction<M>::S,
                                                    MatrixAbstraction<M>::static_rows,
                                                    MatrixAbstraction<M>::static_cols>>>::type
make_field_container_ptr(const M& mat)
{
  static const size_t rows = MatrixAbstraction<M>::static_rows;
  static const size_t cols = MatrixAbstraction<M>::static_cols;
  auto ret = std::make_unique<FieldMatrix<typename MatrixAbstraction<M>::S, rows, cols>>;
  for (size_t ii = 0; ii < rows; ++ii)
    for (size_t jj = 0; jj < cols; ++jj)
      (*ret)[ii][jj] = get_matrix_entry(mat, ii, jj);
  return std::move(ret);
}


template <class M>
typename std::enable_if<is_matrix<M>::value && MatrixAbstraction<M>::has_static_size,
                        FieldMatrix<typename MatrixAbstraction<M>::S,
                                    MatrixAbstraction<M>::static_rows,
                                    MatrixAbstraction<M>::static_cols>>::type
make_field_container(const M& mat)
{
  static const size_t rows = MatrixAbstraction<M>::static_rows;
  static const size_t cols = MatrixAbstraction<M>::static_cols;
  FieldMatrix<typename MatrixAbstraction<M>::S, rows, cols> ret;
  for (size_t ii = 0; ii < rows; ++ii)
    for (size_t jj = 0; jj < cols; ++jj)
      ret[ii][jj] = get_matrix_entry(mat, ii, jj);
  return ret;
}


template <class K, int ROWS, int COLS>
FieldMatrix<K, ROWS, COLS> make_field_container(Dune::FieldMatrix<K, ROWS, COLS>&& vec)
{
  return std::move(vec);
}


template <class K, int ROWS, int COLS>
FieldMatrix<K, ROWS, COLS> make_field_container(FieldMatrix<K, ROWS, COLS>&& vec)
{
  return std::move(vec);
}


template <class K, int ROWS, int COLS>
typename std::enable_if<is_arithmetic<K>::value && !is_complex<K>::value, FieldMatrix<K, ROWS, COLS>>::type
real(const Dune::FieldMatrix<K, ROWS, COLS>& real_mat)
{
  return real_mat;
}

template <class K, int ROWS, int COLS>
typename std::enable_if<is_arithmetic<K>::value && !is_complex<K>::value, FieldMatrix<K, ROWS, COLS>>::type
real(const FieldMatrix<K, ROWS, COLS>& real_mat)
{
  return real_mat;
}

template <class K, int ROWS, int COLS>
typename std::enable_if<is_arithmetic<K>::value && !is_complex<K>::value, FieldMatrix<K, ROWS, COLS>>::type
real(Dune::FieldMatrix<K, ROWS, COLS>&& real_mat)
{
  return std::move(real_mat);
}

template <class K, int ROWS, int COLS>
typename std::enable_if<is_arithmetic<K>::value && !is_complex<K>::value, FieldMatrix<K, ROWS, COLS>>::type
real(FieldMatrix<K, ROWS, COLS>&& real_mat)
{
  return std::move(real_mat);
}

template <class K, int ROWS, int COLS>
typename std::enable_if<is_complex<K>::value, FieldMatrix<real_t<K>, ROWS, COLS>>::type
real(const Dune::FieldMatrix<K, ROWS, COLS>& complex_mat)
{
  FieldMatrix<real_t<K>, ROWS, COLS> real_mat;
  for (size_t ii = 0; ii < ROWS; ++ii)
    for (size_t jj = 0; jj < COLS; ++jj)
      real_mat[ii][jj] = complex_mat[ii][jj].real();
  return real_mat;
}

template <class K, int ROWS, int COLS>
typename std::enable_if<is_complex<K>::value, FieldMatrix<real_t<K>, ROWS, COLS>>::type
real(const FieldMatrix<K, ROWS, COLS>& complex_mat)
{
  FieldMatrix<real_t<K>, ROWS, COLS> real_mat;
  for (size_t ii = 0; ii < ROWS; ++ii)
    for (size_t jj = 0; jj < COLS; ++jj)
      real_mat[ii][jj] = complex_mat[ii][jj].real();
  return real_mat;
}


template <class K, int ROWS, int COLS>
typename std::enable_if<is_arithmetic<K>::value && !is_complex<K>::value, FieldMatrix<K, ROWS, COLS>>::type
imag(const Dune::FieldMatrix<K, ROWS, COLS>& /*real_mat*/)
{
  return FieldMatrix<K, ROWS, COLS>(0);
}

template <class K, int ROWS, int COLS>
typename std::enable_if<is_arithmetic<K>::value && !is_complex<K>::value, FieldMatrix<K, ROWS, COLS>>::type
imag(const FieldMatrix<K, ROWS, COLS>& /*real_mat*/)
{
  return FieldMatrix<K, ROWS, COLS>(0);
}

template <class K, int ROWS, int COLS>
typename std::enable_if<is_complex<K>::value, FieldMatrix<real_t<K>, ROWS, COLS>>::type
imag(const Dune::FieldMatrix<K, ROWS, COLS>& complex_mat)
{
  FieldMatrix<real_t<K>, ROWS, COLS> real_mat;
  for (size_t ii = 0; ii < ROWS; ++ii)
    for (size_t jj = 0; jj < COLS; ++jj)
      real_mat[ii][jj] = complex_mat[ii][jj].imag();
  return real_mat;
}

template <class K, int ROWS, int COLS>
typename std::enable_if<is_complex<K>::value, FieldMatrix<real_t<K>, ROWS, COLS>>::type
imag(const FieldMatrix<K, ROWS, COLS>& complex_mat)
{
  FieldMatrix<real_t<K>, ROWS, COLS> real_mat;
  for (size_t ii = 0; ii < ROWS; ++ii)
    for (size_t jj = 0; jj < COLS; ++jj)
      real_mat[ii][jj] = complex_mat[ii][jj].imag();
  return real_mat;
}


} // namespace Common
} // namespace XT

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


} // namespace Dune

#endif // DUNE_XT_COMMON_FMATRIX_27_HH
