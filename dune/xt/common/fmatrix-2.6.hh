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

#ifndef DUNE_XT_COMMON_FMATRIX_26_HH
#define DUNE_XT_COMMON_FMATRIX_26_HH

#include <initializer_list>

#include <dune/xt/common/debug.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/unused.hh>

namespace Dune {
namespace XT {
namespace Common {


/**
 * \todo We need to implement all operators from the base which return the base, to rather return ourselfes!
 */
template <class K, int ROWS, int COLS>
class FieldMatrix : public Dune::FieldMatrix<K, ROWS, COLS>
{
  typedef Dune::FieldMatrix<K, ROWS, COLS> BaseType;
  typedef FieldMatrix<K, ROWS, COLS> ThisType;

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

  FieldMatrix(DXTC_DEBUG_ONLY const size_t rr,
              DXTC_DEBUG_ONLY const size_t cc,
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

  //! This op is not redundant
  ThisType operator*(const K& scal) const
  {
    ThisType ret(*this);
    ret *= scal;
    return ret;
  }

  field_type determinant() const;

  template <class Func>
  void luDecomposition(ThisType& A, Func func) const;

  void invert();

  template <class V, class W>
  void solve(V& x, const W& b) const;

private:
  // copy from dune/common/densematrix.hh, we have to copy it as it is a private member of Dune::DenseMatrix
  struct ElimPivot
  {
    ElimPivot(std::vector<size_type>& pivot)
      : pivot_(pivot)
    {
      typedef typename std::vector<size_type>::size_type size_type;
      for (size_type i = 0; i < pivot_.size(); ++i)
        pivot_[i] = i;
    }

    void swap(int i, int j)
    {
      pivot_[i] = j;
    }

    template <typename T>
    void operator()(const T&, int, int)
    {}

    std::vector<size_type>& pivot_;
  }; // struct ElimPivot

  template <typename V>
  struct Elim
  {
    Elim(V& rhs)
      : rhs_(&rhs)
    {}

    void swap(int i, int j)
    {
      std::swap((*rhs_)[i], (*rhs_)[j]);
    }

    void operator()(const typename V::field_type& factor, int k, int i)
    {
      (*rhs_)[k] -= factor * (*rhs_)[i];
    }

    V* rhs_;
  };

  struct ElimDet
  {
    ElimDet(field_type& sign)
      : sign_(sign)
    {
      sign_ = 1;
    }

    void swap(int, int)
    {
      sign_ *= -1;
    }

    void operator()(const field_type&, int, int) {}

    field_type& sign_;
  };
}; // class FieldMatrix<...>

// Direct copy of the luDecomposition function in dune/common/densematrix.hh
// The only (functional) change is that this version always performs pivotization (the version in dune-common only
// performs pivotization if the diagonal entry is below a certain threshold)
// See dune/xt/la/test/matrixinverter_for_real_matrix_from_3d_pointsource.tpl for an example where the dune-common
// version fails due to stability issues.
// TODO: Fixed in dune-common master (see MR !449 in dune-common's gitlab), remove this copy once we depend on a
// suitable version of dune-common (probably 2.7).
template <class K, int ROWS, int COLS>
template <typename Func>
inline void FieldMatrix<K, ROWS, COLS>::luDecomposition(FieldMatrix<K, ROWS, COLS>& A, Func func) const
{
  typedef typename FieldTraits<value_type>::real_type real_type;
  real_type norm = A.infinity_norm_real(); // for relative thresholds
  real_type singthres =
      std::max(FMatrixPrecision<real_type>::absolute_limit(), norm * FMatrixPrecision<real_type>::singular_limit());

  // LU decomposition of A in A
  for (size_type i = 0; i < ROWS; i++) // loop over all rows
  {
    typename FieldTraits<value_type>::real_type pivmax = fvmeta::absreal(A[i][i]);

    // compute maximum of column
    size_type imax = i;
    typename FieldTraits<value_type>::real_type abs(0.0);
    for (size_type k = i + 1; k < ROWS; k++)
      if ((abs = fvmeta::absreal(A[k][i])) > pivmax) {
        pivmax = abs;
        imax = k;
      }
    // swap rows
    if (imax != i) {
      for (size_type j = 0; j < ROWS; j++)
        std::swap(A[i][j], A[imax][j]);
      assert(imax < std::numeric_limits<int>::max() && i < std::numeric_limits<int>::max());
      func.swap(static_cast<int>(i), static_cast<int>(imax)); // swap the pivot or rhs
    }

    // singular ?
    if (pivmax < singthres)
      DUNE_THROW(FMatrixError, "matrix is singular");

    // eliminate
    for (size_type k = i + 1; k < ROWS; k++) {
      field_type factor = A[k][i] / A[i][i];
      A[k][i] = factor;
      for (size_type j = i + 1; j < ROWS; j++)
        A[k][j] -= factor * A[i][j];
      assert(k < std::numeric_limits<int>::max() && i < std::numeric_limits<int>::max());
      func(factor, static_cast<int>(k), static_cast<int>(i));
    }
  }
}

// Direct copy of the invert function in dune/common/densematrix.hh
// The only (functional) change is the replacement of the luDecomposition of DenseMatrix by our own version.
// TODO: Fixed in dune-common master (see MR !449 in dune-common's gitlab), remove this copy once we depend on a
// suitable version of dune-common (probably 2.7).
template <class K, int ROWS, int COLS>
inline void FieldMatrix<K, ROWS, COLS>::invert()
{
  // never mind those ifs, because they get optimized away
  if (ROWS != COLS)
    DUNE_THROW(Dune::FMatrixError, "Can't invert a " << ROWS << "x" << COLS << " matrix!");
  if (ROWS <= 3) {
    BaseType::invert();
  } else {
    auto A = *this;
    std::vector<size_type> pivot(ROWS);
    this->luDecomposition(A, ElimPivot(pivot));
    auto& L = A;
    auto& U = A;

    // initialize inverse
    *this = field_type();

    for (size_type i = 0; i < ROWS; ++i)
      (*this)[i][i] = 1;

    // L Y = I; multiple right hand sides
    for (size_type i = 0; i < ROWS; i++)
      for (size_type j = 0; j < i; j++)
        for (size_type k = 0; k < ROWS; k++)
          (*this)[i][k] -= L[i][j] * (*this)[j][k];

    // U A^{-1} = Y
    for (size_type i = ROWS; i > 0;) {
      --i;
      for (size_type k = 0; k < ROWS; k++) {
        for (size_type j = i + 1; j < ROWS; j++)
          (*this)[i][k] -= U[i][j] * (*this)[j][k];
        (*this)[i][k] /= U[i][i];
      }
    }

    for (size_type i = ROWS; i > 0;) {
      --i;
      if (i != pivot[i])
        for (size_type j = 0; j < ROWS; ++j)
          std::swap((*this)[j][pivot[i]], (*this)[j][i]);
    }
  }
}

// Direct copy of the determinant function in dune/common/densematrix.hh
// The only (functional) change is the replacement of the luDecomposition of DenseMatrix by our own version.
// TODO: Fixed in dune-common master (see MR !449 in dune-common's gitlab), remove this copy once we depend on a
// suitable version of dune-common (probably 2.7).
template <class K, int ROWS, int COLS>
inline typename FieldMatrix<K, ROWS, COLS>::field_type FieldMatrix<K, ROWS, COLS>::determinant() const
{
  // never mind those ifs, because they get optimized away
  if (ROWS != COLS)
    DUNE_THROW(FMatrixError, "There is no determinant for a " << ROWS << "x" << COLS << " matrix!");

  if (ROWS == 1)
    return (*this)[0][0];

  if (ROWS == 2)
    return (*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0];

  if (ROWS == 3) {
    // code generated by maple
    field_type t4 = (*this)[0][0] * (*this)[1][1];
    field_type t6 = (*this)[0][0] * (*this)[1][2];
    field_type t8 = (*this)[0][1] * (*this)[1][0];
    field_type t10 = (*this)[0][2] * (*this)[1][0];
    field_type t12 = (*this)[0][1] * (*this)[2][0];
    field_type t14 = (*this)[0][2] * (*this)[2][0];

    return (t4 * (*this)[2][2] - t6 * (*this)[2][1] - t8 * (*this)[2][2] + t10 * (*this)[2][1] + t12 * (*this)[1][2]
            - t14 * (*this)[1][1]);
  }

  auto A = *this;
  field_type det;
  try {
    this->luDecomposition(A, ElimDet(det));
  } catch (FMatrixError&) {
    return 0;
  }
  for (size_type i = 0; i < ROWS; ++i)
    det *= A[i][i];
  return det;
}


// Direct copy of the solve function in dune/common/densematrix.hh
// The only (functional) change is the replacement of the luDecomposition of DenseMatrix by our own version.
// TODO: Fixed in dune-common master (see MR !449 in dune-common's gitlab), remove this copy once we depend on a
// suitable version of dune-common (probably 2.7).
template <class K, int ROWS, int COLS>
template <class V, class W>
inline void FieldMatrix<K, ROWS, COLS>::solve(V& x, const W& b) const
{
  // never mind those ifs, because they get optimized away
  if (ROWS != COLS)
    DUNE_THROW(FMatrixError, "Can't solve for a " << ROWS << "x" << COLS << " matrix!");

  if (ROWS == 1) {

#ifdef DUNE_FMatrix_WITH_CHECKING
    if (fvmeta::absreal((*this)[0][0]) < FMatrixPrecision<>::absolute_limit())
      DUNE_THROW(FMatrixError, "matrix is singular");
#endif
    x[0] = b[0] / (*this)[0][0];

  } else if (ROWS == 2) {

    field_type detinv = (*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0];
#ifdef DUNE_FMatrix_WITH_CHECKING
    if (fvmeta::absreal(detinv) < FMatrixPrecision<>::absolute_limit())
      DUNE_THROW(FMatrixError, "matrix is singular");
#endif
    detinv = 1.0 / detinv;

    x[0] = detinv * ((*this)[1][1] * b[0] - (*this)[0][1] * b[1]);
    x[1] = detinv * ((*this)[0][0] * b[1] - (*this)[1][0] * b[0]);

  } else if (ROWS == 3) {

    field_type d = determinant();
#ifdef DUNE_FMatrix_WITH_CHECKING
    if (fvmeta::absreal(d) < FMatrixPrecision<>::absolute_limit())
      DUNE_THROW(FMatrixError, "matrix is singular");
#endif

    x[0] = (b[0] * (*this)[1][1] * (*this)[2][2] - b[0] * (*this)[2][1] * (*this)[1][2]
            - b[1] * (*this)[0][1] * (*this)[2][2] + b[1] * (*this)[2][1] * (*this)[0][2]
            + b[2] * (*this)[0][1] * (*this)[1][2] - b[2] * (*this)[1][1] * (*this)[0][2])
           / d;

    x[1] = ((*this)[0][0] * b[1] * (*this)[2][2] - (*this)[0][0] * b[2] * (*this)[1][2]
            - (*this)[1][0] * b[0] * (*this)[2][2] + (*this)[1][0] * b[2] * (*this)[0][2]
            + (*this)[2][0] * b[0] * (*this)[1][2] - (*this)[2][0] * b[1] * (*this)[0][2])
           / d;

    x[2] = ((*this)[0][0] * (*this)[1][1] * b[2] - (*this)[0][0] * (*this)[2][1] * b[1]
            - (*this)[1][0] * (*this)[0][1] * b[2] + (*this)[1][0] * (*this)[2][1] * b[0]
            + (*this)[2][0] * (*this)[0][1] * b[1] - (*this)[2][0] * (*this)[1][1] * b[0])
           / d;

  } else {

    V& rhs = x; // use x to store rhs
    rhs = b; // copy data
    Elim<V> elim(rhs);
    auto A = *this;

    this->luDecomposition(A, elim);

    // backsolve
    for (int i = ROWS - 1; i >= 0; i--) {
      for (size_type j = i + 1; j < ROWS; j++)
        rhs[i] -= A[i][j] * x[j];
      x[i] = rhs[i] / A[i][i];
    }
  }
}


/**
 * \todo We need to implement all operators from the base which return the base, to rather return ourselfes!
 */
template <class K>
class FieldMatrix<K, 1, 1> : public Dune::FieldMatrix<K, 1, 1>
{
  static constexpr int ROWS = 1;
  static constexpr int COLS = 1;
  typedef Dune::FieldMatrix<K, ROWS, COLS> BaseType;
  typedef FieldMatrix<K, ROWS, COLS> ThisType;

public:
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

  FieldMatrix(DXTC_DEBUG_ONLY const size_t rr,
              DXTC_DEBUG_ONLY const size_t cc,
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

  ThisType operator*(const K& scal) const
  {
    ThisType ret(*this);
    ret *= scal;
    return ret;
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

  BlockedFieldMatrix(DXTC_DEBUG_ONLY const size_t rows, DXTC_DEBUG_ONLY const size_t cols, const K& val = K(0.))
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
  typedef typename Dune::FieldTraits<K>::field_type ScalarType;
  typedef typename Dune::FieldTraits<K>::real_type RealType;
  typedef ScalarType S;
  typedef RealType R;
  template <size_t rows = N, size_t cols = M, class FieldType = K>
  using MatrixTypeTemplate = Dune::XT::Common::FieldMatrix<FieldType, rows, cols>;

  static constexpr bool is_matrix = true;

  static constexpr bool has_static_size = true;

  static constexpr size_t static_rows = N;

  static constexpr size_t static_cols = M;

  static constexpr StorageLayout storage_layout = StorageLayout::dense_row_major;

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
  typedef typename Dune::FieldTraits<K>::field_type ScalarType;
  typedef typename Dune::FieldTraits<K>::real_type RealType;
  typedef ScalarType S;
  typedef RealType R;

  static constexpr bool is_matrix = true;

  static constexpr bool has_static_size = true;

  static constexpr size_t static_rows = MatrixType::num_rows;

  static constexpr size_t static_cols = MatrixType::num_cols;

  template <size_t rows = static_rows, size_t cols = static_cols, class FieldType = K>
  using MatrixTypeTemplate = Dune::XT::Common::BlockedFieldMatrix<FieldType, rows / block_rows, block_rows, block_cols>;

  static constexpr StorageLayout storage_layout = StorageLayout::other;

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
  static constexpr size_t rows = MatrixAbstraction<M>::static_rows;
  static constexpr size_t cols = MatrixAbstraction<M>::static_cols;
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
  static constexpr size_t rows = MatrixAbstraction<M>::static_rows;
  static constexpr size_t cols = MatrixAbstraction<M>::static_cols;
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
typename std::enable_if<
    !std::is_same<L, R>::value,
    Dune::XT::Common::FieldMatrix<typename PromotionTraits<L, R>::PromotedType, L_ROWS, R_COLS>>::type
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


} // namespace Dune

#endif // DUNE_XT_COMMON_FMATRIX_26_HH
