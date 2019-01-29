// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2014 - 2018)
//   Rene Milk        (2015 - 2016, 2018)
//   Tobias Leibner   (2014, 2017 - 2018)

#ifndef DUNE_XT_LA_CONTAINER_MATRIX_INTERFACE_HH
#define DUNE_XT_LA_CONTAINER_MATRIX_INTERFACE_HH

#include <cmath>
#include <limits>
#include <iostream>
#include <type_traits>

#include <dune/common/ftraits.hh>

#include <dune/xt/common/crtp.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/unused.hh>

#include <dune/xt/la/type_traits.hh>

#include "container-interface.hh"
#include "pattern.hh"
#include "vector-interface.hh"

namespace Dune {
namespace XT {
namespace LA {
namespace internal {


template <class ScalarImp,
          class derived_imp,
          class BackendImp,
          Backends backend_imp,
          Backends vector_imp,
          bool is_sparse = false>
class MatrixTraitsBase
{
public:
  using ScalarType = typename Dune::FieldTraits<ScalarImp>::field_type;
  using RealType = typename Dune::FieldTraits<ScalarImp>::real_type;
  using MutexesType = std::vector<std::mutex>;
  using DataType = ScalarType;
  using BackendType = BackendImp;
  using derived_type = derived_imp;
  static constexpr Backends backend_type = backend_imp;
  static constexpr Backends vector_type = vector_imp;
  static constexpr bool sparse = is_sparse;
};


} // namespace internal


template <class TraitsImp, class ScalarImp = typename TraitsImp::ScalarType>
class MatrixInterface : public ContainerInterface<TraitsImp, ScalarImp>
{
  using BaseType = ContainerInterface<TraitsImp, ScalarImp>;

public:
  using typename BaseType::derived_type;
  using typename BaseType::Traits;
  using ScalarType = typename Traits::ScalarType;
  using RealType = typename Traits::RealType;
  using MutexesType = typename Traits::MutexesType;
  static constexpr Backends vector_type = Traits::vector_type;
  static constexpr bool sparse = Traits::sparse;
  static_assert(std::is_same<ScalarImp, typename Traits::ScalarType>::value, "");

  virtual ~MatrixInterface() = default;

  /// \name Have to be implemented by a derived class in addition to the ones required by ContainerInterface!
  /// \{

  inline size_t rows() const
  {
    CHECK_CRTP(this->as_imp().rows());
    return this->as_imp().rows();
  }

  inline size_t cols() const
  {
    CHECK_CRTP(this->as_imp().cols());
    return this->as_imp().cols();
  }

  template <class XX, class YY>
  inline void mv(const XX& xx, YY& yy) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().mv(xx, yy));
  }

  template <class XX, class YY>
  inline void mtv(const XX& xx, YY& yy) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().mtv(xx, yy));
  }

  inline void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    CHECK_AND_CALL_CRTP(this->as_imp().add_to_entry(ii, jj, value));
  }

  inline void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    CHECK_AND_CALL_CRTP(this->as_imp().set_entry(ii, jj, value));
  }

  inline ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    CHECK_CRTP(this->as_imp().get_entry(ii, jj));
    return this->as_imp().get_entry(ii, jj);
  }

  inline void clear_row(const size_t ii)
  {
    CHECK_AND_CALL_CRTP(this->as_imp().clear_row(ii));
  }

  inline void clear_col(const size_t jj)
  {
    CHECK_AND_CALL_CRTP(this->as_imp().clear_col(jj));
  }

  inline void unit_row(const size_t ii)
  {
    CHECK_AND_CALL_CRTP(this->as_imp().unit_row(ii));
  }

  inline void unit_col(const size_t jj)
  {
    CHECK_AND_CALL_CRTP(this->as_imp().unit_col(jj));
  }

  /**
   * \brief  Checks entries for inf or nan.
   * \return false if any entry is inf or nan, else true
   */
  inline bool valid() const
  {
    CHECK_CRTP(this->as_imp().valid());
    return this->as_imp().valid();
  }

  /// \}
  /// \name Provided by the interface for convenience.
  /// \note Those marked with vitual should be overriden by any devired class that can do better.
  /// \{

  using BaseType::operator*;

  template <class XX>
  typename XX::derived_type operator*(const VectorInterface<XX, ScalarType>& xx) const
  {
    typename XX::derived_type yy(rows());
    mv(xx.as_imp(), yy);
    return yy;
  }

  template <class MM>
  derived_type operator*(const MatrixInterface<MM, ScalarType>& other) const
  {
    return multiply(other);
  }

  virtual derived_type operator*(const derived_type& other) const
  {
    return multiply(other);
  }

  template <class MM>
  derived_type operator+(const MatrixInterface<MM, ScalarType>& other) const
  {
    return add(other);
  }

  virtual derived_type operator+(const derived_type& other) const
  {
    return add(other);
  }

  template <class MM>
  derived_type operator-(const MatrixInterface<MM, ScalarType>& other) const
  {
    return subtract(other);
  }

  virtual derived_type operator-(const derived_type& other) const
  {
    return subtract(other);
  }

  template <class MM>
  derived_type& operator+=(const MatrixInterface<MM, ScalarType>& other)
  {
    return add_assign(other);
  }

  virtual derived_type& operator+=(const derived_type& other)
  {
    return add_assign(other);
  }

  template <class MM>
  derived_type& operator-=(const MatrixInterface<MM, ScalarType>& other)
  {
    return subtract_assign(other);
  }

  virtual derived_type& operator-=(const derived_type& other)
  {
    return subtract_assign(other);
  }

  virtual RealType sup_norm() const
  {
    RealType ret = 0;
    for (size_t ii = 0; ii < rows(); ++ii)
      for (size_t jj = 0; jj < cols(); ++jj)
        ret = std::max(ret, std::abs(get_entry(ii, jj)));
    return ret;
  } // ... sup_norm(...)

  derived_type transposed() const
  {
    derived_type yy(rows(), cols(), 0.);
    for (size_t rr = 0; rr < rows(); ++rr)
      for (size_t cc = 0; cc < cols(); ++cc)
        yy.set_entry(rr, cc, get_entry(cc, rr));
    return yy;
  }

  /**
   * \brief Returns the number of entries in the sparsity pattern of the matrix.
   *
   * This is mainly useful for sparse matrices and returns rows() times cols() for dense matrices.
   *
   * \note Some implementations do not report the correct number here, so use and interpret only if you know what you
   * are doing!
   */
  virtual size_t non_zeros() const
  {
    return rows() * cols();
  }

  /**
   * \brief Computes the sparsity pattern of the matrix.
   *
   * This is mainly useful for sparse matrices and returns a full pattern for dense matrices
   *
   * \param prune If true, treats all entries smaller than eps as zero and does not include these indices in the
   * returned pattern
   */
  virtual SparsityPatternDefault pattern(const bool prune = false,
                                         const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
                                             Common::FloatCmp::DefaultEpsilon<ScalarType>::value()) const
  {
    SparsityPatternDefault ret(rows());
    const ScalarType zero(0);
    if (prune) {
      for (size_t ii = 0; ii < rows(); ++ii)
        for (size_t jj = 0; jj < cols(); ++jj)
          if (Common::FloatCmp::ne<Common::FloatCmp::Style::absolute>(get_entry(ii, jj), zero, eps))
            ret.insert(ii, jj);
    } else {
      for (size_t ii = 0; ii < rows(); ++ii)
        for (size_t jj = 0; jj < cols(); ++jj)
          ret.insert(ii, jj);
    }
    ret.sort();
    return ret;
  } // ... pattern(...)

  /**
   * \brief Returns a pruned variant of this matrix.
   *
   * This is mainly useful for sparse matrices and returns a matrix that should be very close to this matrix, except for
   * very small values, which are set to zero and the entries of which are removed from the sparsity pattern.
   *
   * \sa    pattern
   * \param eps Is forwarded to pattern(true, eps)
   */
  virtual derived_type pruned(const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
                                  Common::FloatCmp::DefaultEpsilon<ScalarType>::value()) const
  {
    const auto pruned_pattern = pattern(true, eps);
    derived_type ret(rows(), cols(), pruned_pattern);
    for (size_t ii = 0; ii < pruned_pattern.size(); ++ii)
      for (const size_t& jj : pruned_pattern.inner(ii))
        ret.set_entry(ii, jj, get_entry(ii, jj));
    return ret;
  } // ... pruned(...)

  virtual bool almost_equal(const derived_type& other,
                            const ScalarType epsilon = Common::FloatCmp::DefaultEpsilon<ScalarType>::value()) const
  {
    if (other.rows() != rows())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "rows(): " << rows() << "\n   "
                            << "other.rows(): " << other.rows());
    if (other.cols() != cols())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "cols(): " << cols() << "\n   "
                            << "other.cols(): " << other.cols());
    auto my_pattern = pattern();
    auto other_pattern = other.pattern();
    for (size_t ii = 0; ii < rows(); ++ii) {
      const auto my_cols = std::set<size_t>(my_pattern.inner(ii).begin(), my_pattern.inner(ii).end());
      const auto other_cols = std::set<size_t>(other_pattern.inner(ii).begin(), other_pattern.inner(ii).end());
      for (const auto& jj : my_cols)
        if (other_cols.count(jj) == 0) {
          if (Common::FloatCmp::ne(get_entry(ii, jj), ScalarType(0.), epsilon))
            return false;
        } else {
          if (Common::FloatCmp::ne(get_entry(ii, jj), other.get_entry(ii, jj), epsilon))
            return false;
        }
      for (const auto& jj : other_cols)
        if (my_cols.count(jj) == 0 && Common::FloatCmp::ne(other.get_entry(ii, jj), ScalarType(0.), epsilon))
          return false;
    }
    return true;
  } // ... almost_equal(...)
  /// \}

  template <class DuneDenseMatrixImp>
  void copy_to_densematrix(DuneDenseMatrixImp& ret) const
  {
    for (size_t ii = 0; ii < rows(); ++ii)
      for (size_t jj = 0; jj < cols(); ++jj)
        ret[ii][jj] = get_entry(ii, jj);
  }

  template <int ROWS, int COLS>
  explicit operator Dune::FieldMatrix<ScalarType, ROWS, COLS>() const
  {
    assert(ROWS == rows() && COLS == cols());
    Dune::FieldMatrix<ScalarType, ROWS, COLS> ret(ScalarType(0));
    CHECK_CRTP(this->as_imp().copy_to_densematrix(ret));
    this->as_imp().copy_to_densematrix(ret);
    return ret;
  }

  template <int ROWS, int COLS>
  explicit operator std::unique_ptr<Dune::FieldMatrix<ScalarType, ROWS, COLS>>() const
  {
    auto ret = XT::Common::make_unique<Dune::FieldMatrix<ScalarType, ROWS, COLS>>(ScalarType(0));
    CHECK_CRTP(this->as_imp().copy_to_densematrix(*ret));
    this->as_imp().copy_to_densematrix(*ret);
    return ret;
  }

  explicit operator Dune::DynamicMatrix<ScalarType>() const
  {
    Dune::DynamicMatrix<ScalarType> ret(rows(), cols(), ScalarType(0));
    CHECK_CRTP(this->as_imp().copy_to_densematrix(ret));
    this->as_imp().copy_to_densematrix(ret);
    return ret;
  }

protected:
  template <class MM>
  derived_type multiply(const MatrixInterface<MM, ScalarType>& other) const
  {
    if (other.rows() != cols())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match, "Dimensions of matrices to be multiplied do not match!");
    const auto new_pattern = multiplication_pattern(pattern(), other.pattern(), other.cols());
    derived_type yy(rows(), other.cols(), new_pattern);
    for (size_t rr = 0; rr < rows(); ++rr)
      for (size_t cc = 0; cc < other.cols(); ++cc)
        for (size_t kk = 0; kk < cols(); ++kk)
          yy.add_to_entry(rr, cc, get_entry(rr, kk) * other.get_entry(kk, cc));
    return yy;
  }

  template <class MM>
  derived_type add(const MatrixInterface<MM, ScalarType>& other) const
  {
    if (other.rows() != rows() || other.cols() != cols())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match, "Dimensions of matrices to be added do not match!");
    const auto new_pattern = pattern() + other.pattern();
    derived_type yy(rows(), other.cols(), new_pattern);
    for (size_t rr = 0; rr < rows(); ++rr)
      for (const auto& cc : new_pattern.inner(rr))
        yy.set_entry(rr, cc, get_entry(rr, cc) + other.get_entry(rr, cc));
    return yy;
  }

  template <class MM>
  derived_type subtract(const MatrixInterface<MM, ScalarType>& other) const
  {
    if (other.rows() != rows() || other.cols() != cols())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match, "Dimensions of matrices to be subtracted do not match!");
    const auto new_pattern = pattern() + other.pattern();
    derived_type yy(rows(), other.cols(), new_pattern);
    for (size_t rr = 0; rr < rows(); ++rr)
      for (const auto& cc : new_pattern.inner(rr))
        yy.set_entry(rr, cc, get_entry(rr, cc) - other.get_entry(rr, cc));
    return yy;
  }

  template <class MM>
  derived_type& add_assign(const MatrixInterface<MM, ScalarType>& other)
  {
    if (other.rows() != rows() || other.cols() != cols())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match, "Dimensions of matrices to be added do not match!");
    const auto this_pattern = pattern();
    auto new_pattern = this_pattern + other.pattern();
    if (new_pattern != this_pattern)
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "The matrix to be added contains entries that are not in this' pattern!");
    for (size_t rr = 0; rr < rows(); ++rr)
      for (const auto& cc : this_pattern.inner(rr))
        add_to_entry(rr, cc, other.get_entry(rr, cc));
    return this->as_imp();
  }

  template <class MM>
  derived_type& subtract_assign(const MatrixInterface<MM, ScalarType>& other)
  {
    if (other.rows() != rows() || other.cols() != cols())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match, "Dimensions of matrices to be added do not match!");
    const auto this_pattern = pattern();
    auto new_pattern = this_pattern + other.pattern();
    if (new_pattern != this_pattern)
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "The matrix to be subtracted contains entries that are not in this' pattern!");
    for (size_t rr = 0; rr < rows(); ++rr)
      for (const auto& cc : this_pattern.inner(rr))
        add_to_entry(rr, cc, -other.get_entry(rr, cc));
    return this->as_imp();
  }

private:
  template <class T, class S>
  friend std::ostream& operator<<(std::ostream& /*out*/, const MatrixInterface<T, S>& /*matrix*/);
}; // class MatrixInterface


template <class T, class S>
std::ostream& operator<<(std::ostream& out, const MatrixInterface<T, S>& matrix)
{
  out << "[";
  const size_t rows = matrix.rows();
  const size_t cols = matrix.cols();
  if (rows > 0 && cols > 0) {
    for (size_t ii = 0; ii < rows; ++ii) {
      if (ii > 0)
        out << "\n ";
      out << "[" << matrix.get_entry(ii, 0);
      for (size_t jj = 1; jj < cols; ++jj)
        out << " " << matrix.get_entry(ii, jj);
      out << "]";
      if (rows > 1 && ii < (rows - 1))
        out << ",";
    }
    out << "]";
  } else
    out << "[ ]]";
  return out;
} // ... operator<<(...)


template <class MatrixType>
typename std::enable_if<is_matrix<MatrixType>::value, MatrixType>::type
create(const std::initializer_list<std::initializer_list<typename Common::MatrixAbstraction<MatrixType>::ScalarType>>&
           initializer_list,
       const SparsityPatternDefault& pattern = SparsityPatternDefault())
{
  const size_t num_rows = initializer_list.size();
  const size_t num_cols = initializer_list.begin()->size();
  MatrixType ret(num_rows, num_cols, pattern.size() == 0 ? dense_pattern(num_rows, num_cols) : pattern);
  size_t ii = 0;
  for (const auto& row : initializer_list) {
    size_t jj = 0;
    for (const auto& entry : row) {
      if (entry != 0.)
        ret.set_entry(ii, jj, entry);
      ++jj;
    } // jj
    ++ii;
  } // ii
  return ret;
}


template <class MatrixType>
typename std::enable_if<Common::is_matrix<MatrixType>::value && !is_matrix<MatrixType>::value, MatrixType>::type
create(const std::initializer_list<std::initializer_list<typename Common::MatrixAbstraction<MatrixType>::ScalarType>>&
           initializer_list,
       const SparsityPatternDefault& /*pattern*/ = SparsityPatternDefault())
{
  using Mat = Common::MatrixAbstraction<MatrixType>;
  auto ret = Mat::create(initializer_list.size(), initializer_list.begin()->size());
  size_t ii = 0;
  for (const auto& row : initializer_list) {
    size_t jj = 0;
    for (const auto& entry : row) {
      Mat::set_entry(ret, ii, jj, entry);
      ++jj;
    } // jj
    ++ii;
  } // ii
  return ret;
}


namespace internal {


template <class MatrixImp>
struct MatrixAbstractionBase
{
  static const bool is_matrix = LA::is_matrix<MatrixImp>::value;

  static const bool has_static_size = false;

  static const size_t static_rows = std::numeric_limits<size_t>::max();

  static const size_t static_cols = std::numeric_limits<size_t>::max();

  static const constexpr Common::StorageLayout storage_layout = Common::StorageLayout::other;

  static const bool has_ostream = true;

  using MatrixType = typename std::conditional<is_matrix, MatrixImp, void>::type;
  using ScalarType = typename std::conditional<is_matrix, typename MatrixImp::ScalarType, void>::type;
  using RealType = typename std::conditional<is_matrix, typename MatrixImp::RealType, void>::type;
  using S = ScalarType;
  using R = RealType;

  template <class SparsityPatternType = XT::Common::FullPattern>
  static inline typename std::enable_if<is_matrix, MatrixType>::type
  create(const size_t rows,
         const size_t cols,
         const ScalarType& val = XT::Common::suitable_default<ScalarType>::value(),
         const SparsityPatternType& pattern = SparsityPatternType())
  {
    SparsityPatternDefault actual_pattern = get_actual_pattern(pattern, rows, cols);
    MatrixType ret(rows, cols, actual_pattern);
    fill_matrix(ret, val, actual_pattern);
    return ret;
  }

  template <class SparsityPatternType = XT::Common::FullPattern>
  static inline typename std::enable_if<is_matrix, std::unique_ptr<MatrixType>>::type
  make_unique(const size_t rows,
              const size_t cols,
              const ScalarType& val = XT::Common::suitable_default<ScalarType>::value(),
              const SparsityPatternType& pattern = SparsityPatternType())
  {
    SparsityPatternDefault actual_pattern = get_actual_pattern(pattern, rows, cols);
    auto ret = std::make_unique<MatrixType>(rows, cols, actual_pattern);
    fill_matrix(*ret, val, actual_pattern);
    return ret;
  }

  static inline typename std::enable_if<is_matrix, size_t>::type rows(const MatrixType& mat)
  {
    return mat.rows();
  }

  static inline typename std::enable_if<is_matrix, size_t>::type cols(const MatrixType& mat)
  {
    return mat.cols();
  }

  static void set_entry(MatrixType& mat, const size_t row, const size_t col, const ScalarType& val)
  {
    mat.set_entry(row, col, val);
  }

  static inline typename std::enable_if<is_matrix, ScalarType>::type
  get_entry(const MatrixType& mat, const size_t row, const size_t col)
  {
    return mat.get_entry(row, col);
  }

  static void add_to_entry(MatrixType& mat, const size_t row, const size_t col, const ScalarType& val)
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

private:
  template <class SparsityPatternType>
  static SparsityPatternDefault
  get_actual_pattern(SparsityPatternType&& pattern, size_t DXTC_DEBUG_ONLY(rows), size_t /*cols*/)
  {
    assert(rows == pattern.size());
    return std::forward<SparsityPatternType>(pattern);
  }

  static SparsityPatternDefault get_actual_pattern(const XT::Common::FullPattern& /*pattern*/, size_t rows, size_t cols)
  {
    return dense_pattern(rows, cols);
  }

  static void fill_matrix(MatrixType& matrix, const ScalarType& val, const SparsityPatternDefault& pattern)
  {
    for (size_t ii = 0; ii < rows(matrix); ++ii)
      for (const auto& jj : pattern.inner(ii))
        matrix.set_entry(ii, jj, val);
  }
}; // struct MatrixAbstractionBase


} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_MATRIX_INTERFACE_HH
