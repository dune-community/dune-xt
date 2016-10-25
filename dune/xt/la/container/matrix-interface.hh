// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2014 - 2016)
//   Rene Milk        (2015)
//   Tobias Leibner   (2014)

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
#include <dune/xt/common/type_traits.hh>

#include "container-interface.hh"
#include "pattern.hh"
#include "vector-interface.hh"

namespace Dune {
namespace XT {
namespace LA {


template <class Traits, class ScalarImp = typename Traits::ScalarType>
class MatrixInterface : public ContainerInterface<Traits, ScalarImp>
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Dune::FieldTraits<ScalarImp>::field_type ScalarType;
  typedef typename Dune::FieldTraits<ScalarImp>::real_type RealType;
  static const constexpr Backends vector_type = Traits::vector_type;
  static_assert(std::is_same<ScalarType, typename Traits::ScalarType>::value, "");

  virtual ~MatrixInterface()
  {
  }

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

  template <class XX>
  typename XX::derived_type operator*(const VectorInterface<XX, ScalarType>& xx) const
  {
    typename XX::derived_type yy(rows());
    mv(xx.as_imp(xx), yy);
    return yy;
  }

  virtual RealType sup_norm() const
  {
    RealType ret = 0;
    for (size_t ii = 0; ii < rows(); ++ii)
      for (size_t jj = 0; jj < cols(); ++jj)
        ret = std::max(ret, std::abs(get_entry(ii, jj)));
    return ret;
  } // ... sup_norm(...)

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

namespace internal {

template <class M>
struct is_matrix_helper
{
  DXTC_has_typedef_initialize_once(Traits);
  DXTC_has_typedef_initialize_once(ScalarType);

  static const bool is_candidate = DXTC_has_typedef(Traits)<M>::value && DXTC_has_typedef(ScalarType)<M>::value;
}; // class is_matrix_helper

} // namespace internal

template <class M, bool candidate = internal::is_matrix_helper<M>::is_candidate>
struct is_matrix : public std::is_base_of<MatrixInterface<typename M::Traits, typename M::ScalarType>, M>
{
};

template <class M>
struct is_matrix<M, false> : public std::false_type
{
};

namespace internal {

template <class MatrixImp>
struct MatrixAbstractionBase
{
  static const bool is_matrix = LA::is_matrix<MatrixImp>::value;

  static const bool has_static_size = false;

  static const size_t static_rows = std::numeric_limits<size_t>::max();

  static const size_t static_cols = std::numeric_limits<size_t>::max();

  typedef typename std::conditional<is_matrix, MatrixImp, void>::type MatrixType;
  typedef typename std::conditional<is_matrix, typename MatrixImp::ScalarType, void>::type ScalarType;
  typedef ScalarType S;

  static inline typename std::enable_if<is_matrix, MatrixType>::type create(const size_t rows, const size_t cols)
  {
    return MatrixType(rows, cols);
  }

  static inline typename std::enable_if<is_matrix, MatrixType>::type
  create(const size_t rows, const size_t cols, const ScalarType& val)
  {
    return MatrixType(rows, cols, val);
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
}; // struct MatrixAbstractionBase

} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_MATRIX_INTERFACE_HH
