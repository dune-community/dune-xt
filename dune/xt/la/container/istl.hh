// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr     (2014)
//   Barbara Verfürth (2015)
//   Felix Schindler  (2013 - 2017)
//   Rene Milk        (2014 - 2016)
//   Tobias Leibner   (2014, 2016)

#ifndef DUNE_XT_LA_CONTAINER_ISTL_HH
#define DUNE_XT_LA_CONTAINER_ISTL_HH

#include <vector>
#include <initializer_list>
#include <complex>
#include <mutex>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/unused.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#endif

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/timings.hh>
#include <dune/xt/common/math.hh>

#include "interfaces.hh"
#include "pattern.hh"

namespace Dune {
namespace XT {
namespace LA {

// forward
template <class ScalarImp>
class IstlDenseVector;

template <class ScalarImp>
class IstlRowMajorSparseMatrix;

#if HAVE_DUNE_ISTL

namespace internal {

/**
 * \brief Traits for IstlDenseVector.
 */
template <class ScalarImp>
class IstlDenseVectorTraits
{
public:
  typedef typename Dune::FieldTraits<ScalarImp>::field_type ScalarType;
  typedef typename Dune::FieldTraits<ScalarImp>::real_type RealType;
  typedef ScalarType DataType;
  typedef IstlDenseVector<ScalarImp> derived_type;
  typedef BlockVector<FieldVector<ScalarType, 1>> BackendType;
  static const constexpr Backends backend_type = Backends::istl_dense;
  static const constexpr Backends dense_matrix_type = Backends::none;
  static const constexpr Backends sparse_matrix_type = Backends::istl_sparse;
}; // class IstlDenseVectorTraits

/**
 * \brief Traits for IstlRowMajorSparseMatrix.
 */
template <class ScalarImp>
class IstlRowMajorSparseMatrixTraits
{
public:
  typedef typename Dune::FieldTraits<ScalarImp>::field_type ScalarType;
  typedef typename Dune::FieldTraits<ScalarImp>::real_type RealType;
  typedef IstlRowMajorSparseMatrix<ScalarType> derived_type;
  typedef BCRSMatrix<FieldMatrix<ScalarType, 1, 1>> BackendType;
  static const constexpr Backends backend_type = Backends::istl_sparse;
  static const constexpr Backends vector_type = Backends::istl_dense;
}; // class RowMajorSparseMatrixTraits

} // namespace internal

/**
 *  \brief A dense vector implementation of VectorInterface using the Dune::BlockVector from dune-istl.
 */
template <class ScalarImp = double>
class IstlDenseVector : public VectorInterface<internal::IstlDenseVectorTraits<ScalarImp>, ScalarImp>,
                        public ProvidesBackend<internal::IstlDenseVectorTraits<ScalarImp>>,
                        public ProvidesDataAccess<internal::IstlDenseVectorTraits<ScalarImp>>
{
  typedef IstlDenseVector<ScalarImp> ThisType;
  typedef VectorInterface<internal::IstlDenseVectorTraits<ScalarImp>, ScalarImp> VectorInterfaceType;

public:
  typedef internal::IstlDenseVectorTraits<ScalarImp> Traits;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::RealType RealType;
  typedef typename Traits::DataType DataType;
  typedef typename Traits::BackendType BackendType;

  explicit IstlDenseVector(const size_t ss = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(ss))
  {
    backend_->operator=(value);
  }

  explicit IstlDenseVector(const std::vector<ScalarType>& other)
    : backend_(new BackendType(other.size()))
  {
    for (size_t ii = 0; ii < other.size(); ++ii)
      backend_->operator[](ii)[0] = other[ii];
  }

  explicit IstlDenseVector(const std::initializer_list<ScalarType>& other)
    : backend_(new BackendType(other.size()))
  {
    size_t ii = 0;
    for (auto element : other) {
      backend_->operator[](ii)[0] = element;
      ++ii;
    }
  } // IstlDenseVector(...)

  IstlDenseVector(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  explicit IstlDenseVector(const BackendType& other,
                           const bool /*prune*/ = false,
                           const ScalarType /*eps*/ = Common::FloatCmp::DefaultEpsilon<ScalarType>::value())
    : backend_(new BackendType(other))
  {
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  explicit IstlDenseVector(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  explicit IstlDenseVector(std::shared_ptr<BackendType> backend_ptr)
    : backend_(backend_ptr)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    backend_ = other.backend_;
    return *this;
  }

  /**
   *  \note Does a deep copy.
   */
  ThisType& operator=(const BackendType& other)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    backend_ = std::make_shared<BackendType>(other);
    return *this;
  }

  /// \name Required by the ProvidesBackend interface.
  /// \{

  ThisType& operator=(const ScalarType& value)
  {
    ensure_uniqueness();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    for (auto& element : *this)
      element = value;
    return *this;
  }

  BackendType& backend()
  {
    ensure_uniqueness();
    return *backend_;
  }

  const BackendType& backend() const
  {
    return *backend_;
  }

  /// \}
  /// \name Required by ProvidesDataAccess.
  /// \{

  /** \attention This makes only sense for scalar data types, not for complex! **/
  DataType* data()
  {
    return &(backend()[0][0]);
  }

  size_t data_size() const
  {
    return size();
  }

  /// \}
  /// \name Required by ContainerInterface.
  /// \{

  ThisType copy() const
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    return ThisType(*backend_);
  }

  void scal(const ScalarType& alpha)
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    backend_ref *= alpha;
  }

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    if (xx.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of x (" << xx.size() << ") does not match the size of this (" << size() << ")!");
    backend_ref.axpy(alpha, xx.backend());
  }

  bool has_equal_shape(const ThisType& other) const
  {
    return size() == other.size();
  }

  /// \}
  /// \name Required by VectorInterface.
  /// \{

  inline size_t size() const
  {
    // as long as we have scalar blocks of size 1 here,
    // using backend's size would give a severe performance hit
    // since that iterates over the entire vector summing up 1's
    return backend_->N();
  }

  void add_to_entry(const size_t ii, const ScalarType& value)
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    assert(ii < size());
    backend_ref[ii][0] += value;
  }

  void set_entry(const size_t ii, const ScalarType& value)
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    assert(ii < size());
    backend_ref[ii][0] = value;
  }

  ScalarType get_entry(const size_t ii) const
  {
    assert(ii < size());
    return backend_->operator[](ii)[0];
  }

protected:
  inline ScalarType& get_entry_ref(const size_t ii)
  {
    return backend_->operator[](ii)[0];
  }

  inline const ScalarType& get_entry_ref(const size_t ii) const
  {
    return backend_->operator[](ii)[0];
  }

public:
  inline ScalarType& operator[](const size_t ii)
  {
    return backend()[ii][0];
  }

  inline const ScalarType& operator[](const size_t ii) const
  {
    return backend()[ii][0];
  }

  /// \}
  /// \name These methods override default implementations from VectorInterface..
  /// \{

  virtual ScalarType dot(const ThisType& other) const override final
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    return backend().dot(other.backend());
  } // ... dot(...)

  virtual RealType l1_norm() const override final
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    return backend().one_norm();
  }

  virtual RealType l2_norm() const override final
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    return backend().two_norm();
  }

  virtual RealType sup_norm() const override final
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    return backend().infinity_norm();
  }

  virtual void iadd(const ThisType& other) override final
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    backend_ref += other.backend();
  } // ... iadd(...)

  virtual void isub(const ThisType& other) override final
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    backend_ref -= other.backend();
  } // ... isub(...)

  /// \}

protected:
  /**
   * \see ContainerInterface
   */
  inline void ensure_uniqueness()
  {
    if (!backend_.unique()) {
      std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
      backend_ = std::make_shared<BackendType>(*backend_);
    }
  } // ... ensure_uniqueness(...)

private:
  friend class VectorInterface<internal::IstlDenseVectorTraits<ScalarType>, ScalarType>;
  friend class IstlRowMajorSparseMatrix<ScalarType>;

  std::shared_ptr<BackendType> backend_;
  mutable std::mutex mutex_;
}; // class IstlDenseVector

/**
 * \brief A sparse matrix implementation of the MatrixInterface using the Dune::BCRSMatrix from dune-istl.
 */
template <class ScalarImp = double>
class IstlRowMajorSparseMatrix : public MatrixInterface<internal::IstlRowMajorSparseMatrixTraits<ScalarImp>, ScalarImp>,
                                 public ProvidesBackend<internal::IstlRowMajorSparseMatrixTraits<ScalarImp>>
{
  typedef IstlRowMajorSparseMatrix<ScalarImp> ThisType;

public:
  typedef internal::IstlRowMajorSparseMatrixTraits<ScalarImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::RealType RealType;

  static std::string static_id()
  {
    return "xt.la.container.istl.istlrowmajorsparsematrix";
  }

  /**
   * \brief This is the constructor of interest which creates a sparse matrix.
   */
  IstlRowMajorSparseMatrix(const size_t rr, const size_t cc, const SparsityPatternDefault& patt)
  {
    if (patt.size() != rr)
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of the pattern (" << patt.size() << ") does not match the number of rows of this (" << rows()
                                             << ")!");
    build_sparse_matrix(rr, cc, patt);
    backend_->operator*=(ScalarType(0));
  } // ... IstlRowMajorSparseMatrix(...)

  explicit IstlRowMajorSparseMatrix(const size_t rr = 0, const size_t cc = 0)
    : backend_(new BackendType(rr, cc, BackendType::row_wise))
  {
  }

  IstlRowMajorSparseMatrix(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  explicit IstlRowMajorSparseMatrix(const BackendType& mat,
                                    const bool prune = false,
                                    const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
                                        Common::FloatCmp::DefaultEpsilon<ScalarType>::value())
  {
    if (prune) {
      const auto pruned_pattern = pruned_pattern_from_backend(mat, eps);
      build_sparse_matrix(mat.N(), mat.M(), pruned_pattern);
      for (size_t ii = 0; ii < pruned_pattern.size(); ++ii) {
        const auto& row_indices = pruned_pattern.inner(ii);
        if (row_indices.size() > 0) {
          const auto& mat_row = mat[ii];
          auto& backend_row = backend_->operator[](ii);
          for (const auto& jj : row_indices)
            backend_row[jj][0][0] = mat_row[jj][0][0];
        }
      }
    } else
      backend_ = std::shared_ptr<BackendType>(new BackendType(mat));
  } // IstlRowMajorSparseMatrix(...)

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  explicit IstlRowMajorSparseMatrix(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  explicit IstlRowMajorSparseMatrix(std::shared_ptr<BackendType> backend_ptr)
    : backend_(backend_ptr)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    backend_ = other.backend_;
    return *this;
  } // ... operator=(...)

  /**
   *  \note Does a deep copy.
   */
  ThisType& operator=(const BackendType& other)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    backend_ = std::make_shared<BackendType>(other);
    return *this;
  } // ... operator=(...)

  /// \name Required by the ProvidesBackend interface.
  /// \{

  BackendType& backend()
  {
    ensure_uniqueness();
    return *backend_;
  }

  const BackendType& backend() const
  {
    return *backend_;
  }

  /// \}
  /// \name Required by ContainerInterface.
  /// \{

  ThisType copy() const
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    return ThisType(*backend_);
  }

  void scal(const ScalarType& alpha)
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    backend_ref *= alpha;
  }

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    if (!has_equal_shape(xx))
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The shape of xx (" << xx.rows() << "x" << xx.cols() << ") does not match the shape of this ("
                                     << rows()
                                     << "x"
                                     << cols()
                                     << ")!");
    backend_ref.axpy(alpha, xx.backend());
  } // ... axpy(...)

  bool has_equal_shape(const ThisType& other) const
  {
    return (rows() == other.rows()) && (cols() == other.cols());
  }

  /// \}
  /// \name Required by MatrixInterface.
  /// \{

  inline size_t rows() const
  {
    return backend_->N();
  }

  inline size_t cols() const
  {
    return backend_->M();
  }

  inline void mv(const IstlDenseVector<ScalarType>& xx, IstlDenseVector<ScalarType>& yy) const
  {
    DUNE_XT_COMMON_TIMING_SCOPE(static_id() + ".mv");
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    backend_ref.mv(xx.backend(), yy.backend());
  }

  void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    assert(these_are_valid_indices(ii, jj));
    backend_ref[ii][jj][0][0] += value;
  }

  void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    assert(these_are_valid_indices(ii, jj));
    backend_ref[ii][jj][0][0] = value;
  }

  ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    assert(ii < rows());
    assert(jj < cols());
    if (these_are_valid_indices(ii, jj))
      return backend_->operator[](ii)[jj][0][0];
    else
      return ScalarType(0);
  } // ... get_entry(...)

  ScalarType& get_entry_ref(const size_t ii, const size_t jj)
  {
    assert(ii < rows());
    assert(jj < cols());
    if (!these_are_valid_indices(ii, jj))
      DUNE_THROW(Dune::RangeError, "Matrix entry not in pattern!");
    return backend_->operator[](ii)[jj][0][0];
  } // ... get_entry(...)

  void clear_row(const size_t ii)
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    if (ii >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    backend_ref[ii] *= ScalarType(0);
  } // ... clear_row(...)

  void clear_col(const size_t jj)
  {
    ensure_uniqueness();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    if (jj >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    for (size_t ii = 0; ii < rows(); ++ii) {
      auto& row = backend_->operator[](ii);
      const auto search_result = row.find(jj);
      if (search_result != row.end())
        row.operator[](jj)[0][0] = ScalarType(0);
    }
  } // ... clear_col(...)

  void unit_row(const size_t ii)
  {
    ensure_uniqueness();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    if (ii >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the cols of this (" << cols() << ")!");
    if (ii >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    if (!backend_->exists(ii, ii))
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Diagonal entry (" << ii << ", " << ii << ") is not contained in the sparsity pattern!");
    backend_->operator[](ii) *= ScalarType(0);
    backend_->operator[](ii)[ii] = ScalarType(1);
  } // ... unit_row(...)

  void unit_col(const size_t jj)
  {
    ensure_uniqueness();
    if (jj >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    if (jj >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the rows of this (" << rows() << ")!");
    if (!backend_->exists(jj, jj))
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Diagonal entry (" << jj << ", " << jj << ") is not contained in the sparsity pattern!");
    mutex_.lock();
    for (size_t ii = 0; (ii < rows()) && (ii != jj); ++ii) {
      auto& row = backend_->operator[](ii);
      const auto search_result = row.find(jj);
      if (search_result != row.end())
        row.operator[](jj)[0][0] = ScalarType(0);
    }
    mutex_.unlock();
    set_entry(jj, jj, ScalarType(1));
  } // ... unit_col(...)

  bool valid() const
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    for (size_t ii = 0; ii < rows(); ++ii) {
      const auto& row_vec = backend_->operator[](ii);
      for (size_t jj = 0; jj < cols(); ++jj)
        if (backend_->exists(ii, jj)) {
          const auto& entry = row_vec[jj][0];
          if (Common::isnan(entry[0]) || Common::isinf(entry[0]))
            return false;
        }
    }
    return true;
  } // ... valid(...)

  /**
   * \attention Use and interprete with care, since the Dune::BCRSMatrix is known to report strange things here,
   * depending on its state!
   */
  virtual size_t non_zeros() const override final
  {
    return backend_->nonzeroes();
  }

  virtual SparsityPatternDefault pattern(const bool prune = false,
                                         const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
                                             Common::FloatCmp::DefaultEpsilon<ScalarType>::value()) const override final
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    SparsityPatternDefault ret(rows());
    if (prune) {
      return pruned_pattern_from_backend(*backend_, eps);
    } else {
      for (size_t ii = 0; ii < rows(); ++ii) {
        if (backend_->getrowsize(ii) > 0) {
          const auto& row = backend_->operator[](ii);
          const auto it_end = row.end();
          for (auto it = row.begin(); it != it_end; ++it)
            ret.insert(ii, it.index());
        }
      }
    }
    ret.sort();
    return ret;
  } // ... pattern(...)

  virtual ThisType pruned(const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
                              Common::FloatCmp::DefaultEpsilon<ScalarType>::value()) const override final
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    return ThisType(*backend_, true, eps);
  }

  /// \}

private:
  void build_sparse_matrix(const size_t rr, const size_t cc, const SparsityPatternDefault& patt)
  {
    DUNE_XT_COMMON_TIMING_SCOPE(static_id() + ".build");
    backend_ = std::make_shared<BackendType>(rr, cc, BackendType::random);
    for (size_t ii = 0; ii < patt.size(); ++ii)
      backend_->setrowsize(ii, patt.inner(ii).size());
    backend_->endrowsizes();
    for (size_t ii = 0; ii < patt.size(); ++ii)
      for (const auto& jj : patt.inner(ii))
        backend_->addindex(ii, jj);
    backend_->endindices();
  } // ... build_sparse_matrix(...)

  SparsityPatternDefault
  pruned_pattern_from_backend(const BackendType& mat,
                              const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
                                  Common::FloatCmp::DefaultEpsilon<ScalarType>::value()) const
  {
    SparsityPatternDefault ret(mat.N());
    for (size_t ii = 0; ii < mat.N(); ++ii) {
      if (mat.getrowsize(ii) > 0) {
        const auto& row = mat[ii];
        const auto it_end = row.end();
        for (auto it = row.begin(); it != it_end; ++it) {
          const auto val = it->operator[](0)[0];
          if (Common::FloatCmp::ne<Common::FloatCmp::Style::absolute>(val, decltype(val)(0), eps))
            ret.insert(ii, it.index());
        }
      }
    }
    ret.sort();
    return ret;
  } // ... pruned_pattern_from_backend(...)

  bool these_are_valid_indices(const size_t ii, const size_t jj) const
  {
    if (ii >= rows())
      return false;
    if (jj >= cols())
      return false;
    return backend_->exists(ii, jj);
  } // ... these_are_valid_indices(...)

protected:
  /**
   * \see ContainerInterface
   */
  inline void ensure_uniqueness()
  {
    if (!backend_.unique())
      backend_ = std::make_shared<BackendType>(*backend_);
  } // ... ensure_uniqueness(...)

private:
  std::shared_ptr<BackendType> backend_;
  mutable std::mutex mutex_;
}; // class IstlRowMajorSparseMatrix

template <class S>
std::ostream& operator<<(std::ostream& out, const IstlRowMajorSparseMatrix<S>& matrix)
{
  out << "[";
  const size_t rows = matrix.rows();
  const size_t cols = matrix.cols();
  if (rows > 0 && cols > 0) {
    for (size_t ii = 0; ii < rows; ++ii) {
      if (ii > 0)
        out << "\n ";
      out << "[";
      if (matrix.backend().exists(ii, 0))
        out << matrix.get_entry(ii, 0);
      else
        out << "0";
      for (size_t jj = 1; jj < cols; ++jj) {
        out << " ";
        if (matrix.backend().exists(ii, jj))
          out << matrix.get_entry(ii, jj);
        else
          out << "0";
      }
      out << "]";
      if (rows > 1 && ii < (rows - 1))
        out << ",";
    }
    out << "]";
  } else
    out << "[ ]]";
  return out;
} // ... operator<<(...)

#else // HAVE_DUNE_ISTL

template <class ScalarImp>
class IstlDenseVector
{
  static_assert(Dune::AlwaysFalse<ScalarImp>::value, "You are missing dune-istl!");
};

template <class ScalarImp>
class IstlRowMajorSparseMatrix
{
  static_assert(Dune::AlwaysFalse<ScalarImp>::value, "You are missing dune-istl!");
};

#endif // HAVE_DUNE_ISTL

} // namespace LA
namespace Common {

#if HAVE_DUNE_ISTL

template <class T>
struct VectorAbstraction<LA::IstlDenseVector<T>> : public LA::internal::VectorAbstractionBase<LA::IstlDenseVector<T>>
{
};

template <class T>
struct MatrixAbstraction<LA::IstlRowMajorSparseMatrix<T>>
    : public LA::internal::MatrixAbstractionBase<LA::IstlRowMajorSparseMatrix<T>>
{
};

#endif // HAVE_DUNE_ISTL

} // namespace Common
} // namespace XT
} // namespace Dune


// begin: this is what we need for the lib
#if DUNE_XT_WITH_PYTHON_BINDINGS && HAVE_DUNE_ISTL


extern template class Dune::XT::LA::IstlDenseVector<double>;
extern template class Dune::XT::LA::IstlRowMajorSparseMatrix<double>;
// extern template std::ostream& operator<<(std::ostream&, const Dune::XT::LA::IstlRowMajorSparseMatrix<double>&);


#endif // DUNE_XT_WITH_PYTHON_BINDINGS && HAVE_DUNE_ISTL
// end: this is what we need for the lib


#endif // DUNE_XT_LA_CONTAINER_ISTL_HH
