// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright Holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_CONTAINER_ISTL_HH
#define DUNE_STUFF_LA_CONTAINER_ISTL_HH

#include <dune/common/static_assert.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#endif // HAVE_DUNE_ISTL

#include "interfaces.hh"
#include "pattern.hh"

namespace Dune {
namespace Stuff {
namespace LA {


// forward
template <class ScalarImp>
class IstlDenseVector;

template <class ScalarImp>
class IstlRowMajorSparseMatrix;


#if HAVE_DUNE_ISTL

/// Traits for IstlDenseVector.
template <class ScalarImp>
class IstlDenseVectorTraits
{
public:
  typedef ScalarImp ScalarType;
  typedef IstlDenseVector<ScalarImp> derived_type;
  typedef BlockVector<FieldVector<ScalarType, 1>> BackendType;
}; // class IstlDenseVectorTraits


/**
 *  \brief A dense vector implementation of VectorInterface using the Dune::BlockVector from dune-istl.
 */
template <class ScalarImp = double>
class IstlDenseVector : public VectorInterface<IstlDenseVectorTraits<ScalarImp>>,
                        public ProvidesBackend<IstlDenseVectorTraits<ScalarImp>>
{
  typedef IstlDenseVector<ScalarImp> ThisType;
  typedef VectorInterface<IstlDenseVectorTraits<ScalarImp>> VectorInterfaceType;
  static_assert(!std::is_same<DUNE_STUFF_SSIZE_T, int>::value,
                "You have to manually disable the constructor below which uses DUNE_STUFF_SSIZE_T!");

public:
  typedef IstlDenseVectorTraits<ScalarImp> Traits;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::BackendType BackendType;

  IstlDenseVector(const size_t ss = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(ss))
  {
    backend_->operator=(value);
  }

  /// This constructor is needed for the python bindings.
  IstlDenseVector(const DUNE_STUFF_SSIZE_T ss, const ScalarType value = ScalarType(0))
    : IstlDenseVector(VectorInterfaceType::assert_is_size_t_compatible_and_convert(ss), value)
  {
  }

  /// This constructor is needed because marking the above one as explicit had no effect.
  IstlDenseVector(const int ss, const ScalarType value = ScalarType(0))
    : IstlDenseVector(VectorInterfaceType::assert_is_size_t_compatible_and_convert(ss), value)
  {
  }

  IstlDenseVector(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  IstlDenseVector(const BackendType& other)
    : backend_(new BackendType(other))
  {
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  IstlDenseVector(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  IstlDenseVector(std::shared_ptr<BackendType> backend_ptr)
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
    backend_ = std::make_shared<BackendType>(other);
    return *this;
  } // ... operator=(...)

  /**
   * \defgroup backend ´´These methods are required by the ProvidesBackend interface.``
   * \{
   */

  BackendType& backend()
  {
    ensure_uniqueness();
    return *backend_;
  } // ... backend(...)

  const BackendType& backend() const
  {
    const_cast<ThisType&>(*this).ensure_uniqueness();
    return *backend_;
  } // ... backend(...)
  /**
   * \}
   */

  /**
   * \defgroup container ´´These methods are required by ContainerInterface.``
   * \{
   */

  ThisType copy() const
  {
    return ThisType(*backend_);
  }

  void scal(const ScalarType& alpha)
  {
    ensure_uniqueness();
    backend_->operator*=(alpha);
  } // ... scal(...)

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    if (xx.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of x (" << xx.size() << ") does not match the size of this (" << size() << ")!");
    ensure_uniqueness();
    backend_->axpy(alpha, *(xx.backend_));
  } // ... axpy(...)

  bool has_equal_shape(const ThisType& other) const
  {
    return size() == other.size();
  }
  /**
   * \}
   */

  /**
   * \defgroup vector_required ´´These methods are required by VectorInterface.``
   * \{
   */

  inline size_t size() const
  {
    return backend_->dim();
  }

  void add_to_entry(const size_t ii, const ScalarType& value)
  {
    assert(ii < size());
    ensure_uniqueness();
    backend_->operator[](ii)[0] += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const ScalarType& value)
  {
    assert(ii < size());
    ensure_uniqueness();
    backend_->operator[](ii)[0] = value;
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii) const
  {
    assert(ii < size());
    return backend_->operator[](ii)[0];
  } // ... get_entry(...)

private:
  inline ScalarType& get_entry_ref(const size_t ii)
  {
    return backend_->operator[](ii)[0];
  }

  inline const ScalarType& get_entry_ref(const size_t ii) const
  {
    return backend_->operator[](ii)[0];
  }

public:
  /**
   * \}
   */

  /**
   * \defgroup vector_overrides ´´These methods override default implementations from VectorInterface.``
   * \{
   */

  virtual ScalarType dot(const ThisType& other) const DS_OVERRIDE
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    return backend_->dot(*(other.backend_));
  } // ... dot(...)

  virtual ScalarType l1_norm() const DS_OVERRIDE
  {
    return backend_->one_norm();
  }

  virtual ScalarType l2_norm() const DS_OVERRIDE
  {
    return backend_->two_norm();
  }

  virtual ScalarType sup_norm() const DS_OVERRIDE
  {
    return backend_->infinity_norm();
  }

  virtual void add(const ThisType& other, ThisType& result) const DS_OVERRIDE
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    if (result.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of result (" << result.size() << ") does not match the size of this (" << size()
                                                   << ")!");
    result.backend() = *(backend_);
    result.backend() += *(other.backend_);
  } // ... add(...)

  virtual ThisType add(const ThisType& other) const DS_OVERRIDE
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    ThisType result = copy();
    result.backend_->operator+=(*(other.backend_));
    return result;
  } // ... add(...)

  virtual void iadd(const ThisType& other) DS_OVERRIDE
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    ensure_uniqueness();
    backend_->operator+=(*(other.backend_));
  } // ... iadd(...)

  virtual void sub(const ThisType& other, ThisType& result) const DS_OVERRIDE
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    if (result.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of result (" << result.size() << ") does not match the size of this (" << size()
                                                   << ")!");
    result.backend() = *(backend_);
    result.backend() -= *(other.backend_);
  } // ... sub(...)

  virtual ThisType sub(const ThisType& other) const DS_OVERRIDE
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    ThisType result = copy();
    result.backend_->operator-=(*(other.backend_));
    return result;
  } // ... sub(...)

  virtual void isub(const ThisType& other) DS_OVERRIDE
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    ensure_uniqueness();
    backend_->operator-=(*(other.backend_));
  } // ... isub(...)

  /**
   * \}
   */

private:
  inline void ensure_uniqueness()
  {
    if (!backend_.unique())
      backend_ = std::make_shared<BackendType>(*backend_);
  } // ... ensure_uniqueness(...)

  friend class VectorInterface<IstlDenseVectorTraits<ScalarType>>;
  friend class IstlRowMajorSparseMatrix<ScalarType>;

  std::shared_ptr<BackendType> backend_;
}; // class IstlDenseVector


/**
 * \brief Traits for IstlRowMajorSparseMatrix.
 */
template <class ScalarImp = double>
class IstlRowMajorSparseMatrixTraits
{
public:
  typedef ScalarImp ScalarType;
  typedef IstlRowMajorSparseMatrix<ScalarType> derived_type;
  typedef BCRSMatrix<FieldMatrix<ScalarType, 1, 1>> BackendType;
}; // class RowMajorSparseMatrixTraits

/**
 * \brief A sparse matrix implementation of the MatrixInterface using the Dune::BCRSMatrix from dune-istl.
 */
template <class ScalarImp = double>
class IstlRowMajorSparseMatrix : public MatrixInterface<IstlRowMajorSparseMatrixTraits<ScalarImp>>,
                                 public ProvidesBackend<IstlRowMajorSparseMatrixTraits<ScalarImp>>
{
  typedef IstlRowMajorSparseMatrix<ScalarImp> ThisType;
  //  typedef MatrixInterface< IstlRowMajorSparseMatrixTraits< ScalarImp > > MatrixInterface;
  static_assert(!std::is_same<DUNE_STUFF_SSIZE_T, int>::value,
                "You have to manually disable the constructor below which uses DUNE_STUFF_SSIZE_T!");

public:
  typedef IstlRowMajorSparseMatrixTraits<ScalarImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ScalarType ScalarType;

  /**
   * \brief This is the constructor of interest which creates a sparse matrix.
   */
  IstlRowMajorSparseMatrix(const size_t rr, const size_t cc, const SparsityPatternDefault& pattern)
    : backend_(new BackendType(rr, cc, BackendType::row_wise))
  {
    if (size_t(pattern.size()) != rr)
      DUNE_THROW_COLORFULLY(
          Exceptions::shapes_do_not_match,
          "The size of the pattern (" << pattern.size() << ") does not match the number of rows of this (" << rows()
                                      << ")!");
    size_t row_index = 0;
    for (auto row = backend_->createbegin(); row != backend_->createend(); ++row) {
      assert(row_index < pattern.size());
      const auto& cols = pattern.inner(row_index);
      for (const auto& col : cols)
        row.insert(col);
      ++row_index;
    }
    backend_->operator*=(ScalarType(0));
  }

  IstlRowMajorSparseMatrix(const size_t rr = 0, const size_t cc = 0)
    : backend_(new BackendType(rr, cc, BackendType::row_wise))
  {
  }

  /// This constructor is needed for the python bindings.
  IstlRowMajorSparseMatrix(const DUNE_STUFF_SSIZE_T rr, const DUNE_STUFF_SSIZE_T cc = 0)
    : backend_(new BackendType(this->assert_is_size_t_compatible_and_convert(rr),
                               this->assert_is_size_t_compatible_and_convert(cc), BackendType::row_wise))
  {
  }

  /// This constructor is needed because marking the above one as explicit had no effect.
  IstlRowMajorSparseMatrix(const int rr, const int cc = 0)
    : backend_(new BackendType(this->assert_is_size_t_compatible_and_convert(rr),
                               this->assert_is_size_t_compatible_and_convert(cc), BackendType::row_wise))
  {
  }

  IstlRowMajorSparseMatrix(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  IstlRowMajorSparseMatrix(const BackendType& other)
    : backend_(new BackendType(other))
  {
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  IstlRowMajorSparseMatrix(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  IstlRowMajorSparseMatrix(std::shared_ptr<BackendType> backend_ptr)
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
    backend_ = std::make_shared<BackendType>(other);
    return *this;
  } // ... operator=(...)

  /**
   * \defgroup backend ´´These methods are required by the ProvidesBackend interface.``
   * \{
   */

  BackendType& backend()
  {
    ensure_uniqueness();
    return *backend_;
  } // ... backend(...)

  const BackendType& backend() const
  {
    const_cast<ThisType&>(*this).ensure_uniqueness();
    return *backend_;
  } // ... backend(...)
  /**
   * \}
   */

  /**
   * \defgroup container ´´These methods are required by ContainerInterface.``
   * \{
   */

  ThisType copy() const
  {
    return ThisType(*backend_);
  }

  void scal(const ScalarType& alpha)
  {
    ensure_uniqueness();
    backend_->operator*=(alpha);
  } // ... scal(...)

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    if (!has_equal_shape(xx))
      DUNE_THROW_COLORFULLY(
          Exceptions::shapes_do_not_match,
          "The shape of xx (" << xx.rows() << "x" << xx.cols() << ") does not match the shape of this (" << rows()
                              << "x"
                              << cols()
                              << ")!");
    ensure_uniqueness();
    backend_->axpy(alpha, *(xx.backend_));
  } // ... axpy(...)

  bool has_equal_shape(const ThisType& other) const
  {
    return (rows() == other.rows()) && (cols() == other.cols());
  }
  /**
   * \}
   */

  /**
   * \defgroup matrix_required ´´These methods are required by MatrixInterface.``
   * \{
   */

  inline size_t rows() const
  {
    return backend_->N();
  }

  inline size_t cols() const
  {
    return backend_->M();
  }

  template <class SourceType, class RangeType>
  inline void mv(const SourceType& /*xx*/, RangeType& /*yy*/) const
  {
    static_assert(Dune::AlwaysFalse<SourceType>::value, "Not available for this combination of xx and yy!");
  }

  inline void mv(const IstlDenseVector<ScalarType>& xx, IstlDenseVector<ScalarType>& yy) const
  {
    backend_->mv(xx.backend(), yy.backend());
  }

  void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(these_are_valid_indices(ii, jj));
    backend_->operator[](ii)[jj][0][0] += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(these_are_valid_indices(ii, jj));
    backend_->operator[](ii)[jj][0][0] = value;
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    assert(ii < rows());
    assert(jj < cols());
    return backend_->operator[](ii)[jj][0][0];
  } // ... get_entry(...)

  void clear_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    backend_->operator[](ii) *= ScalarType(0);
  } // ... clear_row(...)

  void clear_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    for (size_t ii = 0; ii < rows(); ++ii) {
      auto& row                = backend_->operator[](ii);
      const auto search_result = row.find(jj);
      if (search_result != row.end())
        row.operator[](jj)[0][0] = ScalarType(0);
    }
  } // ... clear_col(...)

  void unit_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    if (!backend_->exists(ii, ii))
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Diagonal entry (" << ii << ", " << ii << ") is not contained in the sparsity pattern!");
    backend_->operator[](ii) *= ScalarType(0);
    backend_->operator[](ii)[ii] = ScalarType(1);
  } // ... unit_row(...)

  void unit_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    if (!backend_->exists(jj, jj))
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Diagonal entry (" << jj << ", " << jj << ") is not contained in the sparsity pattern!");
    for (size_t ii = 0; (ii < rows()) && (ii != jj); ++ii) {
      auto& row                = backend_->operator[](ii);
      const auto search_result = row.find(jj);
      if (search_result != row.end())
        row.operator[](jj)[0][0] = ScalarType(0);
    }
    set_entry(jj, jj, ScalarType(1));
  } // ... unit_col(...)
  /**
   * \}
   */

private:
  bool these_are_valid_indices(const size_t ii, const size_t jj)
  {
    if (ii >= rows())
      return false;
    if (jj >= cols())
      return false;
    return backend_->exists(ii, jj);
  } // ... these_are_valid_indices(...)

  inline void ensure_uniqueness()
  {
    if (!backend_.unique())
      backend_ = std::make_shared<BackendType>(*backend_);
  } // ... ensure_uniqueness(...)

  std::shared_ptr<BackendType> backend_;
}; // class IstlRowMajorSparseMatrix


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
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_ISTL_HH
