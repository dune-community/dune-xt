// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright Holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_STUFF_LA_CONTAINER_COMMON_HH
#define DUNE_STUFF_LA_CONTAINER_COMMON_HH

#include <memory>
#include <type_traits>

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/stuff/common/reenable_warnings.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/typetraits.hh>

#include "interfaces.hh"
#include "pattern.hh"

namespace Dune {
namespace Stuff {
namespace LA {


// forwards
template <class ScalarImp>
class CommonDenseVector;

template <class ScalarImp>
class CommonDenseMatrix;


/// Traits for CommonDenseVector
template <class ScalarImp = double>
class CommonDenseVectorTraits
{
public:
  typedef ScalarImp ScalarType;
  typedef CommonDenseVector<ScalarType> derived_type;
  typedef Dune::DynamicVector<ScalarType> BackendType;
};


/**
 *  \brief A dense vector implementation of VectorInterface using the Dune::DynamicVector.
 */
template <class ScalarImp = double>
class CommonDenseVector : public VectorInterface<CommonDenseVectorTraits<ScalarImp>>,
                          public ProvidesBackend<CommonDenseVectorTraits<ScalarImp>>
{
  typedef CommonDenseVector<ScalarImp> ThisType;
  typedef VectorInterface<CommonDenseVectorTraits<ScalarImp>> VectorInterfaceType;
  static_assert(!std::is_same<DUNE_STUFF_SSIZE_T, int>::value,
                "You have to manually disable the constructor below which uses DUNE_STUFF_SSIZE_T!");

public:
  typedef CommonDenseVectorTraits<ScalarImp> Traits;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::BackendType BackendType;

  CommonDenseVector(const size_t ss = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(ss, value))
  {
  }

  /// This constructor is needed for the python bindings.
  CommonDenseVector(const DUNE_STUFF_SSIZE_T ss, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(VectorInterfaceType::assert_is_size_t_compatible_and_convert(ss), value))
  {
  }

  CommonDenseVector(const int ss, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(VectorInterfaceType::assert_is_size_t_compatible_and_convert(ss), value))
  {
  }

  CommonDenseVector(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  CommonDenseVector(const BackendType& other)
    : backend_(new BackendType(other))
  {
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  CommonDenseVector(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  CommonDenseVector(std::shared_ptr<BackendType> backend_ptr)
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
    ensure_uniqueness();
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
    backend() *= alpha;
  } // ... scal(...)

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    if (xx.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of x (" << xx.size() << ") does not match the size of this (" << size() << ")!");
    ensure_uniqueness();
    auto& this_ref     = *backend_;
    const auto& xx_ref = *(xx.backend_);
    for (size_t ii = 0; ii < this_ref.size(); ++ii)
      this_ref[ii] += alpha * xx_ref[ii];
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
    return backend_->size();
  }

  void add_to_entry(const size_t ii, const ScalarType& value)
  {
    assert(ii < size());
    ensure_uniqueness();
    backend_->operator[](ii) += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const ScalarType& value)
  {
    assert(ii < size());
    ensure_uniqueness();
    backend_->operator[](ii) = value;
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii) const
  {
    assert(ii < size());
    return backend_->operator[](ii);
  } // ... get_entry(...)

private:
  inline ScalarType& get_entry_ref(const size_t ii)
  {
    return backend()[ii];
  }

  inline const ScalarType& get_entry_ref(const size_t ii) const
  {
    return backend_->operator[](ii);
  }

public:
  /**
   * \}
   */

  /**
   * \defgroup vector_overrides ´´These methods override default implementations from VectorInterface.``
   * \{
   */

  virtual ScalarType dot(const ThisType& other) const DS_OVERRIDE DS_FINAL
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    return backend_->operator*(*(other.backend_));
  } // ... dot(...)

  virtual ScalarType l1_norm() const DS_OVERRIDE DS_FINAL
  {
    return backend_->one_norm();
  }

  virtual ScalarType l2_norm() const DS_OVERRIDE DS_FINAL
  {
    return backend_->two_norm();
  }

  virtual ScalarType sup_norm() const DS_OVERRIDE DS_FINAL
  {
    return backend_->infinity_norm();
  }

  virtual void add(const ThisType& other, ThisType& result) const DS_OVERRIDE DS_FINAL
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    if (result.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of result (" << result.size() << ") does not match the size of this (" << size()
                                                   << ")!");
    BackendType& result_ref = result.backend();
    for (size_t ii = 0; ii < size(); ++ii)
      result_ref[ii] = backend_->operator[](ii) + other.backend_->operator[](ii);
  } // ... add(...)

  virtual void iadd(const ThisType& other) DS_OVERRIDE DS_FINAL
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    ensure_uniqueness();
    backend_->operator+=(*(other.backend_));
  } // ... iadd(...)

  virtual void sub(const ThisType& other, ThisType& result) const DS_OVERRIDE DS_FINAL
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    if (result.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of result (" << result.size() << ") does not match the size of this (" << size()
                                                   << ")!");
    BackendType& result_ref = result.backend();
    for (size_t ii = 0; ii < size(); ++ii)
      result_ref[ii] = backend_->operator[](ii)-other.backend_->operator[](ii);
  } // ... sub(...)

  virtual void isub(const ThisType& other) DS_OVERRIDE DS_FINAL
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

  /**
   * \defgroup vector_defaults ´´These methods are imported from VectorInterface.``
   * \{
   */
  using VectorInterfaceType::add;
  using VectorInterfaceType::sub;
  /**
   * \}
   */

private:
  inline void ensure_uniqueness() const
  {
    if (!backend_.unique())
      backend_ = std::make_shared<BackendType>(*backend_);
  } // ... ensure_uniqueness(...)

  friend class VectorInterface<CommonDenseVectorTraits<ScalarType>>;
  friend class CommonDenseMatrix<ScalarType>;

  mutable std::shared_ptr<BackendType> backend_;
}; // class CommonDenseVector


template <class ScalarImp = double>
class CommonDenseMatrixTraits
{
public:
  typedef ScalarImp ScalarType;
  typedef CommonDenseMatrix<ScalarType> derived_type;
  typedef Dune::DynamicMatrix<ScalarType> BackendType;
};


/**
 *  \brief  A dense matrix implementation of MatrixInterface using the dune-common.
 */
template <class ScalarImp = double>
class CommonDenseMatrix : public MatrixInterface<CommonDenseMatrixTraits<ScalarImp>>,
                          public ProvidesBackend<CommonDenseMatrixTraits<ScalarImp>>
{
  typedef CommonDenseMatrix<ScalarImp> ThisType;
  typedef MatrixInterface<CommonDenseMatrixTraits<ScalarImp>> MatrixInterfaceType;
  static_assert(!std::is_same<DUNE_STUFF_SSIZE_T, int>::value,
                "You have to manually disable the constructor below which uses DUNE_STUFF_SSIZE_T!");

public:
  typedef CommonDenseMatrixTraits<ScalarImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ScalarType ScalarType;

  CommonDenseMatrix(const size_t rr = 0, const size_t cc = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(rr, cc, value))
  {
  }

  /// This constructor is needed for the python bindings.
  CommonDenseMatrix(const DUNE_STUFF_SSIZE_T rr, const DUNE_STUFF_SSIZE_T cc = 0,
                    const ScalarType value = ScalarType(0))
    : backend_(new BackendType(MatrixInterfaceType::assert_is_size_t_compatible_and_convert(rr),
                               MatrixInterfaceType::assert_is_size_t_compatible_and_convert(cc), value))
  {
  }

  CommonDenseMatrix(const int rr, const int cc = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(MatrixInterfaceType::assert_is_size_t_compatible_and_convert(rr),
                               MatrixInterfaceType::assert_is_size_t_compatible_and_convert(cc), value))
  {
  }

  /// This constructors ignores the given pattern and initializes the matrix with 0.
  CommonDenseMatrix(const size_t rr, const size_t cc, const SparsityPatternDefault& /*pattern*/)
    : backend_(new BackendType(MatrixInterfaceType::assert_is_size_t_compatible_and_convert(rr),
                               MatrixInterfaceType::assert_is_size_t_compatible_and_convert(cc), ScalarType(0)))
  {
  }

  CommonDenseMatrix(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  CommonDenseMatrix(const BackendType& other)
    : backend_(new BackendType(other))
  {
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  CommonDenseMatrix(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  CommonDenseMatrix(std::shared_ptr<BackendType> backend_ptr)
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
    ensure_uniqueness();
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
    return backend_->rows();
  }

  inline size_t cols() const
  {
    return backend_->cols();
  }

  inline void mv(const VectorInterface<CommonDenseVectorTraits<ScalarType>>& xx,
                 VectorInterface<CommonDenseVectorTraits<ScalarType>>& yy) const
  {
    mv(static_cast<const typename CommonDenseVectorTraits<ScalarType>::derived_type&>(xx),
       static_cast<typename CommonDenseVectorTraits<ScalarType>::derived_type&>(yy));
  }

  inline void mv(const CommonDenseVector<ScalarType>& xx, CommonDenseVector<ScalarType>& yy) const
  {
    backend_->mv(*(xx.backend_), yy.backend());
  }

  void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(ii < rows());
    assert(jj < cols());
    ensure_uniqueness();
    backend_->operator[](ii)[jj] += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(ii < rows());
    assert(jj < cols());
    ensure_uniqueness();
    backend_->operator[](ii)[jj] = value;
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    assert(ii < rows());
    assert(jj < cols());
    return backend_->operator[](ii)[jj];
  } // ... get_entry(...)

  void clear_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    backend()[ii] *= ScalarType(0);
  } // ... clear_row(...)

  void clear_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    BackendType& backend_ref = backend();
    for (size_t ii = 0; ii < rows(); ++ii)
      backend_ref[ii][jj] = ScalarType(0);
  } // ... clear_col(...)

  void unit_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    auto& row = backend()[ii];
    for (size_t jj = 0; jj < cols(); ++jj)
      row[jj] = ScalarType(0);
    row[ii]   = ScalarType(1);
  } // ... unit_row(...)

  void unit_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    ensure_uniqueness();
    for (size_t ii = 0; ii < rows(); ++ii)
      backend_->operator[](ii)[jj] = ScalarType(0);
    backend_->operator[](jj)[jj] = ScalarType(1);
  } // ... unit_col(...)

  /**
   * \}
   */

private:
  inline void ensure_uniqueness() const
  {
    if (!backend_.unique())
      backend_ = std::make_shared<BackendType>(*backend_);
  } // ... ensure_uniqueness(...)

  mutable std::shared_ptr<BackendType> backend_;
}; // class CommonDenseMatrix


} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_COMMON_HH
