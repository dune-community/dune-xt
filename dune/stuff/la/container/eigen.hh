// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_CONTAINER_EIGEN_HH
#define DUNE_STUFF_LA_CONTAINER_EIGEN_HH

#include <memory>
#include <type_traits>
#include <vector>
#include <initializer_list>

#if HAVE_EIGEN
#include <dune/stuff/common/disable_warnings.hh>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <dune/stuff/common/reenable_warnings.hh>
#endif // HAVE_EIGEN

#include <dune/common/typetraits.hh>
#include <dune/stuff/common/disable_warnings.hh>
#include <dune/common/densematrix.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/common/float_cmp.hh>

#include "interfaces.hh"
#include "pattern.hh"

namespace Dune {
namespace Stuff {
namespace LA {


// forwards
template <class Traits, class ScalarImp>
class EigenBaseVector;

template <class ScalarImp>
class EigenDenseVector;

template <class T>
class EigenMappedDenseVector;

template <class ScalarImp>
class EigenDenseMatrix;

template <class ScalarType>
class EigenRowMajorSparseMatrix;


class EigenVectorInterfaceDynamic
{
};
class EigenMatrixInterfaceDynamic
{
};


#if HAVE_EIGEN


/**
 *  \brief Base class for all eigen implementations of VectorInterface.
 */
template <class ImpTraits, class ScalarImp = double>
class EigenBaseVector : public VectorInterface<ImpTraits>,
                        public EigenVectorInterfaceDynamic,
                        public ProvidesBackend<ImpTraits>
{
  typedef VectorInterface<ImpTraits> VectorInterfaceType;

public:
  typedef ImpTraits Traits;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::derived_type VectorImpType;

private:
//! disambiguation necessary since it exeists in multiple bases
#ifndef NDEBUG
  using VectorInterfaceType::crtp_mutex_;
#endif

  /**
   * \see ContainerInterface
   */
  void ensure_uniqueness() const
  {
    CHECK_AND_CALL_CRTP(VectorInterfaceType::as_imp().ensure_uniqueness());
    VectorInterfaceType::as_imp().ensure_uniqueness();
  }

public:
  VectorImpType& operator=(const VectorImpType& other)
  {
    backend_ = other.backend_;
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
  VectorImpType copy() const
  {
    return VectorImpType(*backend_);
  }

  void scal(const ScalarType& alpha)
  {
    backend() *= alpha;
  }

  template <class T>
  void axpy(const ScalarType& alpha, const EigenBaseVector<T, ScalarType>& xx)
  {
    if (xx.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of xx (" << xx.size() << ") does not match the size of this (" << size() << ")!");
    backend() += alpha * xx.backend();
  } // ... axpy(...)

  bool has_equal_shape(const VectorImpType& other) const
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
    backend()(ii) += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const ScalarType& value)
  {
    assert(ii < size());
    backend()(ii) = value;
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii) const
  {
    assert(ii < size());
    return backend_->operator[](ii);
  } // ... get_entry(...)

protected:
  inline ScalarType& get_entry_ref(const size_t ii)
  {
    return backend()[ii];
  }

  inline const ScalarType& get_entry_ref(const size_t ii) const
  {
    return backend()[ii];
  }
  /**
   * \}
   */

public:
  /**
   * \defgroup vector_overrides ´´These methods override default implementations from VectorInterface.``
   * \{
   */
  virtual std::pair<size_t, ScalarType> amax() const DS_OVERRIDE DS_FINAL
  {
    auto result              = std::make_pair(size_t(0), ScalarType(0));
    size_t min_index         = 0;
    size_t max_index         = 0;
    const ScalarType minimum = backend_->minCoeff(&min_index);
    const ScalarType maximum = backend_->maxCoeff(&max_index);
    if (std::abs(maximum) < std::abs(minimum)
        || (Common::FloatCmp::eq(std::abs(maximum), std::abs(minimum)) && max_index > min_index)) {
      result.first  = min_index;
      result.second = std::abs(minimum);
    } else {
      result.first  = max_index;
      result.second = std::abs(maximum);
    }
    return result;
  } // ... amax(...)

  template <class T>
  bool almost_equal(const EigenBaseVector<T, ScalarType>& other,
                    const ScalarType epsilon = Dune::FloatCmp::DefaultEpsilon<ScalarType>::value()) const
  {
    if (other.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    for (size_t ii = 0; ii < size(); ++ii)
      if (!Dune::FloatCmp::eq<ScalarType>(get_entry(ii), other.get_entry(ii), epsilon))
        return false;
    return true;
  } // ... almost_equal(...)

  virtual bool almost_equal(const VectorImpType& other,
                            const ScalarType epsilon = Dune::FloatCmp::DefaultEpsilon<ScalarType>::value()) const
      DS_OVERRIDE DS_FINAL
  {
    return this->template almost_equal<Traits>(other, epsilon);
  }

  using VectorInterfaceType::almost_equal;

  template <class T>
  ScalarType dot(const EigenBaseVector<T, ScalarType>& other) const
  {
    if (other.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    return backend_->transpose() * *(other.backend_);
  } // ... dot(...)

  virtual ScalarType dot(const VectorImpType& other) const DS_OVERRIDE DS_FINAL
  {
    return this->template dot<Traits>(other);
  }

  virtual ScalarType l1_norm() const DS_OVERRIDE DS_FINAL
  {
    return backend_->template lpNorm<1>();
  }

  virtual ScalarType l2_norm() const DS_OVERRIDE DS_FINAL
  {
    return backend_->template lpNorm<2>();
  }

  virtual ScalarType sup_norm() const DS_OVERRIDE DS_FINAL
  {
    return backend_->template lpNorm<::Eigen::Infinity>();
  }

  template <class T1, class T2>
  void add(const EigenBaseVector<T1, ScalarType>& other, EigenBaseVector<T2, ScalarType>& result) const
  {
    if (other.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    if (result.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of result (" << result.size() << ") does not match the size of this (" << size() << ")!");
    result.backend() = *backend_ + *(other.backend_);
  } // ... add(...)

  virtual void add(const VectorImpType& other, VectorImpType& result) const DS_OVERRIDE DS_FINAL
  {
    return this->template add<Traits, Traits>(other, result);
  }

  template <class T>
  void iadd(const EigenBaseVector<T, ScalarType>& other)
  {
    if (other.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    backend() += *(other.backend_);
  } // ... iadd(...)

  virtual void iadd(const VectorImpType& other) DS_OVERRIDE DS_FINAL
  {
    return this->template iadd<Traits>(other);
  }

  template <class T1, class T2>
  void sub(const EigenBaseVector<T1, ScalarType>& other, EigenBaseVector<T2, ScalarType>& result) const
  {
    if (other.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    if (result.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of result (" << result.size() << ") does not match the size of this (" << size() << ")!");
    result.backend() = *backend_ - *(other.backend_);
  } // ... sub(...)

  virtual void sub(const VectorImpType& other, VectorImpType& result) const DS_OVERRIDE DS_FINAL
  {
    return this->template sub<Traits, Traits>(other, result);
  }

  template <class T>
  void isub(const EigenBaseVector<T, ScalarType>& other)
  {
    if (other.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    backend() -= *(other.backend_);
  } // ... isub(...)

  virtual void isub(const VectorImpType& other) DS_OVERRIDE DS_FINAL
  {
    this->template isub<Traits>(other);
  }
  /**
   * \}
   */

private:
  friend class VectorInterface<Traits>;
  friend class EigenDenseMatrix<ScalarType>;
  friend class EigenRowMajorSparseMatrix<ScalarType>;

protected:
  mutable std::shared_ptr<BackendType> backend_;
}; // class EigenBaseVector


/**
 *  \brief Traits for EigenDenseVector.
 */
template <class ScalarImp = double>
class EigenDenseVectorTraits
{
public:
  typedef ScalarImp ScalarType;
  typedef EigenDenseVector<ScalarType> derived_type;
  typedef typename ::Eigen::Matrix<ScalarType, ::Eigen::Dynamic, 1> BackendType;
}; // class EigenDenseVectorTraits


/**
 *  \brief A dense vector implementation of VectorInterface using the eigen backend.
 */
template <class ScalarImp = double>
class EigenDenseVector : public EigenBaseVector<EigenDenseVectorTraits<ScalarImp>>,
                         public ProvidesDataAccess<EigenDenseVectorTraits<ScalarImp>>
{
  typedef EigenDenseVector<ScalarImp> ThisType;
  typedef VectorInterface<EigenDenseVectorTraits<ScalarImp>> VectorInterfaceType;
  typedef EigenBaseVector<EigenDenseVectorTraits<ScalarImp>> BaseType;
  static_assert(!std::is_same<DUNE_STUFF_SSIZE_T, int>::value,
                "You have to manually disable the constructor below which uses DUNE_STUFF_SSIZE_T!");

public:
  typedef EigenDenseVectorTraits<ScalarImp> Traits;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::BackendType BackendType;

  EigenDenseVector(const size_t ss = 0, const ScalarType value = ScalarType(0))
  {
    this->backend_ = std::make_shared<BackendType>(ss);
    if (FloatCmp::eq(value, ScalarType(0)))
      this->backend_->setZero();
    else {
      this->backend_->setOnes();
      this->backend_->operator*=(value);
    }
  }

  /// This constructor is needed for the python bindings.
  EigenDenseVector(const DUNE_STUFF_SSIZE_T ss, const ScalarType value = ScalarType(0))
  {
    this->backend_ = std::make_shared<BackendType>(VectorInterfaceType::assert_is_size_t_compatible_and_convert(ss));
    if (FloatCmp::eq(value, ScalarType(0)))
      this->backend_->setZero();
    else {
      this->backend_->setOnes();
      this->backend_->operator*=(value);
    }
  }

  EigenDenseVector(const int ss, const ScalarType value = ScalarType(0))
  {
    this->backend_ = std::make_shared<BackendType>(VectorInterfaceType::assert_is_size_t_compatible_and_convert(ss));
    if (FloatCmp::eq(value, ScalarType(0)))
      this->backend_->setZero();
    else {
      this->backend_->setOnes();
      this->backend_->operator*=(value);
    }
  }

  EigenDenseVector(const std::vector<ScalarType>& other)
  {
    this->backend_ = std::make_shared<BackendType>(other.size());
    for (size_t ii = 0; ii < other.size(); ++ii)
      this->backend_->operator[](ii) = other[ii];
  }

  EigenDenseVector(const std::initializer_list<ScalarType>& other)
  {
    this->backend_ = std::make_shared<BackendType>(other.size());
    size_t ii = 0;
    for (auto element : other) {
      this->backend_->operator[](ii) = element;
      ++ii;
    }
  }

  EigenDenseVector(const BackendType& other)
  {
    this->backend_ = std::make_shared<BackendType>(other);
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  EigenDenseVector(BackendType* backend_ptr)
  {
    this->backend_ = std::shared_ptr<BackendType>(backend_ptr);
  }

  EigenDenseVector(std::shared_ptr<BackendType> backend_ptr)
  {
    this->backend_ = backend_ptr;
  }

  /**
   *  \note Does a deep copy.
   */
  ThisType& operator=(const BackendType& other)
  {
    this->backend_ = std::make_shared<BackendType>(other);
    return *this;
  } // ... operator=(...)

  using VectorInterfaceType::add;
  using VectorInterfaceType::sub;
  using BaseType::backend;

private:
  inline void ensure_uniqueness() const
  {
    if (!this->backend_.unique())
      this->backend_ = std::make_shared<BackendType>(*(this->backend_));
  } // ... ensure_uniqueness(...)

  friend class EigenBaseVector<EigenDenseVectorTraits<ScalarType>>;
}; // class EigenDenseVector


/**
 *  \brief Traits for EigenMappedDenseVector.
 */
template <class ScalarImp = double>
class EigenMappedDenseVectorTraits
{
  typedef typename ::Eigen::Matrix<ScalarImp, ::Eigen::Dynamic, 1> PlainBackendType;

public:
  typedef ScalarImp ScalarType;
  typedef EigenMappedDenseVector<ScalarType> derived_type;
  typedef Eigen::Map<PlainBackendType> BackendType;
};


/**
 *  \brief  A dense vector implementation of VectorInterface using the eigen backend which wrappes a raw array.
 */
template <class ScalarImp = double>
class EigenMappedDenseVector : public EigenBaseVector<EigenMappedDenseVectorTraits<ScalarImp>>,
                               public ProvidesBackend<EigenMappedDenseVectorTraits<ScalarImp>>
{
  typedef EigenMappedDenseVector<ScalarImp> ThisType;
  typedef VectorInterface<EigenMappedDenseVectorTraits<ScalarImp>> VectorInterfaceType;
  typedef EigenBaseVector<EigenMappedDenseVectorTraits<ScalarImp>> BaseType;
  static_assert(std::is_same<ScalarImp, double>::value, "Undefined behaviour for non-double data!");
  static_assert(!std::is_same<DUNE_STUFF_SSIZE_T, int>::value,
                "You have to manually disable the constructor below which uses DUNE_STUFF_SSIZE_T!");

public:
  typedef EigenMappedDenseVectorTraits<ScalarImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ScalarType ScalarType;

  /**
   *  \brief  This is the constructor of interest which wrappes a raw array.
   */
  EigenMappedDenseVector(ScalarType* data, size_t data_size)
  {
    this->backend_ = std::make_shared<BackendType>(data, data_size);
  }

  /**
   *  \brief  This constructor allows to create an instance of this type just like any other vector.
   */
  EigenMappedDenseVector(const size_t ss = 0, const ScalarType value = ScalarType(0))
  {
    this->backend_ = std::make_shared<BackendType>(new ScalarType[ss], ss);
    if (FloatCmp::eq(value, ScalarType(0)))
      this->backend_->setZero();
    else {
      this->backend_->setOnes();
      this->backend_->operator*=(value);
    }
  }

  /// This constructor is needed for the python bindings.
  EigenMappedDenseVector(const DUNE_STUFF_SSIZE_T ss, const ScalarType value = ScalarType(0))
  {
    const size_t ss_size_t = VectorInterfaceType::assert_is_size_t_compatible_and_convert(ss);
    this->backend_ = std::make_shared<BackendType>(new ScalarType[ss_size_t], ss_size_t);
    if (FloatCmp::eq(value, ScalarType(0)))
      this->backend_->setZero();
    else {
      this->backend_->setOnes();
      this->backend_->operator*=(value);
    }
  }

  EigenMappedDenseVector(const int ss, const ScalarType value = ScalarType(0))
  {
    const size_t ss_size_t = VectorInterfaceType::assert_is_size_t_compatible_and_convert(ss);
    this->backend_ = std::make_shared<BackendType>(new ScalarType[ss_size_t], ss_size_t);
    if (FloatCmp::eq(value, ScalarType(0)))
      this->backend_->setZero();
    else {
      this->backend_->setOnes();
      this->backend_->operator*=(value);
    }
  }

  EigenMappedDenseVector(const std::vector<ScalarType>& other)
  {
    this->backend_ = std::make_shared<BackendType>(new ScalarType[other.size()], other.size());
    for (size_t ii = 0; ii < other.size(); ++ii)
      this->backend_->operator[](ii) = other[ii];
  }

  EigenMappedDenseVector(const std::initializer_list<ScalarType>& other)
  {
    this->backend_ = std::make_shared<BackendType>(new ScalarType[other.size()], other.size());
    size_t ii = 0;
    for (auto element : other) {
      this->backend_->operator[](ii) = element;
      ++ii;
    }
  }

  /**
   *  \brief  This constructor does not do a deep copy.
   */
  EigenMappedDenseVector(const ThisType& other)
  {
    this->backend_ = other.backend_;
  }

  /**
   * \brief This constructor does a deep copy.
   */
  EigenMappedDenseVector(const BackendType& other)
  {
    this->backend_          = std::make_shared<BackendType>(new ScalarType[other.size()], other.size());
    this->backend_->operator=(other);
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  EigenMappedDenseVector(BackendType* backend_ptr)
  {
    this->backend_ = std::shared_ptr<BackendType>(backend_ptr);
  }

  EigenMappedDenseVector(std::shared_ptr<BackendType> backend_ptr)
  {
    this->backend_ = backend_ptr;
  }

  /**
   * \brief does a deep copy;
   */
  ThisType& operator=(const BackendType& other)
  {
    this->backend_          = std::make_shared<BackendType>(new ScalarType[other.size()], other.size());
    this->backend_->operator=(other);
    return *this;
  }

  using VectorInterfaceType::add;
  using VectorInterfaceType::sub;
  using BaseType::backend;

private:
  inline void ensure_uniqueness() const
  {
    if (!this->backend_.unique()) {
      auto new_backend     = std::make_shared<BackendType>(new ScalarType[this->backend_->size()], this->backend_->size());
      new_backend->operator=(*(this->backend_));
      this->backend_       = new_backend;
    }
  } // ... ensure_uniqueness(...)

  friend class EigenBaseVector<EigenMappedDenseVectorTraits<ScalarType>>;
}; // class EigenMappedDenseVector


/**
 *  \brief Traits for EigenDenseMatrix.
 */
template <class ScalarImp = double>
class EigenDenseMatrixTraits
{
public:
  typedef ScalarImp ScalarType;
  typedef EigenDenseMatrix<ScalarType> derived_type;
  typedef typename ::Eigen::Matrix<ScalarType, ::Eigen::Dynamic, ::Eigen::Dynamic> BackendType;
}; // class DenseMatrixTraits


/**
 *  \brief  A dense metrix implementation of MatrixInterface using the eigen backend.
 */
template <class ScalarImp = double>
class EigenDenseMatrix : public MatrixInterface<EigenDenseMatrixTraits<ScalarImp>>,
                         public EigenMatrixInterfaceDynamic,
                         public ProvidesBackend<EigenDenseMatrixTraits<ScalarImp>>,
                         public ProvidesDataAccess<EigenDenseMatrixTraits<ScalarImp>>
{
  typedef EigenDenseMatrix<ScalarImp> ThisType;
  typedef MatrixInterface<EigenDenseMatrixTraits<ScalarImp>> MatrixInterfaceType;
  static_assert(!std::is_same<DUNE_STUFF_SSIZE_T, int>::value,
                "You have to manually disable the constructor below which uses DUNE_STUFF_SSIZE_T!");

public:
  typedef EigenDenseMatrixTraits<ScalarImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ScalarType ScalarType;

  EigenDenseMatrix(const size_t rr = 0, const size_t cc = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(rr, cc))
  {
    if (FloatCmp::eq(value, ScalarType(0)))
      backend_->setZero();
    else {
      backend_->setOnes();
      backend_->operator*=(value);
    }
  }

  /// This constructor is needed for the python bindings.
  EigenDenseMatrix(const DUNE_STUFF_SSIZE_T rr, const DUNE_STUFF_SSIZE_T cc = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(MatrixInterfaceType::assert_is_size_t_compatible_and_convert(rr),
                               MatrixInterfaceType::assert_is_size_t_compatible_and_convert(cc)))
  {
    if (FloatCmp::eq(value, ScalarType(0)))
      backend_->setZero();
    else {
      backend_->setOnes();
      backend_->operator*=(value);
    }
  }

  EigenDenseMatrix(const int rr, const int cc = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(MatrixInterfaceType::assert_is_size_t_compatible_and_convert(rr),
                               MatrixInterfaceType::assert_is_size_t_compatible_and_convert(cc)))
  {
    if (FloatCmp::eq(value, ScalarType(0)))
      backend_->setZero();
    else {
      backend_->setOnes();
      backend_->operator*=(value);
    }
  }

  /// This constructors ignores the given pattern and initializes the matrix with 0.
  EigenDenseMatrix(const size_t rr, const size_t cc, const SparsityPatternDefault& /*pattern*/)
    : backend_(new BackendType(MatrixInterfaceType::assert_is_size_t_compatible_and_convert(rr),
                               MatrixInterfaceType::assert_is_size_t_compatible_and_convert(cc)))
  {
    backend_->setZero();
  }

  EigenDenseMatrix(const ThisType& other) = default;

  EigenDenseMatrix(const BackendType& other)
    : backend_(new BackendType(other))
  {
  }

  template <class M>
  EigenDenseMatrix(const MatrixInterface<M>& other)
    : backend_(new BackendType(other.rows(), other.cols()))
  {
    for (size_t ii = 0; ii < other.rows(); ++ii)
      for (size_t jj = 0; jj < other.cols(); ++jj)
        set_entry(ii, jj, other.get_entry(ii, jj));
  } // EigenDenseMatrix(...)

  template <class T>
  EigenDenseMatrix(const DenseMatrix<T>& other)
    : backend_(new BackendType(other.rows(), other.cols()))
  {
    for (size_t ii = 0; ii < other.rows(); ++ii)
      for (size_t jj = 0; jj < other.cols(); ++jj)
        set_entry(ii, jj, other[ii][jj]);
  } // EigenDenseMatrix(...)

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  EigenDenseMatrix(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  EigenDenseMatrix(std::shared_ptr<BackendType> backend_ptr)
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
   * \defgroup data ´´These methods are required by the ProvidesDataAccess interface.``
   * \{
   */

  ScalarType* data()
  {
    return backend_->data();
  }

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
    if (!has_equal_shape(xx))
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The shape of xx (" << xx.rows() << "x" << xx.cols() << ") does not match the shape of this ("
                                     << rows()
                                     << "x"
                                     << cols()
                                     << ")!");
    const auto& xx_ref = *(xx.backend_);
    backend() += alpha * xx_ref;
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

  template <class T1, class T2>
  inline void mv(const EigenBaseVector<T1, ScalarType>& xx, EigenBaseVector<T2, ScalarType>& yy) const
  {
    yy.backend().transpose() = backend_->operator*(*xx.backend_);
  }

  void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(ii < rows());
    assert(jj < cols());
    backend()(ii, jj) += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(ii < rows());
    assert(jj < cols());
    backend()(ii, jj) = value;
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    assert(ii < rows());
    assert(jj < cols());
    return backend_->operator()(ii, jj);
  } // ... get_entry(...)

  void clear_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    ensure_uniqueness();
    for (size_t jj = 0; jj < cols(); ++jj)
      backend_->operator()(ii, jj) = ScalarType(0);
  } // ... clear_row(...)

  void clear_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    ensure_uniqueness();
    for (size_t ii = 0; ii < rows(); ++ii)
      backend_->operator()(ii, jj) = ScalarType(0);
  } // ... clear_col(...)

  void unit_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    ensure_uniqueness();
    for (size_t jj = 0; jj < cols(); ++jj)
      backend_->operator()(ii, jj) = ScalarType(0);
    backend_->operator()(ii, ii) = ScalarType(1);
  } // ... unit_row(...)

  void unit_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    ensure_uniqueness();
    for (size_t ii = 0; ii < rows(); ++ii)
      backend_->operator()(ii, jj) = ScalarType(0);
    backend_->operator()(jj, jj) = ScalarType(1);
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
}; // class EigenDenseMatrix


/**
 * \brief Traits for EigenRowMajorSparseMatrix.
 */
template <class ScalarImp = double>
class EigenRowMajorSparseMatrixTraits
{
public:
  typedef ScalarImp ScalarType;
  typedef EigenRowMajorSparseMatrix<ScalarType> derived_type;
  typedef typename ::Eigen::SparseMatrix<ScalarType, ::Eigen::RowMajor> BackendType;
}; // class RowMajorSparseMatrixTraits


/**
 * \brief A sparse matrix implementation of the MatrixInterface with row major memory layout.
 */
template <class ScalarImp = double>
class EigenRowMajorSparseMatrix : public MatrixInterface<EigenRowMajorSparseMatrixTraits<ScalarImp>>,
                                  public EigenMatrixInterfaceDynamic,
                                  public ProvidesBackend<EigenRowMajorSparseMatrixTraits<ScalarImp>>
{
  typedef EigenRowMajorSparseMatrix<ScalarImp> ThisType;
  typedef MatrixInterface<EigenRowMajorSparseMatrixTraits<ScalarImp>> MatrixInterfaceType;
  static_assert(!std::is_same<DUNE_STUFF_SSIZE_T, int>::value,
                "You have to manually disable the constructor below which uses DUNE_STUFF_SSIZE_T!");

public:
  typedef EigenRowMajorSparseMatrixTraits<ScalarImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ScalarType ScalarType;

  /**
   * \brief This is the constructor of interest which creates a sparse matrix.
   */
  EigenRowMajorSparseMatrix(const size_t rr, const size_t cc, const SparsityPatternDefault& pattern)
  {
    backend_ = std::make_shared<BackendType>(assert_is_IndexType_compatible_and_convert(rr),
                                             assert_is_IndexType_compatible_and_convert(cc));
    if (size_t(pattern.size()) != rr)
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of the pattern (" << pattern.size() << ") does not match the number of rows of this ("
                                             << rows()
                                             << ")!");
    for (size_t row = 0; row < size_t(pattern.size()); ++row) {
      backend_->startVec(assert_is_IndexType_compatible_and_convert(row));
      const auto& columns = pattern.inner(row);
      for (auto& column : columns) {
        backend_->insertBackByOuterInner(assert_is_IndexType_compatible_and_convert(row),
                                         assert_is_IndexType_compatible_and_convert(column));
      }
      // create diagonal entry (insertBackByOuterInner() can not handle empty rows)
      if (columns.size() == 0)
        backend_->insertBackByOuterInner(assert_is_IndexType_compatible_and_convert(row),
                                         assert_is_IndexType_compatible_and_convert(row));
    }
    backend_->finalize();
    backend_->makeCompressed();
  }

  EigenRowMajorSparseMatrix(const size_t rr = 0, const size_t cc = 0)
  {
    backend_ = std::make_shared<BackendType>(assert_is_IndexType_compatible_and_convert(rr),
                                             assert_is_IndexType_compatible_and_convert(cc));
  }

  /// This constructor is needed for the python bindings.
  EigenRowMajorSparseMatrix(const DUNE_STUFF_SSIZE_T rr, const DUNE_STUFF_SSIZE_T cc = 0)
    : backend_(new BackendType(MatrixInterfaceType::assert_is_size_t_compatible_and_convert(rr),
                               MatrixInterfaceType::assert_is_size_t_compatible_and_convert(cc)))
  {
  }

  EigenRowMajorSparseMatrix(const int rr, const int cc = 0)
    : EigenRowMajorSparseMatrix(MatrixInterfaceType::assert_is_size_t_compatible_and_convert(rr),
                                MatrixInterfaceType::assert_is_size_t_compatible_and_convert(cc))
  {
    //    backend_ = std::make_shared<BackendType>(assert_is_IndexType_compatible_and_convert(rr),
    //                                             assert_is_IndexType_compatible_and_convert(cc));
  }

  EigenRowMajorSparseMatrix(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  EigenRowMajorSparseMatrix(const BackendType& other)
    : backend_(new BackendType(other))
  {
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  EigenRowMajorSparseMatrix(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  EigenRowMajorSparseMatrix(std::shared_ptr<BackendType> backend_ptr)
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
    if (!has_equal_shape(xx))
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The shape of xx (" << xx.rows() << "x" << xx.cols() << ") does not match the shape of this ("
                                     << rows()
                                     << "x"
                                     << cols()
                                     << ")!");
    const auto& xx_ref = *(xx.backend_);
    backend() += alpha * xx_ref;
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

  template <class T1, class T2>
  inline void mv(const EigenBaseVector<T1, ScalarType>& xx, EigenBaseVector<T2, ScalarType>& yy) const
  {
    yy.backend().transpose() = backend_->operator*(*xx.backend_);
  }

  void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(these_are_valid_indices(ii, jj));
    backend().coeffRef(assert_is_IndexType_compatible_and_convert(ii),
                       assert_is_IndexType_compatible_and_convert(jj)) += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(these_are_valid_indices(ii, jj));
    backend().coeffRef(assert_is_IndexType_compatible_and_convert(ii), assert_is_IndexType_compatible_and_convert(jj)) =
        value;
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    assert(ii < rows());
    assert(jj < cols());
    return backend_->coeff(assert_is_IndexType_compatible_and_convert(ii),
                           assert_is_IndexType_compatible_and_convert(jj));
  } // ... get_entry(...)

  void clear_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    backend().row(assert_is_IndexType_compatible_and_convert(ii)) *= ScalarType(0);
  } // ... clear_row(...)

  void clear_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    ensure_uniqueness();
    for (size_t row = 0; assert_is_IndexType_compatible_and_convert(row) < backend_->outerSize(); ++row) {
      for (typename BackendType::InnerIterator row_it(*backend_, assert_is_IndexType_compatible_and_convert(row));
           row_it;
           ++row_it) {
        const size_t col = row_it.col();
        if (col == jj) {
          backend_->coeffRef(assert_is_IndexType_compatible_and_convert(row),
                             assert_is_IndexType_compatible_and_convert(jj)) = ScalarType(0);
          break;
        } else if (col > jj)
          break;
      }
    }
  } // ... clear_col(...)

  void unit_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    if (!these_are_valid_indices(ii, ii))
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Diagonal entry (" << ii << ", " << ii << ") is not contained in the sparsity pattern!");
    backend().row(assert_is_IndexType_compatible_and_convert(ii)) *= ScalarType(0);
    set_entry(ii, ii, ScalarType(1));
  } // ... unit_row(...)

  void unit_col(const size_t jj)
  {
    if (jj >= rows())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    ensure_uniqueness();
    for (size_t row = 0; assert_is_IndexType_compatible_and_convert(row) < backend_->outerSize(); ++row) {
      for (typename BackendType::InnerIterator row_it(*backend_, assert_is_IndexType_compatible_and_convert(row));
           row_it;
           ++row_it) {
        const size_t col = row_it.col();
        if (col == jj) {
          if (col == row)
            backend_->coeffRef(assert_is_IndexType_compatible_and_convert(row),
                               assert_is_IndexType_compatible_and_convert(col)) = ScalarType(1);
          else
            backend_->coeffRef(assert_is_IndexType_compatible_and_convert(row),
                               assert_is_IndexType_compatible_and_convert(jj)) = ScalarType(0);
          break;
        } else if (col > jj)
          break;
      }
    }
  } // ... unit_col(...)
  /**
   * \}
   */

private:
  typedef typename BackendType::Index IndexType;

  bool these_are_valid_indices(const size_t ii, const size_t jj) const
  {
    if (ii >= rows())
      return false;
    if (jj >= cols())
      return false;
    for (size_t row = ii; assert_is_IndexType_compatible_and_convert(row) < backend_->outerSize(); ++row) {
      for (typename BackendType::InnerIterator row_it(*backend_, assert_is_IndexType_compatible_and_convert(row));
           row_it;
           ++row_it) {
        const size_t col = row_it.col();
        if ((ii == row) && (jj == col))
          return true;
        else if ((row > ii) && (col > jj))
          return false;
      }
    }
    return false;
  } // ... these_are_valid_indices(...)

  template <class SizeType>
  static inline IndexType assert_is_IndexType_compatible_and_convert(const SizeType& size)
  {
    assert(size < std::numeric_limits<IndexType>::max());
    return IndexType(size);
  } // ... assert_is_IndexType_compatible_and_convert(...)

  inline void ensure_uniqueness() const
  {
    if (!backend_.unique())
      backend_ = std::make_shared<BackendType>(*backend_);
  } // ... ensure_uniqueness(...)

  mutable std::shared_ptr<BackendType> backend_;
}; // class EigenRowMajorSparseMatrix


#else // HAVE_EIGEN


template <class Traits, class ScalarImp>
class EigenBaseVector
{
  static_assert(Dune::AlwaysFalse<ScalarImp>::value, "You are missing Eigen!");
};

template <class ScalarImp>
class EigenDenseVector
{
  static_assert(Dune::AlwaysFalse<ScalarImp>::value, "You are missing Eigen!");
};

template <class ScalarImp>
class EigenMappedDenseVector
{
  static_assert(Dune::AlwaysFalse<ScalarImp>::value, "You are missing Eigen!");
};

template <class ScalarImp>
class EigenDenseMatrix
{
  static_assert(Dune::AlwaysFalse<ScalarImp>::value, "You are missing Eigen!");
};

template <class ScalarImp>
class EigenRowMajorSparseMatrix
{
  static_assert(Dune::AlwaysFalse<ScalarImp>::value, "You are missing Eigen!");
};


#endif // HAVE_EIGEN

} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_EIGEN_HH
