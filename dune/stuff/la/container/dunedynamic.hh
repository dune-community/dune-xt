// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright Holders: Felix Albrecht, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_CONTAINER_DYNAMICVECTOR_HH
#define DUNE_STUFF_LA_CONTAINER_DYNAMICVECTOR_HH

#include <memory>
#include <type_traits>

#include <dune/common/dynvector.hh>
#include <dune/stuff/common/disable_warnings.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/stuff/common/reenable_warnings.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/typetraits.hh>

#include "interfaces.hh"
#include "pattern.hh"

namespace Dune {
namespace Pymor {
namespace Operators {

// forwards, needed for friendlyness
template <class ScalarImp>
class DuneDynamic;

template <class ScalarImp>
class DuneDynamicInverse;

} // namespace Operators
} // namespace Pymor

namespace Stuff {
namespace LA {


// forwards
template <class ScalarImp>
class DuneDynamicVector;

template <class ScalarImp>
class DuneDynamicMatrix;


/// Traits for DuneDynamicVector
template <class ScalarImp = double>
class DuneDynamicVectorTraits
{
public:
  typedef ScalarImp ScalarType;
  typedef DuneDynamicVector<ScalarType> derived_type;
  typedef Dune::DynamicVector<ScalarType> BackendType;
};


/**
 *  \brief A dense vector implementation of VectorInterface using the Dune::DynamicVector.
 */
template <class ScalarImp = double>
class DuneDynamicVector : public VectorInterface<DuneDynamicVectorTraits<ScalarImp>>,
                          public ProvidesBackend<DuneDynamicVectorTraits<ScalarImp>>
{
  typedef DuneDynamicVector<ScalarImp> ThisType;
  typedef VectorInterface<DuneDynamicVectorTraits<ScalarImp>> VectorInterfaceType;
  static_assert(!std::is_same<DUNE_STUFF_SSIZE_T, int>::value,
                "You have to manually disable the constructor below which uses DUNE_STUFF_SSIZE_T!");

public:
  typedef DuneDynamicVectorTraits<ScalarImp> Traits;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::BackendType BackendType;

  DuneDynamicVector(const size_t ss = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(ss, value))
  {
  }

  /// This constructor is needed for the python bindings.
  DuneDynamicVector(const DUNE_STUFF_SSIZE_T ss, const ScalarType value = ScalarType(0))
    : DuneDynamicVector(VectorInterfaceType::assert_is_size_t_compatible_and_convert(ss), value)
  {
  }

  /// This constructor is needed because marking the above one as explicit had no effect.
  DuneDynamicVector(const int ss, const ScalarType value = ScalarType(0))
    : DuneDynamicVector(VectorInterfaceType::assert_is_size_t_compatible_and_convert(ss), value)
  {
  }

  DuneDynamicVector(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  DuneDynamicVector(const BackendType& other)
    : backend_(new BackendType(other))
  {
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  DuneDynamicVector(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  DuneDynamicVector(std::shared_ptr<BackendType> backend_ptr)
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
      DUNE_THROW_COLORFULLY(Exception::shapes_do_not_match,
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
    return backend_->operator[](ii);
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

  virtual ScalarType dot(const ThisType& other) const DS_OVERRIDE
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exception::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    return backend_->operator*(*(other.backend_));
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
      DUNE_THROW_COLORFULLY(Exception::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    if (result.size() != size())
      DUNE_THROW_COLORFULLY(Exception::shapes_do_not_match,
                            "The size of result (" << result.size() << ") does not match the size of this (" << size()
                                                   << ")!");
    for (size_t ii = 0; ii < size(); ++ii)
      result.backend_->operator[](ii) = backend_->operator[](ii) + other.backend_->operator[](ii);
  } // ... add(...)

  virtual void iadd(const ThisType& other) DS_OVERRIDE
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exception::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    ensure_uniqueness();
    backend_->operator+=(*(other.backend_));
  } // ... iadd(...)

  virtual void sub(const ThisType& other, ThisType& result) const DS_OVERRIDE
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exception::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    if (result.size() != size())
      DUNE_THROW_COLORFULLY(Exception::shapes_do_not_match,
                            "The size of result (" << result.size() << ") does not match the size of this (" << size()
                                                   << ")!");
    for (size_t ii = 0; ii < size(); ++ii)
      result.backend_->operator[](ii) = backend_->operator[](ii)-other.backend_->operator[](ii);
  } // ... sub(...)

  virtual void isub(const ThisType& other) DS_OVERRIDE
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exception::shapes_do_not_match,
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
  inline void ensure_uniqueness()
  {
    if (!backend_.unique())
      backend_ = std::make_shared<BackendType>(*backend_);
  } // ... ensure_uniqueness(...)

  friend class VectorInterface<DuneDynamicVectorTraits<ScalarType>>;
  friend class DuneDynamicMatrix<ScalarType>;
  friend class Dune::Pymor::Operators::DuneDynamicInverse<ScalarType>;
  friend class Dune::Pymor::Operators::DuneDynamic<ScalarType>;

  std::shared_ptr<BackendType> backend_;
}; // class DuneDynamicVector


template <class ScalarImp = double>
class DuneDynamicMatrixTraits
{
public:
  typedef ScalarImp ScalarType;
  typedef DuneDynamicMatrix<ScalarType> derived_type;
  typedef Dune::DynamicMatrix<ScalarType> BackendType;
};


/**
 *  \brief  A dense metrix implementation of MatrixInterface using the eigen backend.
 */
template <class ScalarImp = double>
class DuneDynamicMatrix : public MatrixInterface<DuneDynamicMatrixTraits<ScalarImp>>,
                          public ProvidesBackend<DuneDynamicMatrixTraits<ScalarImp>>
{
  typedef DuneDynamicMatrix<ScalarImp> ThisType;
  typedef MatrixInterface<DuneDynamicMatrixTraits<ScalarImp>> MatrixInterfaceType;
  static_assert(!std::is_same<DUNE_STUFF_SSIZE_T, int>::value,
                "You have to manually disable the constructor below which uses DUNE_STUFF_SSIZE_T!");

public:
  typedef DuneDynamicMatrixTraits<ScalarImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ScalarType ScalarType;

  DuneDynamicMatrix(const size_t rr = 0, const size_t cc = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(rr, cc, value))
  {
  }

  /// This constructor is needed for the python bindings.
  DuneDynamicMatrix(const DUNE_STUFF_SSIZE_T rr, const DUNE_STUFF_SSIZE_T cc = 0,
                    const ScalarType value = ScalarType(0))
    : DuneDynamicMatrix(MatrixInterfaceType::assert_is_size_t_compatible_and_convert(rr),
                        MatrixInterfaceType::assert_is_size_t_compatible_and_convert(cc), value)
  {
  }

  /// This constructor is needed because marking the above one as explicit had no effect.
  DuneDynamicMatrix(const int rr, const int cc = 0, const ScalarType value = ScalarType(0))
    : DuneDynamicMatrix(MatrixInterfaceType::assert_is_size_t_compatible_and_convert(rr),
                        MatrixInterfaceType::assert_is_size_t_compatible_and_convert(cc), value)
  {
  }

  /// This constructors ignores the given pattern and initializes the matrix with 0.
  DuneDynamicMatrix(const size_t rr, const size_t cc, const SparsityPatternDefault& /*pattern*/)
    : DuneDynamicMatrix(rr, cc)
  {
  }

  DuneDynamicMatrix(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  DuneDynamicMatrix(const BackendType& other)
    : backend_(new BackendType(other))
  {
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  DuneDynamicMatrix(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  DuneDynamicMatrix(std::shared_ptr<BackendType> backend_ptr)
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
          Exception::shapes_do_not_match,
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

  template <class SourceType, class RangeType>
  inline void mv(const SourceType& /*xx*/, RangeType& /*yy*/) const
  {
    static_assert(Dune::AlwaysFalse<SourceType>::value, "Not available for this combination of xx and yy!");
  }

  inline void mv(const DuneDynamicVector<ScalarType>& xx, DuneDynamicVector<ScalarType>& yy) const
  {
    backend_->mv(xx.backend(), yy.backend());
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
      DUNE_THROW_COLORFULLY(Exception::index_out_of_range,
                            "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    backend_->operator[](ii) *= ScalarType(0);
  } // ... clear_row(...)

  void clear_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW_COLORFULLY(Exception::index_out_of_range,
                            "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    for (size_t ii = 0; ii < rows(); ++ii)
      backend_->operator[](ii)[jj] = ScalarType(0);
  } // ... clear_col(...)

  void unit_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW_COLORFULLY(Exception::index_out_of_range,
                            "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    auto& row = backend_->operator[](ii);
    for (size_t jj = 0; jj < cols(); ++jj)
      row[jj] = ScalarType(0);
    row[ii]   = ScalarType(1);
  } // ... unit_row(...)

  void unit_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW_COLORFULLY(Exception::index_out_of_range,
                            "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    for (size_t ii = 0; ii < rows(); ++ii)
      backend_->operator[](ii)[jj] = ScalarType(0);
    backend_->operator[](jj)[jj] = ScalarType(1);
  } // ... unit_col(...)

  /**
   * \}
   */

private:
  inline void ensure_uniqueness()
  {
    if (!backend_.unique())
      backend_ = std::make_shared<BackendType>(*backend_);
  } // ... ensure_uniqueness(...)

  friend class Dune::Pymor::Operators::DuneDynamicInverse<ScalarType>;
  friend class Dune::Pymor::Operators::DuneDynamic<ScalarType>;

  std::shared_ptr<BackendType> backend_;
}; // class DuneDynamicMatrix


} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_DYNAMICVECTOR_HH
