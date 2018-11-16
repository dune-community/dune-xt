// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2017 - 2018)

#ifndef DUNE_XT_LA_CONTAINER_COMMON_VECTOR_SPARSE_HH
#define DUNE_XT_LA_CONTAINER_COMMON_VECTOR_SPARSE_HH

#include <cmath>
#include <initializer_list>
#include <memory>
#include <type_traits>
#include <vector>
#include <complex>
#include <mutex>
#include <numeric>

#include <dune/common/dynvector.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/unused.hh>

#include <dune/xt/common/exceptions.hh>

#include <dune/xt/la/container/vector-interface.hh>

namespace Dune {
namespace XT {
namespace LA {


// forwards
template <class ScalarImp>
class CommonSparseVector;


namespace internal {


template <class ScalarImp>
struct CommonSparseVectorTraits : VectorTraitsBase<ScalarImp,
                                                   CommonSparseVector<ScalarImp>,
                                                   void,
                                                   Backends::common_dense,
                                                   Backends::common_dense,
                                                   Backends::common_sparse>
{
  using EntriesVectorType = std::vector<ScalarImp>;
  using IndicesVectorType = std::vector<size_t>;
};


} // namespace internal


/**
 *  \brief A sparse vector implementation of VectorInterface
 */
template <class ScalarImp = double>
class CommonSparseVector : public VectorInterface<internal::CommonSparseVectorTraits<ScalarImp>, ScalarImp>
{
  using ThisType = CommonSparseVector;
  using InterfaceType = VectorInterface<internal::CommonSparseVectorTraits<ScalarImp>, ScalarImp>;

public:
  using typename InterfaceType::RealType;
  using typename InterfaceType::ScalarType;
  using typename InterfaceType::Traits;
  using IndicesVectorType = typename Traits::IndicesVectorType;
  using EntriesVectorType = typename Traits::EntriesVectorType;
  // needed to fix gcc compilation error due to ambiguous lookup of derived type
  using derived_type = typename Traits::derived_type;

private:
  using MutexesType = typename Traits::MutexesType;

public:
  explicit CommonSparseVector(const size_t sz = 0, const size_t num_mutexes = 1)
    : size_(sz)
    , entries_(new EntriesVectorType())
    , indices_(new IndicesVectorType())
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
  }

  explicit CommonSparseVector(const size_t sz, const ScalarType value, const size_t num_mutexes = 1)
    : size_(sz)
    , entries_(new EntriesVectorType(size_, value))
    , indices_(new IndicesVectorType(size_))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    for (size_t ii = 0; ii < size_; ++ii)
      (*indices_)[ii] = ii;
  }

  explicit CommonSparseVector(const std::vector<ScalarType>& other, const size_t num_mutexes = 1)
    : size_(other.size())
    , entries_(new EntriesVectorType(size_))
    , indices_(new IndicesVectorType(size_))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    for (size_t ii = 0; ii < size_; ++ii) {
      (*entries_)[ii] = other[ii];
      (*indices_)[ii] = ii;
    }
  }

  template <class VectorType>
  explicit CommonSparseVector(
      const typename std::enable_if_t<XT::Common::VectorAbstraction<VectorType>::is_vector, size_t> sz,
      const VectorType& entries,
      const IndicesVectorType& indices,
      const size_t num_mutexes = 1)
    : size_(sz)
    , entries_(new EntriesVectorType(entries.size()))
    , indices_(new IndicesVectorType(indices.size()))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    assert(entries.size() == indices.size());
    for (size_t ii = 0; ii < entries.size(); ++ii) {
      (*entries_)[ii] = entries[ii];
      (*indices_)[ii] = indices[ii];
    }
  }

  explicit CommonSparseVector(const std::initializer_list<ScalarType>& other, const size_t num_mutexes = 1)
    : size_(other.size())
    , entries_(new EntriesVectorType(size_))
    , indices_(new IndicesVectorType(size_))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    size_t ii = 0;
    for (auto element : other) {
      entries_->operator[](ii) = element;
      (*indices_)[ii] = ii;
      ++ii;
    }
  }

  CommonSparseVector(const ThisType& other)
    : size_(other.size_)
    , entries_(std::make_shared<EntriesVectorType>(*other.entries_))
    , indices_(std::make_shared<IndicesVectorType>(*other.indices_))
    , mutexes_(std::make_unique<MutexesType>(other.mutexes_->size()))
  {
  }

  template <class OtherVectorType>
  explicit CommonSparseVector(
      const OtherVectorType& other,
      const typename std::enable_if_t<XT::Common::VectorAbstraction<OtherVectorType>::is_vector, bool> prune = true,
      const ScalarType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value(),
      const size_t num_mutexes = 1)
    : size_(other.size())
    , entries_(new EntriesVectorType())
    , indices_(new IndicesVectorType())
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    for (size_t ii = 0; ii < size_; ++ii) {
      const ScalarType& value = XT::Common::VectorAbstraction<OtherVectorType>::get_entry(other, ii);
      if (!prune || XT::Common::FloatCmp::ne(value, ScalarType(0), eps)) {
        entries_->push_back(value);
        indices_->push_back(ii);
      }
    } // ii
  }

  template <class OtherVectorType>
  explicit CommonSparseVector(
      const OtherVectorType& other,
      const typename std::enable_if_t<XT::Common::VectorAbstraction<OtherVectorType>::is_vector, bool> prune,
      const size_t num_mutexes)
    : CommonSparseVector(other, prune, Common::FloatCmp::DefaultEpsilon<ScalarType>::value(), num_mutexes)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      size_ = other.size_;
      *entries_ = *other.entries_;
      *indices_ = *other.indices_;
      mutexes_ = std::make_unique<MutexesType>(other.mutexes_->size());
    }
    return *this;
  }

  template <int size>
  ThisType& operator=(const Dune::FieldVector<ScalarType, size>& other)
  {
    clear();
    size_ = size;
    for (size_t kk = 0; kk < size; ++kk)
      if (XT::Common::FloatCmp::ne(other[kk], 0.)) {
        entries_->push_back(other[kk]);
        indices_->push_back(kk);
      }
  }

  void deep_copy(const ThisType& other)
  {
    size_ = other.size_;
    *entries_ = *other.entries_;
    *indices_ = *other.indices_;
  }

  ThisType& operator=(const ScalarType& value)
  {
    for (auto& element : *entries_)
      element = value;
    return *this;
  } // ... operator=(...)

  /// \name Required by ContainerInterface.
  /// \{

  ThisType copy() const
  {
    return ThisType(size_, *entries_, *indices_, mutexes_->size());
  }

  void scal(const ScalarType& alpha)
  {
    const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
    auto& entries_ref = *entries_;
    for (size_t ii = 0; ii < entries_ref.size(); ++ii)
      entries_ref[ii] *= alpha;
  }

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    const auto& xx_entries = *xx.entries_;
    const auto& xx_indices = *xx.indices_;
    if (xx.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of x (" << xx.size() << ") does not match the size of this (" << size() << ")!");
    for (size_t ii = 0; ii < xx_indices.size(); ++ii)
      add_to_entry(xx_indices[ii], alpha * xx_entries[ii]);
  } // ... axpy(...)

  bool has_equal_shape(const ThisType& other) const
  {
    return size() == other.size();
  }

  /// \}
  /// \name Required by VectorInterface.
  /// \{

  inline size_t size() const
  {
    return size_;
  }

  inline void resize(const size_t new_size)
  {
    if (new_size != size()) {
      EntriesVectorType new_entries;
      IndicesVectorType new_indices;
      for (size_t ii = 0; ii < indices_->size(); ++ii) {
        if ((*indices_)[ii] < new_size) {
          new_entries.push_back((*entries_)[ii]);
          new_indices.push_back((*indices_)[ii]);
        }
      }
      *entries_ = new_entries;
      *indices_ = new_indices;
      size_ = new_size;
    }
  }

  void add_to_entry(const size_t ii, const ScalarType& value)
  {
    internal::LockGuard DUNE_UNUSED(lock)(*mutexes_, ii, size());
    get_unchecked_ref(ii) += value;
  }

  void set_entry(const size_t ii, const ScalarType& value)
  {
    get_unchecked_ref(ii) = value;
  }

  // set entry without checking if entry already exists
  void set_new_entry(const size_t ii, const ScalarType& value, bool front = false)
  {
    if (front) {
      indices_->insert(indices_->begin(), ii);
      entries_->insert(entries_->begin(), value);
    } else {
      indices_->push_back(ii);
      entries_->push_back(value);
    }
  }

  ScalarType get_entry(const size_t ii) const
  {
    auto it = std::lower_bound(indices_->begin(), indices_->end(), ii);
    if (it == indices_->end() || *it != ii)
      return ScalarType(0.);
    else
      return (*entries_)[std::distance(indices_->begin(), it)];
  }

protected:
  inline ScalarType& get_unchecked_ref(const size_t ii)
  {
    auto it = std::lower_bound(indices_->begin(), indices_->end(), ii);
    size_t index = std::distance(indices_->begin(), it);
    if (it == indices_->end()) {
      indices_->push_back(ii);
      entries_->push_back(0.);
      return entries_->back();
    } else if (*it != ii) {
      indices_->insert(it, ii);
      entries_->insert(entries_->begin() + index, 0.);
      return (*entries_)[index];
    } else {
      return (*entries_)[index];
    }
  }

  inline const ScalarType& get_unchecked_ref(const size_t ii) const
  {
    auto it = std::lower_bound(indices_->begin(), indices_->end(), ii);
    size_t index = std::distance(indices_->begin(), it);
    if (it == indices_->end()) {
      indices_->push_back(ii);
      entries_->push_back(0.);
      return entries_->back();
    } else if (*it != ii) {
      indices_->insert(it, ii);
      entries_->insert(entries_->begin() + index, 0.);
      return (*entries_)[index];
    } else {
      return (*entries_)[index];
    }
  }

public:
  inline ScalarType& operator[](const size_t ii)
  {
    return get_unchecked_ref(ii);
  }

  inline const ScalarType& operator[](const size_t ii) const
  {
    return get_unchecked_ref(ii);
  }

  /// \}

  EntriesVectorType& entries()
  {
    return *entries_;
  }

  const EntriesVectorType& entries() const
  {
    return *entries_;
  }

  IndicesVectorType& indices()
  {
    return *indices_;
  }

  const IndicesVectorType& indices() const
  {
    return *indices_;
  }

  void clear()
  {
    entries_->clear();
    indices_->clear();
  }

  /// \name These methods override default implementations from VectorInterface.
  /// \{

  virtual ScalarType dot(const ThisType& other) const override final
  {
    assert(other.size() == size());
    const auto& indices = *indices_;
    const auto& other_indices = *other.indices_;
    const auto& entries = *entries_;
    const auto& other_entries = *other.entries_;
    ScalarType ret(0);
    size_t kk = 0;
    size_t other_kk = 0;
    while (kk < indices.size() && other_kk < other_indices.size()) {
      if (indices[kk] < other_indices[other_kk])
        ++kk;
      else if (indices[kk] > other_indices[other_kk])
        ++other_kk;
      else
        ret += entries[kk++] * other_entries[other_kk++];
    } // kk
    return ret;
  } // ... dot(...)

  virtual RealType l1_norm() const override final
  {
    return std::accumulate(entries_->begin(),
                           entries_->end(),
                           RealType(0.),
                           [](const RealType& a, const ScalarType& b) { return a + std::abs(b); });
  }

  virtual RealType l2_norm() const override final
  {
    return std::sqrt(
        std::accumulate(entries_->begin(), entries_->end(), RealType(0.), [](const RealType& a, const ScalarType& b) {
          return a + std::pow(std::abs(b), 2);
        }));
  }

  virtual RealType sup_norm() const override final
  {
    auto it = std::max_element(entries_->begin(), entries_->end(), [](const ScalarType& a, const ScalarType& b) {
      return std::abs(a) < std::abs(b);
    });
    return entries_->size() > 0 ? std::abs(*it) : 0.;
  }

  virtual ThisType add(const ThisType& other) const override final
  {
    ThisType ret = this->copy();
    ret.entries_->reserve(ret.entries_->size() + other.entries_->size());
    for (size_t kk = 0; kk < other.indices_->size(); ++kk)
      ret.add_to_entry(other.indices()[kk], other.entries()[kk]);
    return ret;
  } // ... add(...)

  virtual void add(const ThisType& other, ThisType& result) const override final
  {
    result.deep_copy(add(other));
  } // ... add(...)

  virtual ThisType sub(const ThisType& other) const override final
  {
    ThisType ret = this->copy();
    ret.entries_->reserve(ret.entries_->size() + other.entries_->size());
    for (size_t kk = 0; kk < other.indices_->size(); ++kk)
      ret.add_to_entry(other.indices()[kk], -other.entries()[kk]);
    return ret;
  } // ... sub(...)

  virtual void sub(const ThisType& other, ThisType& result) const override final
  {
    result.deep_copy(sub(other));
  } // ... add(...)

  virtual void iadd(const ThisType& other) override final
  {
    axpy(ScalarType(1), other);
  } // ... iadd(...)

  virtual void isub(const ThisType& other) override final
  {
    axpy(ScalarType(-1), other);
  } // ... isub(...)

  /// \}

  // without these using declarations, the free operator+/* function in xt/common/vector.hh is chosen instead of the
  // member function
  using InterfaceType::operator+;
  using InterfaceType::operator-;
  using InterfaceType::operator*;

private:
  friend class VectorInterface<internal::CommonSparseVectorTraits<ScalarType>, ScalarType>;

  size_t size_;
  std::shared_ptr<EntriesVectorType> entries_;
  std::shared_ptr<IndicesVectorType> indices_;
  std::unique_ptr<MutexesType> mutexes_;
}; // class CommonSparseVector


} // namespace LA
namespace Common {


template <class T>
struct VectorAbstraction<LA::CommonSparseVector<T>>
    : public LA::internal::VectorAbstractionBase<LA::CommonSparseVector<T>>
{
  static const bool contiguous = false;
};


} // namespace Common
} // namespace XT


template <class ScalarType, int size>
FieldVector<ScalarType, size>& operator+=(FieldVector<ScalarType, size>& lhs,
                                          const XT::LA::CommonSparseVector<ScalarType>& rhs)
{
  const auto& indices = rhs.indices();
  const auto& entries = rhs.entries();
  for (size_t kk = 0; kk < indices.size(); ++kk)
    lhs[indices[kk]] += entries[kk];
  return lhs;
}


} // namespace Dune


// begin: this is what we need for the lib
#if DUNE_XT_WITH_PYTHON_BINDINGS


extern template class Dune::XT::LA::CommonSparseVector<double>;


#endif // DUNE_XT_WITH_PYTHON_BINDINGS
// end: this is what we need for the lib


#endif // DUNE_XT_LA_CONTAINER_COMMON_VECTOR_SPARSE_HH
