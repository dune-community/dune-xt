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

#ifndef DUNE_XT_LA_CONTAINER_VECTOR_VIEW_HH
#define DUNE_XT_LA_CONTAINER_VECTOR_VIEW_HH

#include <dune/xt/common/exceptions.hh>

#include "vector-interface.hh"

namespace Dune {
namespace XT {
namespace LA {


// forwards
template <class VectorImp>
class ConstVectorView;

template <class VectorImp>
class VectorView;


namespace internal {


template <class VectorImp>
class ConstVectorViewTraits
  : public VectorTraitsBase<typename VectorImp::ScalarType,
                            ConstVectorView<VectorImp>,
                            typename VectorImp::Traits::BackendType,
                            VectorImp::Traits::backend_type,
                            VectorImp::Traits::dense_matrix_type,
                            VectorImp::Traits::sparse_matrix_type>
{};

template <class VectorImp>
class VectorViewTraits
  : public VectorTraitsBase<typename VectorImp::ScalarType,
                            VectorView<VectorImp>,
                            typename VectorImp::Traits::BackendType,
                            VectorImp::Traits::backend_type,
                            VectorImp::Traits::dense_matrix_type,
                            VectorImp::Traits::sparse_matrix_type>
{};

template <class VectorImp>
VectorImp& empty_vector_ref()
{
  static VectorImp vector_;
  return vector_;
}


} // namespace internal


template <class VectorImp>
class ConstVectorView
  : public VectorInterface<internal::ConstVectorViewTraits<VectorImp>, typename VectorImp::ScalarType>
{
  using BaseType = VectorInterface<internal::ConstVectorViewTraits<VectorImp>, typename VectorImp::ScalarType>;
  using ThisType = ConstVectorView;

public:
  using ScalarType = typename BaseType::ScalarType;
  using RealType = typename BaseType::RealType;
  using Vector = VectorImp;

  // This constructor is only here for the interface to compile
  explicit ConstVectorView(const size_t /*sz*/ = 0,
                           const ScalarType /*value*/ = ScalarType(0),
                           const size_t /*num_mutexes*/ = 1)
    : vector_(internal::empty_vector_ref<VectorImp>())
    , first_entry_(0)
    , past_last_entry_(0)
  {
    DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "This constructor does not make sense for VectorView");
  }

  // This is the actual constructor
  ConstVectorView(const Vector& vector, const size_t first_entry, const size_t past_last_entry)
    : vector_(vector)
    , first_entry_(first_entry)
    , past_last_entry_(past_last_entry)
  {}

  size_t index(const size_t ii) const
  {
    return first_entry_ + ii;
  }

  inline size_t size() const
  {
    return past_last_entry_ - first_entry_;
  }

  inline void resize(const size_t /*new_size*/)
  {
    DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "You cannot use non-const methods on ConstVectorView");
  }

  inline void add_to_entry(const size_t /*ii*/, const ScalarType& /*value*/)
  {
    DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "You cannot use non-const methods on ConstVectorView");
  }

  inline void set_entry(const size_t /*ii*/, const ScalarType& /*value*/)
  {
    DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "You cannot use non-const methods on ConstVectorView");
  }

  inline ScalarType get_entry(const size_t ii) const
  {
    return vector_.get_entry(index(ii));
  }

  inline ScalarType& operator[](const size_t /*ii*/)
  {
    DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "You cannot use non-const methods on ConstVectorView");
    static ScalarType ret;
    return ret;
  }

  inline const ScalarType& operator[](const size_t ii) const
  {
    return vector_[index(ii)];
  }

  inline void scal(const ScalarType& /*alpha*/)
  {
    DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "You cannot use non-const methods on ConstVectorView");
  }

  inline void axpy(const ScalarType& /*alpha*/, const ThisType& /*xx*/)
  {
    DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "You cannot use non-const methods on ConstVectorView");
  }

  inline void add_to_entry(const size_t /*ii*/, const size_t /*jj*/, const ScalarType& /*value*/)
  {
    DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "You cannot use non-const methods on ConstVectorView");
  }

  inline void set_entry(const size_t /*ii*/, const size_t /*jj*/, const ScalarType& /*value*/)
  {
    DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "You cannot use non-const methods on ConstVectorView");
  }

  const Vector& vec() const
  {
    return vector_;
  }

  size_t first_entry() const
  {
    return first_entry_;
  }

  size_t past_last_entry() const
  {
    return past_last_entry_;
  }

protected:
  inline const ScalarType& get_unchecked_ref(const size_t ii) const
  {
    return vector_[first_entry_ + ii];
  }

private:
  friend class VectorInterface<internal::ConstVectorViewTraits<VectorImp>, ScalarType>;

  const Vector& vector_;
  const size_t first_entry_;
  const size_t past_last_entry_;
}; // class ConstVectorView


template <class VectorImp>
class VectorView : public VectorInterface<internal::VectorViewTraits<VectorImp>, typename VectorImp::ScalarType>
{
  using BaseType = VectorInterface<internal::VectorViewTraits<VectorImp>, typename VectorImp::ScalarType>;
  using ConstVectorViewType = ConstVectorView<VectorImp>;
  using ThisType = VectorView;

public:
  using ScalarType = typename BaseType::ScalarType;
  using RealType = typename BaseType::RealType;
  using Vector = VectorImp;

  // This constructor is only here for the interface to compile
  explicit VectorView(const size_t /*sz*/ = 0,
                      const ScalarType /*value*/ = ScalarType(0),
                      const size_t /*num_mutexes*/ = 1)
    : const_vector_view_()
    , vector_(internal::empty_vector_ref<VectorImp>())
  {
    DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "This constructor does not make sense for VectorView");
  }

  // This is the actual constructor
  VectorView(Vector& vector, const size_t first_entry, const size_t past_last_entry)
    : const_vector_view_(vector, first_entry, past_last_entry)
    , vector_(vector)
  {}

  size_t index(const size_t ii) const
  {
    return const_vector_view_.index(ii);
  }

  inline size_t size() const
  {
    return const_vector_view_.size();
  }

  inline void scal(const ScalarType& alpha)
  {
    for (size_t ii = 0; ii < size(); ++ii)
      set_entry(ii, get_entry(ii) * alpha);
  }

  inline void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    for (size_t ii = 0; ii < size(); ++ii)
      add_to_entry(ii, alpha * xx[ii]);
  }

  inline void add_to_entry(const size_t ii, const ScalarType& value)
  {
    assert(ii < size());
    vector_.add_to_entry(index(ii), value);
  }

  inline void set_entry(const size_t ii, const ScalarType& value)
  {
    assert(ii < size());
    vector_.set_entry(index(ii), value);
  }

  inline ScalarType get_entry(const size_t ii) const
  {
    return const_vector_view_.get_entry(ii);
  }

  inline ScalarType& operator[](const size_t ii)
  {
    return vector_[index(ii)];
  }

  inline const ScalarType& operator[](const size_t ii) const
  {
    return vector_[index(ii)];
  }

  template <class T>
  ThisType& operator+=(const VectorInterface<T, ScalarType>& other)
  {
    iadd(other);
    return *this;
  }

  template <class T>
  ThisType& operator-=(const VectorInterface<T, ScalarType>& other)
  {
    isub(other);
    return *this;
  }

  template <class T>
  void iadd(const VectorInterface<T, ScalarType>& other)
  {
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    for (size_t ii = 0; ii < size(); ++ii)
      add_to_entry(ii, other.get_entry(ii));
  } // ... iadd(...)

  template <class T>
  void isub(const VectorInterface<T, ScalarType>& other)
  {
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    for (size_t ii = 0; ii < size(); ++ii)
      add_to_entry(ii, -other.get_entry(ii));
  } // ... isub(...)

protected:
  inline ScalarType& get_unchecked_ref(const size_t ii)
  {
    return vector_[const_vector_view_.first_entry() + ii];
  }

  inline const ScalarType& get_unchecked_ref(const size_t ii) const
  {
    return vector_[const_vector_view_.first_entry() + ii];
  }

private:
  friend class VectorInterface<internal::VectorViewTraits<VectorImp>, ScalarType>;

  ConstVectorView<VectorImp> const_vector_view_;
  Vector& vector_;
}; // class VectorView


} // namespace LA
namespace Common {


template <class VectorImp>
struct VectorAbstraction<LA::ConstVectorView<VectorImp>>
  : public LA::internal::VectorAbstractionBase<LA::ConstVectorView<VectorImp>>
{
  static const bool is_contiguous = LA::internal::VectorAbstractionBase<VectorImp>::is_contiguous;
};


template <class VectorImp>
struct VectorAbstraction<LA::VectorView<VectorImp>>
  : public LA::internal::VectorAbstractionBase<LA::VectorView<VectorImp>>
{
  static const bool is_contiguous = LA::internal::VectorAbstractionBase<VectorImp>::is_contiguous;
};


} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_VECTOR_VIEW_HH
