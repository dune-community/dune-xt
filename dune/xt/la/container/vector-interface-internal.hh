// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2014 - 2016, 2018)
//   Tobias Leibner  (2014, 2017)

#ifndef DUNE_XT_LA_CONTAINER_VECTOR_INTERFACE_INTERNAL_HH
#define DUNE_XT_LA_CONTAINER_VECTOR_INTERFACE_INTERNAL_HH

#include <iterator>
#include <type_traits>

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/crtp.hh>
#include <dune/xt/common/exceptions.hh>

namespace Dune {
namespace XT {
namespace LA {

// forward
template <class Traits, class ScalarImp>
class VectorInterface;

namespace internal {

template <class Traits, class ScalarImp>
class VectorInputIterator : public std::iterator<std::input_iterator_tag, typename Traits::ScalarType>
{
  typedef VectorInputIterator<Traits, ScalarImp> ThisType;

public:
  typedef VectorInterface<Traits, ScalarImp> VectorType;
  typedef typename VectorType::ScalarType ScalarType;

public:
  explicit VectorInputIterator(const VectorType& vec, const bool end = false)
    : const_vec_(vec)
    , position_(0)
    , end_(end)
  {}

  ThisType& operator++()
  {
    if (!end_ && position_ < const_vec_.size() - 1)
      ++position_;
    else
      end_ = true;
    return *this;
  } // ... operator++()

  bool operator==(const ThisType& other)
  {
    return (end_ && other.end_) || ((!end_ && !other.end_) && (position_ == other.position_));
  }

  bool operator!=(const ThisType& other)
  {
    return !operator==(other);
  }

  const ScalarType& operator*() const
  {
    if (end_)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "This is the end!");
    return const_vec_[position_];
  }

private:
  const VectorType& const_vec_;

protected:
  size_t position_;
  bool end_;
}; // class VectorInputIterator

template <class Traits, class ScalarImp>
class VectorOutputIterator
  : public VectorInputIterator<Traits, ScalarImp>
  , public std::iterator<std::output_iterator_tag, typename Traits::ScalarType>
{
  typedef VectorInputIterator<Traits, ScalarImp> BaseType;
  typedef VectorOutputIterator<Traits, ScalarImp> ThisType;

public:
  typedef VectorInterface<Traits, ScalarImp> VectorType;
  typedef typename VectorType::ScalarType ScalarType;

private:
  static_assert(std::is_same<ScalarImp, ScalarType>::value, "");

public:
  explicit VectorOutputIterator(VectorType& vec, const bool end = false)
    : BaseType(vec, end)
    , vec_(vec)
  {}

  ScalarType& operator*()
  {
    if (this->end_)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "This is the end!");
    return vec_[this->position_];
  } // ... operator*()

private:
  VectorType& vec_;
}; // class VectorOutputIterator

} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_VECTOR_INTERFACE_INTERNAL_HH
