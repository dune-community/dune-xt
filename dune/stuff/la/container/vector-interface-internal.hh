// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_STUFF_LA_CONTAINER_VECTOR_INTERFACE_INTERNAL_HH
#define DUNE_STUFF_LA_CONTAINER_VECTOR_INTERFACE_INTERNAL_HH

#include <iterator>
#include <type_traits>

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/common/exceptions.hh>

namespace Dune {
namespace Stuff {
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

private:
  struct ConstHolder
  {
    explicit ConstHolder(const VectorType& vec)
      : element(vec)
    {
    }

    const VectorType& element;
  }; // struct ConstHolder

public:
  explicit VectorInputIterator(const VectorType& vec, const bool end = false)
    : const_holder_(std::make_shared<ConstHolder>(vec))
    , position_(0)
    , end_(end)
  {
  }

  ThisType& operator++()
  {
    if (!end_ && position_ < (const_holder_->element.size() - 1))
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
      DUNE_THROW(Exceptions::you_are_using_this_wrong, "This is the end!");
    return const_holder_->element[position_];
  }

private:
  std::shared_ptr<ConstHolder> const_holder_;

protected:
  size_t position_;
  bool end_;
}; // class VectorInputIterator

template <class Traits, class ScalarImp>
class VectorOutputIterator : public VectorInputIterator<Traits, ScalarImp>,
                             public std::iterator<std::output_iterator_tag, typename Traits::ScalarType>
{
  typedef VectorInputIterator<Traits, ScalarImp> BaseType;
  typedef VectorOutputIterator<Traits, ScalarImp> ThisType;

public:
  typedef VectorInterface<Traits, ScalarImp> VectorType;
  typedef typename VectorType::ScalarType ScalarType;

private:
  static_assert(std::is_same<ScalarImp, ScalarType>::value, "");

  struct Holder
  {
    explicit Holder(VectorType& vec)
      : element(vec)
    {
    }

    VectorType& element;
  }; // struct Holder

public:
  explicit VectorOutputIterator(VectorType& vec, const bool end = false)
    : BaseType(vec, end)
    , holder_(std::make_shared<Holder>(vec))
  {
  }

  ScalarType& operator*()
  {
    if (this->end_)
      DUNE_THROW(Exceptions::you_are_using_this_wrong, "This is the end!");
    return holder_->element[this->position_];
  } // ... operator*()

private:
  std::shared_ptr<Holder> holder_;
}; // class VectorOutputIterator

} // namespace internal
} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_VECTOR_INTERFACE_INTERNAL_HH
