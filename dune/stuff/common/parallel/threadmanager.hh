// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_COMMON_THREADMANAGER_HH
#define DUNE_STUFF_COMMON_THREADMANAGER_HH

#include <vector>
#include <algorithm>
#include <type_traits>

namespace Dune {
namespace Stuff {

/** abstractions of threading functionality
 *  currently forwads to dune-fem if possible, falls back to single-thread dummy imp
 **/
struct ThreadManager
{
  //! return maximal number of threads possbile in the current run
  static unsigned int max_threads();

  //! return number of current threads
  static unsigned int current_threads();

  //! return thread number
  static unsigned int thread();

  //! set maximal number of threads available during run
  static void set_max_threads(const unsigned int count);

private:
  //  std::unique_ptr<tbb::task_scheduler_init> tbb_init_;
};

template <class ValueImp>
class PerThreadValue
{
public:
  typedef ValueImp ValueType;

private:
  typedef PerThreadValue<ValueImp> ThisType;
  typedef std::vector<ValueType> ContainerType;

public:
  //! Initialization by copy construction of ValueType
  explicit PerThreadValue(const ValueType& value)
    : values_(ThreadManager::max_threads(), value)
  {
  }

  //! Initialization by in-place construction ValueType with \param ctor_args
  template <class... InitTypes>
  PerThreadValue(InitTypes&&... ctor_args)
    : values_(ThreadManager::max_threads(), ValueType(std::forward<InitTypes>(ctor_args)...))
  {
  }

  ValueType& operator*()
  {
    return values_[ThreadManager::thread()];
  }

  const ValueType& operator*() const
  {
    return values_[ThreadManager::thread()];
  }

  ValueType* operator->()
  {
    return &values_[ThreadManager::thread()];
  }

  const ValueType* operator->() const
  {
    return &values_[ThreadManager::thread()];
  }

  typename ContainerType::size_type size() const
  {
    return values_.size();
  }

  template <class BinaryOperation>
  ValueType accumulate(ValueType init, BinaryOperation op) const
  {
    return std::accumulate(values_.begin(), values_.end(), init);
  }

  ValueType sum() const
  {
    return accumulate(ValueType(0), std::plus<ValueType>());
  }

private:
  ContainerType values_;
};
}
}

#endif // DUNE_STUFF_COMMON_THREADMANAGER_HH