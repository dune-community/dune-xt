// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2019)
//   René Fritze     (2013 - 2016, 2018 - 2019)
//   Tobias Leibner  (2018 - 2020)

#ifndef DUNE_XT_COMMON_PARALLEL_THREADSTORAGE_HH
#define DUNE_XT_COMMON_PARALLEL_THREADSTORAGE_HH

#if HAVE_TBB
#  include <tbb/enumerable_thread_specific.h>
#endif

#include <algorithm>
#include <list>
#include <mutex>
#include <numeric>
#include <system_error>
#include <type_traits>

#include <dune/xt/common/exceptions.hh>

// TODO: the following includes can be removed when UnsafePerThreadValue is removed
#include <deque>
#include <memory>
#include <boost/noncopyable.hpp>
#include <dune/xt/common/parallel/threadmanager.hh>

namespace Dune {
namespace XT {
namespace Common {
namespace internal {


#if HAVE_TBB

template <class ValueImp>
class EnumerableThreadSpecificWrapper
{
  // enumerable_thread_specific does not compile with ConstValueType as template param
  using BackendType = typename tbb::enumerable_thread_specific<std::remove_const_t<ValueImp>>;

public:
  using ValueType = ValueImp;
  using ConstValueType = std::add_const_t<ValueType>;
  using iterator = typename BackendType::iterator;
  using const_iterator = typename BackendType::const_iterator;

  template <class... InitTypes>
  explicit EnumerableThreadSpecificWrapper(InitTypes&&... ctor_args)
    : values_(std::forward<InitTypes>(ctor_args)...)
  {}

  ValueType& local()
  {
    return values_.local();
  }

  // tbb does not provide a const version of local (as elements may be inserted when a new thread accesses values_), so
  // values_ has to be mutable
  const ValueType& local() const
  {
    return values_.local();
  }

  typename BackendType::iterator begin()
  {
    return values_.begin();
  }

  typename BackendType::iterator end()
  {
    return values_.end();
  }

  typename BackendType::const_iterator begin() const
  {
    return values_.begin();
  }

  typename BackendType::const_iterator end() const
  {
    return values_.end();
  }

  template <class BinaryOperation>
  ValueType combine(BinaryOperation op) const
  {
    return values_.combine(op);
  }

private:
  mutable BackendType values_;
}; // class EnumerableThreadSpecificWrapper<ValueImp>

#else // HAVE_TBB

template <class ValueImp>
class EnumerableThreadSpecificWrapper
{
  using BackendType = std::unique_ptr<std::remove_const_t<ValueImp>>;

public:
  using ValueType = ValueImp;
  using ConstValueType = std::add_const_t<ValueType>;
  using iterator = ValueType*;
  using const_iterator = const ValueType*;

  //! Initialization by copy construction of ValueType
  explicit EnumerableThreadSpecificWrapper(ConstValueType& value)
    : values_(std::make_unique<std::remove_const_t<ValueType>>(value))
  {}

  //! Initialization by in-place construction ValueType with \param ctor_args
  template <class... InitTypes>
  explicit EnumerableThreadSpecificWrapper(InitTypes&&... ctor_args)
    : values_(std::make_unique<std::remove_const_t<ValueType>>(std::forward<InitTypes>(ctor_args)...))
  {}

  ValueType& local()
  {
    return *values_;
  }

  const ValueType& local() const
  {
    return *values_;
  }

  iterator begin()
  {
    return values_.get();
  }

  iterator end()
  {
    return values_.get() + 1;
  }

  const_iterator begin() const
  {
    return values_.get();
  }

  const_iterator end() const
  {
    return values_.get() + 1;
  }

  template <class BinaryOperation>
  ValueType combine(BinaryOperation /*op*/) const
  {
    return *values_;
  }

private:
  BackendType values_;
}; // class EnumerableThreadSpecificWrapper<ValueImp>

#endif // HAVE_TBB


} // namespace internal


/** Automatic Storage of non-static, N thread-local values
 **/
template <class ValueImp>
class PerThreadValue
{
  using ContainerType = internal::EnumerableThreadSpecificWrapper<ValueImp>;

public:
  using ValueType = typename ContainerType::ValueType;
  using ConstValueType = typename ContainerType::ConstValueType;

  //! Initialization by copy construction of ValueType
  explicit PerThreadValue(ConstValueType& value)
    : values_(value)
  {}

  //! Initialization by in-place construction ValueType with \param ctor_args
  template <class... InitTypes>
  explicit PerThreadValue(InitTypes&&... ctor_args)
    : values_(std::forward<InitTypes>(ctor_args)...)
  {}

  operator ValueType() const
  {
    return values_.local();
  }

  ValueType& operator*()
  {
    return values_.local();
  }

  ConstValueType& operator*() const
  {
    return values_.local();
  }

  ValueType* operator->()
  {
    return &values_.local();
  }

  ConstValueType* operator->() const
  {
    return &values_.local();
  }

  template <class BinaryOperation>
  ValueType accumulate(ValueType init, BinaryOperation op) const
  {
    return op(init, values_.combine(op));
  }

  ValueType sum() const
  {
    return accumulate(ValueType(0), std::plus<ValueType>());
  }

  typename ContainerType::iterator begin()
  {
    return values_.begin();
  }
  typename ContainerType::iterator end()
  {
    return values_.end();
  }

  typename ContainerType::const_iterator begin() const
  {
    return values_.begin();
  }
  typename ContainerType::const_iterator end() const
  {
    return values_.end();
  }

private:
  ContainerType values_;
}; // class PerThreadValue<ValueImp>


/**
 * Previous implementation of PerThreadValue. This implementation suffers from the fact that it is not possible (or
 * at least we did not find a way yet) to set a hard upper limit on the number of threads TBB uses. Setting max_threads
 * via tbb::task_scheduler_init apparently only sets a soft limit on the number of threads. In addition, even if TBB
 * uses only N threads at a time, it might be possible that a thread is destroyed and later in the program another
 * thread with a different id replaces it, which will then get a number greater than or equal to N in our
 * implementation (see ThreadManager::thread()). This occasionally leads to segfaults.
 * We keep this implementation around as it is currently used by TimingData (see dune/xt/common/timings.hh) and the
 * new implementation can't replace it in that context, as the new implementation based on
 * tbb::enumerable_thread_specific lazily initalizes the values in each thread.
 * \todo Either fix TimingData and remove this class or fix this class.
 **/
template <class ValueImp>
class UnsafePerThreadValue : public boost::noncopyable
{
public:
  typedef ValueImp ValueType;
  typedef typename std::conditional<std::is_const<ValueImp>::value, ValueImp, const ValueImp>::type ConstValueType;

private:
  typedef UnsafePerThreadValue<ValueImp> ThisType;
  typedef std::deque<std::unique_ptr<ValueType>> ContainerType;

public:
  //! Initialization by copy construction of ValueType
  explicit UnsafePerThreadValue(ConstValueType& value)
    : values_(threadManager().max_threads())
  {
    std::generate(values_.begin(), values_.end(), [=]() { return std::make_unique<ValueType>(value); });
  }

  //! Initialization by in-place construction ValueType with \param ctor_args
  template <class... InitTypes>
  explicit UnsafePerThreadValue(InitTypes&&... ctor_args)
    : values_(threadManager().max_threads())
  {
    for (auto&& val : values_)
      val = std::make_unique<ValueType>(ctor_args...);
  }

  ThisType& operator=(ConstValueType&& value)
  {
    std::generate(values_.begin(), values_.end(), [=]() { return std::make_unique<ValueType>(value); });
    return *this;
  }

  operator ValueType() const
  {
    return this->operator*();
  }

  ValueType& operator*()
  {
    return *values_[threadManager().thread()];
  }

  ConstValueType& operator*() const
  {
    return *values_[threadManager().thread()];
  }

  ValueType* operator->()
  {
    return values_[threadManager().thread()].get();
  }

  ConstValueType* operator->() const
  {
    return values_[threadManager().thread()].get();
  }

  auto& get_pointer()
  {
    return values_[threadManager().thread()];
  }

  template <class BinaryOperation>
  ValueType accumulate(ValueType init, BinaryOperation op) const
  {
    typedef const typename ContainerType::value_type ptr;
    auto l = [&](ConstValueType& a, ptr& b) { return op(a, *b); };
    return std::accumulate(values_.begin(), values_.end(), init, l);
  }

  ValueType sum() const
  {
    return accumulate(ValueType(0), std::plus<ValueType>());
  }

  typename ContainerType::iterator begin()
  {
    return values_.begin();
  }
  typename ContainerType::iterator end()
  {
    return values_.end();
  }

  typename ContainerType::const_iterator begin() const
  {
    return values_.begin();
  }
  typename ContainerType::const_iterator end() const
  {
    return values_.end();
  }

private:
  ContainerType values_;
}; // class UnsafePerThreadValue<...>


/**
 * \brief Propagates results from copies in each thread (usually stemming from PerThreadValue) to base functor on
 *calling finalize.
 * \param imp_ Pointer to implementation
 * \param base_ Pointer to initial base functor
 **/
template <class Imp, typename Result, class Reduction = std::plus<Result>>
class ThreadResultPropagator
{
public:
  ThreadResultPropagator(Imp* imp)
    : imp_(imp)
    , base_(imp)
  {}

  ThreadResultPropagator(const ThreadResultPropagator& other)
    : imp_(other.imp_)
    , base_(other.base_)
    , mutex_()
  {}

  Imp* copy_imp()
  {
    auto* cpy = new Imp(*imp_);
    cpy->base_ = imp_->base_;
    cpy->imp_ = cpy;
    return cpy;
  }

  void finalize_imp()
  {
    Reduction reduce;
    base_->set_result(reduce(base_->result(), imp_->result()));
  }

protected:
  Imp* imp_;
  Imp* base_;
  std::mutex mutex_;
};


template <typename T>
struct concatenate_container
{
  T operator()(const T& a, const T& b)
  {
    T result = a;
    result.insert(b.begin(), b.end());
    return result;
  }
}; // struct concatenate_container


} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_COMMON_PARALLEL_THREADSTORAGE_HH
