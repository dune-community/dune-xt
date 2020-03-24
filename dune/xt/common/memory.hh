// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014 - 2019)
//   René Fritze     (2011 - 2012, 2014 - 2019)
//   Tobias Leibner  (2017, 2019)

#ifndef DUNE_XT_COMMON_MEMORY_HH
#define DUNE_XT_COMMON_MEMORY_HH

#include <cassert>
#include <memory>
#include <boost/noncopyable.hpp>

#include <dune/xt/common/debug.hh>

namespace Dune {
namespace XT {
namespace Common {


//! just like boost::noncopyable, but for move assign/ctor
struct nonmoveable
{
  constexpr nonmoveable() = default;
  ~nonmoveable() = default;
  nonmoveable& operator=(nonmoveable&& source) = delete;
  explicit nonmoveable(nonmoveable&& source) = delete;
};


namespace internal {


template <class T>
class ConstAccessInterface
{
public:
  virtual ~ConstAccessInterface() = default;

  ConstAccessInterface<T>& operator=(const ConstAccessInterface<T>& other) = delete;
  ConstAccessInterface<T>& operator=(ConstAccessInterface<T>&& source) = delete;

  virtual const T& access() const = 0;

  virtual ConstAccessInterface<T>* copy() const = 0;
}; // class ConstAccessInterface


template <class T>
class ConstAccessByReference : public ConstAccessInterface<T>
{
public:
  explicit ConstAccessByReference(const T& tt)
    : tt_(tt)
  {}

  const T& access() const override final
  {
    return tt_;
  }

  ConstAccessInterface<T>* copy() const override final
  {
    return new ConstAccessByReference<T>(tt_);
  }

private:
  const T& tt_;
}; // class ConstAccessByReference


template <class T>
class ConstAccessByPointer : public ConstAccessInterface<T>
{
public:
  /**
   * \attention This ctor transfers ownership to ConstAccessByPointer, do not delete tt manually!
   */
  explicit ConstAccessByPointer(const T*&& tt)
    : tt_(tt)
  {}

  explicit ConstAccessByPointer(std::unique_ptr<const T>&& tt)
    : tt_(tt)
  {}

  explicit ConstAccessByPointer(std::shared_ptr<const T> tt)
    : tt_(tt)
  {}

  const T& access() const override final
  {
    return *tt_;
  }

  ConstAccessInterface<T>* copy() const override final
  {
    return new ConstAccessByPointer<T>(tt_);
  }

private:
  std::shared_ptr<const T> tt_;
}; // class ConstAccessByPointer


template <class T>
class AccessInterface
{
public:
  virtual ~AccessInterface() = default;

  AccessInterface<T>& operator=(const AccessInterface<T>& other) = delete;
  AccessInterface<T>& operator=(AccessInterface<T>&& source) = delete;

  virtual T& access() = 0;

  virtual const T& access() const = 0;

  virtual AccessInterface<T>* copy() = 0;
}; // class AccessInterface


template <class T>
class AccessByReference : public AccessInterface<T>
{
public:
  explicit AccessByReference(T& tt)
    : tt_(tt)
  {}

  T& access() override final
  {
    return tt_;
  }

  const T& access() const override final
  {
    return tt_;
  }

  AccessInterface<T>* copy() override final
  {
    return new AccessByReference<T>(tt_);
  }

private:
  T& tt_;
}; // class AccessByReference


template <class T>
class AccessByPointer : public AccessInterface<T>
{
public:
  /**
   * \attention This ctor transfers ownership to AccessByPointer, do not delete tt manually!
   */
  explicit AccessByPointer(T*&& tt)
    : tt_(tt)
  {}

  explicit AccessByPointer(std::unique_ptr<T>&& tt)
    : tt_(tt)
  {}

  explicit AccessByPointer(std::shared_ptr<T> tt)
    : tt_(tt)
  {}

  T& access() override final
  {
    return *tt_;
  }

  const T& access() const override final
  {
    return *tt_;
  }

  AccessInterface<T>* copy() override final
  {
    return new AccessByPointer<T>(tt_);
  }

private:
  std::shared_ptr<T> tt_;
}; // class AccessByPointer


} // namespace internal


/**
 * \brief Provides generic (const) access to objects of different origins.
 *
 * Cosider the following base class which always requires a vector:
\code
struct Base
{
  Base(const std::vector<double>& vec);
};
\endcode
 * Consider further a derived class which should be constructible with a given vector as well as without:
\code
struct Derived : public Base
{
  Derived(const std::vector<double>& vec);
  Derived();
};
\endcode
 * In the latter case, a vector should be created automatically, which is problematic due to the requirements of Base.
 * The solution is to first derive from ConstStorageProvider or StorageProvider, which handles the management of the
 * vector:
\code
struct Derived
 : private ConstStorageProvider<std::vector<double>>
 , public Base
{
  typedef ConstStorageProvider<std::vector<double>> VectorProvider;

  Derived(const std::vector<double>& vec)
    : VectorProvider(vec)
    , Base(VectorProvider::access())
  {}

  Derived()
    : VectorProvider(new std::vector<double>())
    , Base(VectorProvider::access())
};
\endcode
 * For the latter to work, ConstStorageProvider (as well as StorageProvider) needs to take ownership of the provided raw
 * pointer.
 * \attention ConstStorageProvider (as well as StorageProvider) takes ownership of the provided raw pointer. Thus, the
 *            following code is supposed to fail:
\code
const T* tt = new T();
{
  ConstStorageProvider<T> provider(tt);
}
const T& derefed_tt = *tt;
\endcode
 *            Do the following instead:
\code
const T* tt = new T();
{
  ConstStorageProvider<T> provider(*tt);
}
const T& derefed_tt = *tt;
 */
template <class T>
class ConstStorageProvider
{
public:
  explicit ConstStorageProvider()
    : storage_(nullptr)
  {}

  explicit ConstStorageProvider(const T& tt)
    : storage_(new internal::ConstAccessByReference<T>(tt))
  {}

  /**
   * \attention This ctor transfers ownership to ConstStorageProvider, do not delete tt manually!
   */
  explicit ConstStorageProvider(const T*&& tt)
    : storage_(new internal::ConstAccessByPointer<T>(std::move(tt)))
  {}

  /**
   * \attention This ctor transfers ownership to ConstStorageProvider, do not delete tt manually!
   */
  explicit ConstStorageProvider(T*&& tt)
    : storage_(new internal::ConstAccessByPointer<T>(std::move(tt)))
  {}

  explicit ConstStorageProvider(std::shared_ptr<const T> tt)
    : storage_(new internal::ConstAccessByPointer<T>(tt))
  {}

  explicit ConstStorageProvider(std::unique_ptr<const T>&& tt)
    : storage_(tt)
  {}

  ConstStorageProvider(const ConstStorageProvider<T>& other)
    : storage_(other.storage_ ? other.storage_->copy() : nullptr)
  {}

  ConstStorageProvider(ConstStorageProvider<T>&& source) = default;

  ConstStorageProvider<T>& operator=(const ConstStorageProvider<T>& other)
  {
    if (&other != this)
      storage_ = std::unique_ptr<internal::ConstAccessInterface<T>>(other.storage_->copy());
    return *this;
  }

  ConstStorageProvider<T>& operator=(ConstStorageProvider<T>&& source) = default;

  bool valid() const
  {
    return !(storage_ == nullptr);
  }

  const T& access() const
  {
    DXT_ASSERT(valid());
    return storage_->access();
  }

  template <typename... Args>
  static ConstStorageProvider<T> make(Args&&... args)
  {
    return ConstStorageProvider<T>(new T(std::forward<Args>(args)...));
  }

private:
  std::unique_ptr<const internal::ConstAccessInterface<T>> storage_;
}; // class ConstStorageProvider


template <typename T, typename... Args>
ConstStorageProvider<T> make_const_storage(Args&&... args)
{
  return ConstStorageProvider<T>::make(std::forward<Args>(args)...);
}


/**
 * \brief Provides generic access to objects of different origins.
 * \sa ConstStorageProvider
 */
template <class T>
class StorageProvider
{
public:
  explicit StorageProvider()
    : storage_(nullptr)
  {}

  explicit StorageProvider(T& tt)
    : storage_(new internal::AccessByReference<T>(tt))
  {}

  /**
   * \attention This ctor transfers ownership to StorageProvider, do not delete it manually!
   */
  explicit StorageProvider(T*&& tt)
    : storage_(new internal::AccessByPointer<T>(std::move(tt)))
  {}

  explicit StorageProvider(std::shared_ptr<T> tt)
    : storage_(new internal::AccessByPointer<T>(tt))
  {}

  StorageProvider(const StorageProvider<T>& other) = delete;

  StorageProvider(StorageProvider<T>& other)
    : storage_(other.storage_ ? other.storage_->copy() : nullptr)
  {}

  StorageProvider(StorageProvider<T>&& source) = default;

  StorageProvider<T>& operator=(const StorageProvider<T>& other) = delete;

  StorageProvider<T>& operator=(StorageProvider<T>& other)
  {
    if (&other != this)
      storage_ = std::unique_ptr<internal::AccessInterface<T>>(other.storage_->copy());
    return *this;
  }

  StorageProvider<T>& operator=(StorageProvider<T>&& source) = default;

  bool valid() const
  {
    return !(storage_ == nullptr);
  }

  T& access()
  {
    DXT_ASSERT(valid());
    return storage_->access();
  }

  const T& access() const
  {
    DXT_ASSERT(valid());
    return storage_->access();
  }

  template <typename... Args>
  static StorageProvider<T> make(Args&&... args)
  {
    return StorageProvider<T>(new T(std::forward<Args>(args)...));
  }

private:
  std::unique_ptr<internal::AccessInterface<T>> storage_;
}; // class StorageProvider


template <typename T, typename... Args>
StorageProvider<T> make_storage(Args&&... args)
{
  return StorageProvider<T>::make(std::forward<Args>(args)...);
}


/**
 * \brief Provides generic (const) access to objects of different origins (if applicable).
 *
 * \sa ConstStorageProvider
 */
template <class T>
class ConstSharedStorageProvider
{
public:
  /**
   * \attention This ctor transfers ownership to ConstSharedStorageProvider, do not delete tt manually!
   */
  explicit ConstSharedStorageProvider(const T*&& tt)
    : storage_(std::move(tt))
  {}

  /**
   * \attention This ctor transfers ownership to ConstSharedStorageProvider, do not delete tt manually!
   */
  explicit ConstSharedStorageProvider(T*&& tt)
    : storage_(std::move(tt))
  {}

  explicit ConstSharedStorageProvider(std::shared_ptr<const T> tt)
    : storage_(tt)
  {}

  explicit ConstSharedStorageProvider(std::shared_ptr<T> tt)
    : storage_(tt)
  {}

  explicit ConstSharedStorageProvider(std::unique_ptr<const T>&& tt)
    : storage_(tt.release())
  {}

  explicit ConstSharedStorageProvider(std::unique_ptr<T>&& tt)
    : storage_(tt.release())
  {}

  ConstSharedStorageProvider(const ConstSharedStorageProvider<T>& other) = default;
  ConstSharedStorageProvider(ConstSharedStorageProvider<T>&& source) = default;

  ConstSharedStorageProvider<T>& operator=(const ConstSharedStorageProvider<T>& other) = default;
  ConstSharedStorageProvider<T>& operator=(ConstSharedStorageProvider<T>&& source) = default;

  std::shared_ptr<const T> access() const
  {
    DXT_ASSERT(storage_);
    return storage_;
  }

  template <typename... Args>
  static ConstSharedStorageProvider<T> make(Args&&... args)
  {
    return ConstSharedStorageProvider<T>(new T(std::forward<Args>(args)...));
  }

private:
  std::shared_ptr<const T> storage_;
}; // class ConstSharedStorageProvider


template <typename T, typename... Args>
ConstSharedStorageProvider<T> make_const_shared_storage(Args&&... args)
{
  return ConstSharedStorageProvider<T>::make(std::forward<Args>(args)...);
}


//! dumps kernel stats into a file
void mem_usage(std::string filename);

void mem_usage();


} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_COMMON_MEMORY_HH
