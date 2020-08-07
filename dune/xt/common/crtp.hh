// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2014, 2016 - 2018)
//   René Fritze     (2014 - 2016, 2018 - 2019)
//   Tobias Leibner  (2018, 2020)

#ifndef DUNE_XT_COMMON_CRTP_HH
#define DUNE_XT_COMMON_CRTP_HH

#include <atomic>
#include <mutex>

#include <dune/xt/common/deprecated.hh>

#include "exceptions.hh"

#ifdef CHECK_CRTP
#  undef CHECK_CRTP
#endif
#ifdef CHECK_AND_CALL_CRTP
#  undef CHECK_AND_CALL_CRTP
#endif

#ifdef NDEBUG
#  define CHECK_CRTP(dummy)
#else

/**
 * This macro is essentially a thread safe variant of the CHECK_INTERFACE_IMPLEMENTATION macro from dune-common.
 */
#  define CHECK_CRTP(__interface_method_to_call__)                                                                     \
    {                                                                                                                  \
      std::lock_guard<std::recursive_mutex> crtp_mutex_guard(this->crtp_mutex_);                                       \
      static std::atomic<bool> call(false);                                                                            \
      if (call)                                                                                                        \
        DUNE_THROW(Dune::XT::Common::Exceptions::CRTP_check_failed,                                                    \
                   "The derived class does not implement the required method!");                                       \
      call = true;                                                                                                     \
      try {                                                                                                            \
        (__interface_method_to_call__);                                                                                \
        call = false;                                                                                                  \
      } catch (...) {                                                                                                  \
        call = false;                                                                                                  \
        throw;                                                                                                         \
      }                                                                                                                \
    }
// CHECK_CRTP
#endif // NDEBUG


#ifdef NDEBUG
#  define CHECK_AND_CALL_CRTP(__interface_method_to_call__) (__interface_method_to_call__)
#else
/**
 * This macro is essentially a slightly modified copy of the CHECK_AND_CALL_INTERFACE_IMPLEMENTATION macro.
 */
#  define CHECK_AND_CALL_CRTP(__interface_method_to_call__) CHECK_CRTP(__interface_method_to_call__)
#endif


namespace Dune::XT {
namespace Common {


template <class InterfaceType, class TraitsImp>
class CRTPInterface
{
public:
  using Traits = TraitsImp;
  using derived_type = typename Traits::derived_type;

  DXT_DEPRECATED_MSG("Use that.as_imp() instead (10.04.2018)!")
  static inline derived_type& as_imp(InterfaceType& that)
  {
    return static_cast<derived_type&>(that);
  }

  DXT_DEPRECATED_MSG("Use that.as_imp() instead (10.04.2018)!")
  static inline const derived_type& as_imp(const InterfaceType& that)
  {
    return static_cast<const derived_type&>(that);
  }

  inline derived_type& as_imp()
  {
    return static_cast<derived_type&>(*this);
  }

  inline const derived_type& as_imp() const
  {
    return static_cast<const derived_type&>(*this);
  }

protected:
  // nicely avoid warning about non-virtual dtor when derived classes have vfunc
  ~CRTPInterface() = default;

#ifndef NDEBUG
  // needs to be static so diff instances don't clash in function local crtp check
  static std::recursive_mutex crtp_mutex_;
#endif
}; // CRTPInterface


#ifndef NDEBUG
template <class I, class T>
std::recursive_mutex CRTPInterface<I, T>::crtp_mutex_;
#endif


} // namespace Common


template <class I, class T>
using CRTPInterface DXT_DEPRECATED_MSG("CRTPInterface was moved into Dune::XT:Common namespace (2018/3/16)") =
    Common::CRTPInterface<I, T>;


} // namespace Dune::XT

#endif // DUNE_XT_COMMON_CRTP_HH
