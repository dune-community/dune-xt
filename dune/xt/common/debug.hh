// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014, 2016 - 2017)
//   René Fritze     (2010 - 2019)
//   Tobias Leibner  (2018 - 2020)

#ifndef DUNE_XT_COMMON_DEBUG_HH
#define DUNE_XT_COMMON_DEBUG_HH

#include <cstring>

#include <dune/xt/common/unused.hh>
#include <dune/xt/common/exceptions.hh>


#define SEGFAULT                                                                                                       \
  {                                                                                                                    \
    int* J = 0;                                                                                                        \
    *J = 9;                                                                                                            \
  }


inline char* charcopy(const char* s)
{
  size_t l = strlen(s) + 1;
  char* t = new char[l];
  for (size_t i = 0; i < l; i++) {
    t[i] = s[i];
  }
  return t;
} // ... charcopy(...)


//! try to ensure var is not optimized out
#define DXTC_DEBUG_AUTO(name) volatile auto DXTC_UNUSED(name)

#ifndef NDEBUG
#  define DXT_ASSERT(condition)                                                                                        \
    DUNE_THROW_IF(!(condition),                                                                                        \
                  Dune::XT::Common::Exceptions::debug_assertion,                                                       \
                  __PRETTY_FUNCTION__ << "\nAssertion failed: \n"                                                      \
                                      << #condition)
#else
#  define DXT_ASSERT(condition)
#endif

#endif // DUNE_XT_COMMON_DEBUG_HH
