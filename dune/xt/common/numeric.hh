// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019 - 2020)

#ifndef DUNE_XT_COMMON_NUMERIC_HH
#define DUNE_XT_COMMON_NUMERIC_HH

// C++17 parallel TS features
#if defined(__has_include) && defined(__cpp_lib_execution) && defined(__cpp_lib_parallel_algorithm)
#  define CPP17_PARALLELISM_TS_SUPPORTED                                                                               \
    __has_include(<execution>) && __cpp_lib_execution >= 201603 && __cpp_lib_parallel_algorithm >= 201603
#else
#  define CPP17_PARALLELISM_TS_SUPPORTED 0
#endif

#include <numeric>

namespace Dune {
namespace XT {
namespace Common {


template <class... Args>
decltype(auto) reduce(Args&&... args)
{
#if CPP17_PARALLELISM_TS_SUPPORTED
  return std::reduce(std::forward<Args>(args)...);
#else
  return std::accumulate(std::forward<Args>(args)...);
#endif
}

template <class... Args>
decltype(auto) transform_reduce(Args&&... args)
{
#if CPP17_PARALLELISM_TS_SUPPORTED
  return std::transform_reduce(std::forward<Args>(args)...);
#else
  return std::inner_product(std::forward<Args>(args)...);
#endif
}


} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_COMMON_NUMERIC_HH
