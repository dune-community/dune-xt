// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   René Fritze    (2020)
//   Tobias Leibner (2019 - 2020)

#ifndef DUNE_XT_COMMON_NUMERIC_HH
#define DUNE_XT_COMMON_NUMERIC_HH

#include <numeric>

#if defined(__cpp_lib_parallel_algorithm) && __cpp_lib_parallel_algorithm >= 201603
#  define CPP17_PARALLELISM_TS_SUPPORTED 1
#else
#  define CPP17_PARALLELISM_TS_SUPPORTED 0
#endif

namespace Dune::XT::Common {


// Uses std::reduce if available, and falls back to std::accumulate on older compilers.
// The std::reduce versions with an execution policy as first argument are not supported.
template <class... Args>
decltype(auto) reduce(Args&&... args)
{
#if CPP17_PARALLELISM_TS_SUPPORTED
  return std::reduce(std::forward<Args>(args)...);
#else
  return std::accumulate(std::forward<Args>(args)...);
#endif
}

// Uses std::transform_reduce if available, and falls back to std::inner_product on older compilers.
// The std::transform_reduce versions with an execution policy as first argument are not supported.
template <class... Args>
decltype(auto) transform_reduce(Args&&... args)
{
#if CPP17_PARALLELISM_TS_SUPPORTED
  return std::transform_reduce(std::forward<Args>(args)...);
#else
  return std::inner_product(std::forward<Args>(args)...);
#endif
}


} // namespace Dune::XT::Common

#endif // DUNE_XT_COMMON_NUMERIC_HH
