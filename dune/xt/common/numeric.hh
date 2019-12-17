// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#ifndef DUNE_XT_COMMON_NUMERIC_HH
#define DUNE_XT_COMMON_NUMERIC_HH

// C++17 parallel TS features
#if defined __has_include
#  define CPP17_PARALLELISM_TS_SUPPORTED __has_include(<execution>)
#else
#  define CPP17_PARALLELISM_TS_SUPPORTED 0
#endif

#include <numeric>

namespace Dune {
namespace XT {
namespace Common {


template <class InputIt, class T>
T reduce(InputIt first, InputIt last, T init)
{
#if CPP17_PARALLELISM_TS_SUPPORTED
  return std::reduce(first, last, init);
#else
  return std::accumulate(first, last, init);
#endif
}

template <class InputIt, class T, class BinaryOp>
T reduce(InputIt first, InputIt last, T init, BinaryOp binary_op)
{
#if CPP17_PARALLELISM_TS_SUPPORTED
  return std::reduce(first, last, init, binary_op);
#else
  return std::accumulate(first, last, init, binary_op);
#endif
}

template <class InputIt1, class InputIt2, class T>
T transform_reduce(InputIt1 first1, InputIt1 last1, InputIt2 first2, T init)
{
#if CPP17_PARALLELISM_TS_SUPPORTED
  return std::transform_reduce(first1, last1, first2, init);
#else
  return std::inner_product(first1, last1, first2, init);
#endif
}

template <class InputIt1, class InputIt2, class T, class BinaryOp1, class BinaryOp2>
T transform_reduce(InputIt1 first1, InputIt1 last1, InputIt2 first2, T init, BinaryOp1 binary_op1, BinaryOp2 binary_op2)
{
#if CPP17_PARALLELISM_TS_SUPPORTED
  return std::transform_reduce(first1, last1, first2, init, binary_op1, binary_op2);
#else
  return std::inner_product(first1, last1, first2, init, binary_op1, binary_op2);
#endif
}


} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_COMMON_NUMERIC_HH
