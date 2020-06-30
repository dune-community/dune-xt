// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2017 - 2019)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_COMMON_VECTOR_STATISTICS_HH
#define DUNE_XT_COMMON_VECTOR_STATISTICS_HH

#include <algorithm>

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/vector.hh>

#include <boost/accumulators/statistics/variance.hpp>

namespace Dune {
namespace XT {
namespace Common {

template <class VectorType>
typename std::enable_if<is_vector<VectorType>::value, typename VectorAbstraction<VectorType>::S>::type
standard_deviation(const VectorType& vector)
{
  using namespace boost::accumulators;
  typedef typename VectorAbstraction<VectorType>::S FieldType;
  accumulator_set<FieldType, stats<tag::lazy_variance>> acc;
  std::for_each(vector.begin(), vector.end(), [&](FieldType value) { acc(value); });
  return std::sqrt(variance(acc));
}


} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_COMMON_VECTOR_STATISTICS_HH
