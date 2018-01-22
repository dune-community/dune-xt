// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_LA_EIGEN_CONTAINER_CONVERSION_HH
#define DUNE_XT_LA_EIGEN_CONTAINER_CONVERSION_HH

#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/la/type_traits.hh>
#include <dune/xt/la/container/matrix-interface.hh>

namespace Dune {
namespace XT {
namespace LA {


template <class RangeType, class SourceType>
typename std::enable_if<Common::is_matrix<SourceType>::value && Common::is_matrix<RangeType>::value, RangeType>::type
convert_to(const SourceType& source)
{
  return Common::convert_to<RangeType>(source);
}


template <class RangeType, class S>
typename std::enable_if<Common::is_matrix<RangeType>::value && !is_matrix<RangeType>::value, RangeType>::type
convert_to(const MatrixInterface<S>& source)
{
  const size_t rows = source.rows();
  const size_t cols = source.cols();
  auto ret = Common::create<RangeType>(rows, cols, 0);
  if (source.sparse) {
    const auto pattern = source.pattern();
    for (size_t ii = 0; ii < rows; ++ii)
      for (const auto& jj : pattern.inner(ii))
        Common::set_matrix_entry(ret,
                                 ii,
                                 jj,
#ifndef DXT_DISABLE_CHECKS
                                 Common::numeric_cast<typename Common::MatrixAbstraction<RangeType>::S>(
#endif
                                     source.get_entry(ii, jj)
#ifndef DXT_DISABLE_CHECKS
                                         )
#endif
                                     );
  } else {
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        Common::set_matrix_entry(ret,
                                 ii,
                                 jj,
#ifndef DXT_DISABLE_CHECKS
                                 Common::numeric_cast<typename Common::MatrixAbstraction<RangeType>::S>(
#endif
                                     source.get_entry(ii, jj)
#ifndef DXT_DISABLE_CHECKS
                                         )
#endif
                                     );
  }
  return ret;
} // ... convert_to(...)


template <class RangeType, class S>
typename std::enable_if<is_matrix<RangeType>::value, RangeType>::type convert_to(const MatrixInterface<S>& source)
{
  const size_t rows = source.rows();
  const size_t cols = source.cols();
  if (RangeType::sparse || source.sparse) {
    const auto pattern = source.pattern();
    RangeType ret(rows, cols, pattern);
    for (size_t ii = 0; ii < rows; ++ii)
      for (const auto& jj : pattern.inner(ii))
        ret.set_entry(ii,
                      jj,
#ifndef DXT_DISABLE_CHECKS
                      Common::numeric_cast<typename Common::MatrixAbstraction<RangeType>::S>(
#endif
                          source.get_entry(ii, jj)
#ifndef DXT_DISABLE_CHECKS
                              )
#endif
                          );
    return ret;
  } else {
    RangeType ret(rows, cols);
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        ret.set_entry(ii,
                      jj,
#ifndef DXT_DISABLE_CHECKS
                      Common::numeric_cast<typename Common::MatrixAbstraction<RangeType>::S>(
#endif
                          source.get_entry(ii, jj)
#ifndef DXT_DISABLE_CHECKS
                              )
#endif
                          );
    return ret;
  }
} // ... convert_to(...)


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_CONTAINER_CONVERSION_HH
