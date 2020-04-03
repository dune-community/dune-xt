// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_COMMON_DENSEVECTOR_HH
#define DUNE_XT_COMMON_DENSEVECTOR_HH

#include <dune/common/densevector.hh>

#include <dune/xt/common/type_traits.hh>

namespace Dune {


template <class S, class V>
typename std::enable_if<Dune::XT::Common::is_arithmetic<S>::value || Dune::XT::Common::is_complex<S>::value,
                        typename DenseVector<V>::derived_type>::type
operator*(const S& alpha, const DenseVector<V>& vec)
{
  typename DenseVector<V>::derived_type ret = vec;
  ret *= alpha;
  return ret;
}


template <class V, class S>
typename std::enable_if<Dune::XT::Common::is_arithmetic<S>::value || Dune::XT::Common::is_complex<S>::value,
                        typename DenseVector<V>::derived_type>::type
operator*(const DenseVector<V>& vec, const S& alpha)
{
  typename DenseVector<V>::derived_type ret = vec;
  ret *= alpha;
  return ret;
}


template <class V, class S>
typename std::enable_if<Dune::XT::Common::is_arithmetic<S>::value || Dune::XT::Common::is_complex<S>::value,
                        typename DenseVector<V>::derived_type>::type
operator/(const DenseVector<V>& vec, const S& alpha)
{
  typename DenseVector<V>::derived_type ret = vec;
  ret /= alpha;
  return ret;
}


} // namespace Dune

#endif // DUNE_XT_COMMON_DENSEVECTOR_HH
