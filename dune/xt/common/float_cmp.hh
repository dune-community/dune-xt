// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2013 - 2018)
//   René Fritze      (2013 - 2016, 2018 - 2019)
//   Tobias Leibner   (2014 - 2015, 2018, 2020)

#ifndef DUNE_XT_COMMON_FLOAT_CMP_HH
#define DUNE_XT_COMMON_FLOAT_CMP_HH

#include <type_traits>
#include <complex>

#include <dune/xt/common/math.hh> // <- This include needs to be before the one from dune-common, otherwise
#include <dune/common/float_cmp.hh> //  std::abs(long unsinged int) is indefined in dune-common!

#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/type_traits.hh>

#include "float_cmp_style.hh"
#include "float_cmp_internal.hh"

namespace Dune {
namespace XT {
namespace Common {
namespace FloatCmp {


template <class T, Style style = Style::defaultStyle, bool = std::is_integral<T>::value>
struct DefaultEpsilon
{
  typedef typename Dune::FloatCmp::EpsilonType<T>::Type Type;
  static Type value()
  {
    return Dune::FloatCmp::DefaultEpsilon<T, internal::ConvertStyle<style>::value>::value();
  }
};

template <class T, Style style>
struct DefaultEpsilon<std::complex<T>, style, false>
{
  typedef typename Dune::FloatCmp::EpsilonType<std::complex<T>>::Type Type;
  static Type value()
  {
    const auto val = Dune::FloatCmp::DefaultEpsilon<T, internal::ConvertStyle<style>::value>::value();
    return std::complex<T>(val, val);
  }
};

template <class T>
struct DefaultEpsilon<T, Style::numpy, false>
{
  typedef typename Dune::FloatCmp::EpsilonType<T>::Type Type;
  static Type value()
  {
    return Dune::FloatCmp::DefaultEpsilon<T, Dune::FloatCmp::relativeStrong>::value();
  }
};

//! necessary to avoid ambig. partial specialization
template <class T>
struct DefaultEpsilon<std::complex<T>, Style::numpy, false>
{
  typedef typename Dune::FloatCmp::EpsilonType<std::complex<T>>::Type Type;
  static Type value()
  {
    const auto val = Dune::FloatCmp::DefaultEpsilon<T, Dune::FloatCmp::relativeStrong>::value();
    return std::complex<T>(val, val);
  }
};

//! these integral specialisations are only necessary because dune-common's DefaultEpsilon
//! tries to instantiate invalid code for integers
template <class T, Style style>
struct DefaultEpsilon<T, style, true>
{
  typedef T Type;
  static T value()
  {
    return std::numeric_limits<T>::epsilon();
  }
};

template <class T>
struct DefaultEpsilon<T, Style::numpy, true>
{
  typedef T Type;
  static T value()
  {
    return std::numeric_limits<T>::epsilon();
  }
};

template <class V>
struct MT
{
  typedef typename std::conditional<Dune::XT::Common::is_matrix<V>::value,
                                    typename Dune::XT::Common::MatrixAbstraction<V>::S,
                                    typename Dune::XT::Common::VectorAbstraction<V>::S>::type T;
  typedef typename Dune::FloatCmp::EpsilonType<typename MT<V>::T>::Type Eps;
};

#include "float_cmp_generated.hxx"

} // namespace FloatCmp
} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_COMMON_FLOAT_CMP_HH
