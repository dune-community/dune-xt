// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_COMMON_FLOAT_CMP_HH
#define DUNE_STUFF_COMMON_FLOAT_CMP_HH

#include <type_traits>
#include <complex>

#include <dune/common/float_cmp.hh>

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/vector.hh>
#include "float_cmp_style.hh"
#include "float_cmp_internal.hh"

namespace Dune {
namespace Stuff {
namespace Common {
namespace FloatCmp {

template <class T, Style style = Style::defaultStyle>
struct DefaultEpsilon
{
  typedef typename Dune::FloatCmp::EpsilonType<T>::Type Type;
  static Type value()
  {
    return Dune::FloatCmp::DefaultEpsilon<T, internal::ConvertStyle<style>::value>::value();
  }
};

//! since we treat complex like a vector its epsilon is (eps,eps) of its scalar type
template <class T, Style style>
struct DefaultEpsilon<std::complex<T>, style>
{
  typedef typename Dune::FloatCmp::EpsilonType<std::complex<T>>::Type Type;
  static Type value()
  {
    const auto val = Dune::FloatCmp::DefaultEpsilon<T, internal::ConvertStyle<style>::value>::value();
    return std::complex<T>(val, val);
  }
};

template <class T>
struct DefaultEpsilon<T, Style::numpy>
{
  typedef typename Dune::FloatCmp::EpsilonType<T>::Type Type;
  static Type value()
  {
    return Dune::FloatCmp::DefaultEpsilon<T, Dune::FloatCmp::relativeStrong>::value();
  }
};

//! necessary to avoid ambig. partial specialization
template <class T>
struct DefaultEpsilon<std::complex<T>, Style::numpy>
{
  typedef typename Dune::FloatCmp::EpsilonType<std::complex<T>>::Type Type;
  static Type value()
  {
    const auto val = Dune::FloatCmp::DefaultEpsilon<T, Dune::FloatCmp::relativeStrong>::value();
    return std::complex<T>(val, val);
  }
};

template <class V>
struct MT
{
  typedef typename VectorAbstraction<V>::S T;
  typedef typename Dune::FloatCmp::EpsilonType<typename MT<V>::T>::Type Eps;
};

#define DUNE_STUFF_COMMON_FLOAT_CMP_GENERATOR(id)                                                                      \
  template <Style style, class FirstType, class SecondType>                                                            \
  typename std::enable_if<internal::cmp_type_check<FirstType, SecondType, typename MT<FirstType>::T>::value,           \
                          bool>::type                                                                                  \
  id(const FirstType& first,                                                                                           \
     const SecondType& second,                                                                                         \
     const typename MT<FirstType>::Eps& rtol = DefaultEpsilon<typename MT<FirstType>::T, style>::value(),              \
     const typename MT<FirstType>::Eps& atol = DefaultEpsilon<typename MT<FirstType>::T, style>::value())              \
  {                                                                                                                    \
    return internal::Call<FirstType, SecondType, typename MT<FirstType>::Eps, style>::id(first, second, rtol, atol);   \
  }                                                                                                                    \
                                                                                                                       \
  template <class FirstType, class SecondType>                                                                         \
  typename std::enable_if<internal::cmp_type_check<FirstType, SecondType, typename MT<FirstType>::T>::value,           \
                          bool>::type                                                                                  \
  id(const FirstType& first,                                                                                           \
     const SecondType& second,                                                                                         \
     const typename MT<FirstType>::Eps& rtol =                                                                         \
         DefaultEpsilon<typename MT<FirstType>::T, Style::defaultStyle>::value(),                                      \
     const typename MT<FirstType>::Eps& atol =                                                                         \
         DefaultEpsilon<typename MT<FirstType>::T, Style::defaultStyle>::value())                                      \
  {                                                                                                                    \
    return id<Style::defaultStyle>(first, second, rtol, atol);                                                         \
  }

DUNE_STUFF_COMMON_FLOAT_CMP_GENERATOR(eq)
DUNE_STUFF_COMMON_FLOAT_CMP_GENERATOR(ne)
DUNE_STUFF_COMMON_FLOAT_CMP_GENERATOR(gt)
DUNE_STUFF_COMMON_FLOAT_CMP_GENERATOR(lt)
DUNE_STUFF_COMMON_FLOAT_CMP_GENERATOR(ge)
DUNE_STUFF_COMMON_FLOAT_CMP_GENERATOR(le)
#undef DUNE_STUFF_COMMON_FLOAT_CMP_GENERATOR

} // namespace FloatCmp
} // namespace Common
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_COMMON_FLOAT_CMP_HH
