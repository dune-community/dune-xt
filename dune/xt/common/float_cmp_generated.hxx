// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2016, 2018 - 2019)

template <Style style, class FirstType, class SecondType, class ToleranceType = typename MT<FirstType>::Eps>
typename std::enable_if<internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value, bool>::type
eq(const FirstType& first,
   const SecondType& second,
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value(),
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value())
{
  return internal::Call<FirstType, SecondType, typename Dune::FloatCmp::EpsilonType<ToleranceType>::Type, style>::eq(
      first, second, rtol, atol);
}

template <class FirstType, class SecondType, class ToleranceType = typename MT<FirstType>::Eps>
typename std::enable_if<internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value, bool>::type
eq(const FirstType& first,
   const SecondType& second,
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType>::value(),
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType>::value())
{
  return eq<Style::defaultStyle>(first, second, rtol, atol);
}

template <Style style, class FirstType, class SecondType, class ToleranceType = typename MT<FirstType>::Eps>
typename std::enable_if<internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value, bool>::type
ne(const FirstType& first,
   const SecondType& second,
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value(),
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value())
{
  return internal::Call<FirstType, SecondType, typename Dune::FloatCmp::EpsilonType<ToleranceType>::Type, style>::ne(
      first, second, rtol, atol);
}

template <class FirstType, class SecondType, class ToleranceType = typename MT<FirstType>::Eps>
typename std::enable_if<internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value, bool>::type
ne(const FirstType& first,
   const SecondType& second,
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType>::value(),
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType>::value())
{
  return ne<Style::defaultStyle>(first, second, rtol, atol);
}

template <Style style, class FirstType, class SecondType, class ToleranceType = typename MT<FirstType>::Eps>
typename std::enable_if<internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value, bool>::type
gt(const FirstType& first,
   const SecondType& second,
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value(),
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value())
{
  return internal::Call<FirstType, SecondType, typename Dune::FloatCmp::EpsilonType<ToleranceType>::Type, style>::gt(
      first, second, rtol, atol);
}

template <class FirstType, class SecondType, class ToleranceType = typename MT<FirstType>::Eps>
typename std::enable_if<internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value, bool>::type
gt(const FirstType& first,
   const SecondType& second,
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType>::value(),
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType>::value())
{
  return gt<Style::defaultStyle>(first, second, rtol, atol);
}

template <Style style, class FirstType, class SecondType, class ToleranceType = typename MT<FirstType>::Eps>
typename std::enable_if<internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value, bool>::type
lt(const FirstType& first,
   const SecondType& second,
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value(),
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value())
{
  return internal::Call<FirstType, SecondType, typename Dune::FloatCmp::EpsilonType<ToleranceType>::Type, style>::lt(
      first, second, rtol, atol);
}

template <class FirstType, class SecondType, class ToleranceType = typename MT<FirstType>::Eps>
typename std::enable_if<internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value, bool>::type
lt(const FirstType& first,
   const SecondType& second,
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType>::value(),
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType>::value())
{
  return lt<Style::defaultStyle>(first, second, rtol, atol);
}

template <Style style, class FirstType, class SecondType, class ToleranceType = typename MT<FirstType>::Eps>
typename std::enable_if<internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value, bool>::type
ge(const FirstType& first,
   const SecondType& second,
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value(),
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value())
{
  return internal::Call<FirstType, SecondType, typename Dune::FloatCmp::EpsilonType<ToleranceType>::Type, style>::ge(
      first, second, rtol, atol);
}

template <class FirstType, class SecondType, class ToleranceType = typename MT<FirstType>::Eps>
typename std::enable_if<internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value, bool>::type
ge(const FirstType& first,
   const SecondType& second,
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType>::value(),
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType>::value())
{
  return ge<Style::defaultStyle>(first, second, rtol, atol);
}

template <Style style, class FirstType, class SecondType, class ToleranceType = typename MT<FirstType>::Eps>
typename std::enable_if<internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value, bool>::type
le(const FirstType& first,
   const SecondType& second,
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value(),
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value())
{
  return internal::Call<FirstType, SecondType, typename Dune::FloatCmp::EpsilonType<ToleranceType>::Type, style>::le(
      first, second, rtol, atol);
}

template <class FirstType, class SecondType, class ToleranceType = typename MT<FirstType>::Eps>
typename std::enable_if<internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value, bool>::type
le(const FirstType& first,
   const SecondType& second,
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType>::value(),
   const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
       Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType>::value())
{
  return le<Style::defaultStyle>(first, second, rtol, atol);
}
