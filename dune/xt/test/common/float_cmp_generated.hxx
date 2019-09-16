// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2016, 2018)

template <Dune::XT::Common::FloatCmp::Style style,
          class FirstType,
          class SecondType,
          class ToleranceType = typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps>
typename std::enable_if<
    Dune::XT::Common::FloatCmp::internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value,
    void>::type
DXTC_EXPECT_FLOAT_EQ(const FirstType& expected,
                     const SecondType& actual,
                     const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
                         Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value(),
                     const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
                         Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value())
{
  const auto eq = Dune::XT::Common::FloatCmp::internal::
      Call<FirstType, SecondType, typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps, style>::eq;
  EXPECT_PRED4(eq, expected, actual, rtol, atol);
}

template <class FirstType,
          class SecondType,
          class ToleranceType = typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps>
typename std::enable_if<
    Dune::XT::Common::FloatCmp::internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value,
    void>::type
DXTC_EXPECT_FLOAT_EQ(
    const FirstType& first,
    const SecondType& second,
    const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
        Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType,
                                                   Dune::XT::Common::FloatCmp::Style::defaultStyle>::value(),
    const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
        Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType,
                                                   Dune::XT::Common::FloatCmp::Style::defaultStyle>::value())
{
  DXTC_EXPECT_FLOAT_EQ<Dune::XT::Common::FloatCmp::Style::defaultStyle>(first, second, rtol, atol);
}

template <Dune::XT::Common::FloatCmp::Style style,
          class FirstType,
          class SecondType,
          class ToleranceType = typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps>
typename std::enable_if<
    Dune::XT::Common::FloatCmp::internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value,
    void>::type
DXTC_EXPECT_FLOAT_NE(const FirstType& expected,
                     const SecondType& actual,
                     const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
                         Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value(),
                     const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
                         Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value())
{
  const auto ne = Dune::XT::Common::FloatCmp::internal::
      Call<FirstType, SecondType, typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps, style>::ne;
  EXPECT_PRED4(ne, expected, actual, rtol, atol);
}

template <class FirstType,
          class SecondType,
          class ToleranceType = typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps>
typename std::enable_if<
    Dune::XT::Common::FloatCmp::internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value,
    void>::type
DXTC_EXPECT_FLOAT_NE(
    const FirstType& first,
    const SecondType& second,
    const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
        Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType,
                                                   Dune::XT::Common::FloatCmp::Style::defaultStyle>::value(),
    const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
        Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType,
                                                   Dune::XT::Common::FloatCmp::Style::defaultStyle>::value())
{
  DXTC_EXPECT_FLOAT_NE<Dune::XT::Common::FloatCmp::Style::defaultStyle>(first, second, rtol, atol);
}

template <Dune::XT::Common::FloatCmp::Style style,
          class FirstType,
          class SecondType,
          class ToleranceType = typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps>
typename std::enable_if<
    Dune::XT::Common::FloatCmp::internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value,
    void>::type
DXTC_EXPECT_FLOAT_GT(const FirstType& expected,
                     const SecondType& actual,
                     const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
                         Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value(),
                     const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
                         Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value())
{
  const auto gt = Dune::XT::Common::FloatCmp::internal::
      Call<FirstType, SecondType, typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps, style>::gt;
  EXPECT_PRED4(gt, expected, actual, rtol, atol);
}

template <class FirstType,
          class SecondType,
          class ToleranceType = typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps>
typename std::enable_if<
    Dune::XT::Common::FloatCmp::internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value,
    void>::type
DXTC_EXPECT_FLOAT_GT(
    const FirstType& first,
    const SecondType& second,
    const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
        Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType,
                                                   Dune::XT::Common::FloatCmp::Style::defaultStyle>::value(),
    const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
        Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType,
                                                   Dune::XT::Common::FloatCmp::Style::defaultStyle>::value())
{
  DXTC_EXPECT_FLOAT_GT<Dune::XT::Common::FloatCmp::Style::defaultStyle>(first, second, rtol, atol);
}

template <Dune::XT::Common::FloatCmp::Style style,
          class FirstType,
          class SecondType,
          class ToleranceType = typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps>
typename std::enable_if<
    Dune::XT::Common::FloatCmp::internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value,
    void>::type
DXTC_EXPECT_FLOAT_LT(const FirstType& expected,
                     const SecondType& actual,
                     const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
                         Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value(),
                     const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
                         Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value())
{
  const auto lt = Dune::XT::Common::FloatCmp::internal::
      Call<FirstType, SecondType, typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps, style>::lt;
  EXPECT_PRED4(lt, expected, actual, rtol, atol);
}

template <class FirstType,
          class SecondType,
          class ToleranceType = typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps>
typename std::enable_if<
    Dune::XT::Common::FloatCmp::internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value,
    void>::type
DXTC_EXPECT_FLOAT_LT(
    const FirstType& first,
    const SecondType& second,
    const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
        Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType,
                                                   Dune::XT::Common::FloatCmp::Style::defaultStyle>::value(),
    const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
        Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType,
                                                   Dune::XT::Common::FloatCmp::Style::defaultStyle>::value())
{
  DXTC_EXPECT_FLOAT_LT<Dune::XT::Common::FloatCmp::Style::defaultStyle>(first, second, rtol, atol);
}

template <Dune::XT::Common::FloatCmp::Style style,
          class FirstType,
          class SecondType,
          class ToleranceType = typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps>
typename std::enable_if<
    Dune::XT::Common::FloatCmp::internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value,
    void>::type
DXTC_EXPECT_FLOAT_GE(const FirstType& expected,
                     const SecondType& actual,
                     const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
                         Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value(),
                     const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
                         Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value())
{
  const auto ge = Dune::XT::Common::FloatCmp::internal::
      Call<FirstType, SecondType, typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps, style>::ge;
  EXPECT_PRED4(ge, expected, actual, rtol, atol);
}

template <class FirstType,
          class SecondType,
          class ToleranceType = typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps>
typename std::enable_if<
    Dune::XT::Common::FloatCmp::internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value,
    void>::type
DXTC_EXPECT_FLOAT_GE(
    const FirstType& first,
    const SecondType& second,
    const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
        Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType,
                                                   Dune::XT::Common::FloatCmp::Style::defaultStyle>::value(),
    const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
        Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType,
                                                   Dune::XT::Common::FloatCmp::Style::defaultStyle>::value())
{
  DXTC_EXPECT_FLOAT_GE<Dune::XT::Common::FloatCmp::Style::defaultStyle>(first, second, rtol, atol);
}

template <Dune::XT::Common::FloatCmp::Style style,
          class FirstType,
          class SecondType,
          class ToleranceType = typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps>
typename std::enable_if<
    Dune::XT::Common::FloatCmp::internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value,
    void>::type
DXTC_EXPECT_FLOAT_LE(const FirstType& expected,
                     const SecondType& actual,
                     const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
                         Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value(),
                     const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
                         Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType, style>::value())
{
  const auto le = Dune::XT::Common::FloatCmp::internal::
      Call<FirstType, SecondType, typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps, style>::le;
  EXPECT_PRED4(le, expected, actual, rtol, atol);
}

template <class FirstType,
          class SecondType,
          class ToleranceType = typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps>
typename std::enable_if<
    Dune::XT::Common::FloatCmp::internal::cmp_type_check<FirstType, SecondType, ToleranceType>::value,
    void>::type
DXTC_EXPECT_FLOAT_LE(
    const FirstType& first,
    const SecondType& second,
    const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& rtol =
        Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType,
                                                   Dune::XT::Common::FloatCmp::Style::defaultStyle>::value(),
    const typename Dune::XT::Common::FloatCmp::MT<FirstType>::Eps& atol =
        Dune::XT::Common::FloatCmp::DefaultEpsilon<ToleranceType,
                                                   Dune::XT::Common::FloatCmp::Style::defaultStyle>::value())
{
  DXTC_EXPECT_FLOAT_LE<Dune::XT::Common::FloatCmp::Style::defaultStyle>(first, second, rtol, atol);
}
