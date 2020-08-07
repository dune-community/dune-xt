// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2016 - 2017)
//   René Fritze      (2015 - 2019)
//   Tim Keil         (2018)
//   Tobias Leibner   (2017 - 2020)

#ifndef DUNE_XT_COMMON_FLOAT_CMP_INTERNAL_HH
#define DUNE_XT_COMMON_FLOAT_CMP_INTERNAL_HH

#include <type_traits>
#include <complex>

#include <dune/xt/common/math.hh> // <- This include needs to be before the one from dune-common, otherwise
#include <dune/common/float_cmp.hh> //  std::abs(long unsinged int) is indefined in dune-common!

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/float_cmp_style.hh>

namespace Dune {
namespace FloatCmp {
template <class T>
struct EpsilonType<std::complex<T>>
{
  //! The epsilon type corresponding to value type std::complex<T>
  using Type = std::complex<T>;
};
} // namespace FloatCmp
namespace XT {
namespace Common {

// forward, include is below
template <class VecType>
struct VectorAbstraction;

template <class MatType>
struct MatrixAbstraction;

namespace FloatCmp {
namespace internal {

/**
 * Note: The organization of this file is not optimal. Currently we require a non dune_ based method for the numpy
 * compare style and a dune_ based method elsewhere. Furthermore, we have 4 specializations for each compare function.
 */

/**
 * Taken from http://docs.pymor.org/en/latest/_modules/pymor/tools/floatcmp.html#float_cmp,
 *            https://github.com/numpy/numpy/blob/v1.9.1/numpy/core/numeric.py#L2238
 */
template <class T>
typename std::enable_if<is_arithmetic<T>::value, bool>::type
float_cmp_eq(const T& xx, const T& yy, const T& rtol, const T& atol)
{
  return absolute_difference(xx, yy) <= atol + Dune::XT::Common::abs(yy) * rtol;
}

template <class T>
typename std::enable_if<is_complex<T>::value, bool>::type
float_cmp_eq(const T& xx, const T& yy, const T& rtol, const T& atol)
{
  using namespace std;
  return (float_cmp_eq(real(xx), real(yy), real(rtol), real(atol))
          && float_cmp_eq(imag(xx), imag(yy), imag(rtol), imag(atol)));
}

template <class XType, class YType, class TolType>
typename std::enable_if<is_vector<XType>::value && is_vector<YType>::value
                            && std::is_same<typename VectorAbstraction<XType>::S, TolType>::value
                            && std::is_same<typename VectorAbstraction<YType>::S, TolType>::value,
                        bool>::type
float_cmp_eq(const XType& xx, const YType& yy, const TolType& rtol, const TolType& atol)
{
  const auto sz = xx.size();
  if (yy.size() != sz)
    return false;
  for (size_t ii = 0; ii < sz; ++ii)
    if (!float_cmp_eq(
            VectorAbstraction<XType>::get_entry(xx, ii), VectorAbstraction<YType>::get_entry(yy, ii), rtol, atol))
      return false;
  return true;
} // ... float_cmp(...)

template <class XType, class YType, class TolType>
typename std::enable_if<is_matrix<XType>::value && is_matrix<YType>::value
                            && std::is_same<typename MatrixAbstraction<XType>::S, TolType>::value
                            && std::is_same<typename MatrixAbstraction<YType>::S, TolType>::value,
                        bool>::type
float_cmp_eq(const XType& xx, const YType& yy, const TolType& rtol, const TolType& atol)
{
  const auto rows = MatrixAbstraction<XType>::rows(xx);
  const auto cols = MatrixAbstraction<XType>::cols(xx);
  if (MatrixAbstraction<YType>::rows(yy) != rows || MatrixAbstraction<YType>::cols(yy) != cols)
    return false;
  for (size_t ii = 0; ii < rows; ++ii)
    for (size_t jj = 0; jj < cols; ++jj)
      if (!float_cmp_eq(MatrixAbstraction<XType>::get_entry(xx, ii, jj),
                        MatrixAbstraction<YType>::get_entry(yy, ii, jj),
                        rtol,
                        atol))
        return false;
  return true;
} // ... float_cmp(...)

template <Dune::FloatCmp::CmpStyle style, class T>
typename std::enable_if<is_arithmetic<T>::value, bool>::type dune_float_cmp_eq(const T& xx, const T& yy, const T& eps)
{
  return Dune::FloatCmp::eq<T, style>(xx, yy, eps);
}

template <Dune::FloatCmp::CmpStyle style, class T>
typename std::enable_if<is_complex<T>::value, bool>::type dune_float_cmp_eq(const T& xx, const T& yy, const T& eps)
{
  using namespace std;
  return Dune::FloatCmp::eq<typename T::value_type, style>(real(xx), real(yy), real(eps))
         && Dune::FloatCmp::eq<typename T::value_type, style>(imag(xx), imag(yy), imag(eps));
}

template <Dune::FloatCmp::CmpStyle style, class XType, class YType, class EpsType>
typename std::enable_if<is_vector<XType>::value && is_vector<YType>::value
                            && std::is_same<typename VectorAbstraction<XType>::S, EpsType>::value
                            && std::is_same<typename VectorAbstraction<YType>::S, EpsType>::value,
                        bool>::type
dune_float_cmp_eq(const XType& xx, const YType& yy, const EpsType& eps)
{
  const auto sz = xx.size();
  if (yy.size() != sz)
    return false;
  for (size_t ii = 0; ii < sz; ++ii)
    if (!dune_float_cmp_eq<style>(xx[ii], yy[ii], eps))
      return false;
  return true;
} // ... dune_float_cmp(...)

template <Dune::FloatCmp::CmpStyle style, class XType, class YType, class EpsType>
typename std::enable_if<is_matrix<XType>::value && is_matrix<YType>::value
                            && std::is_same<typename MatrixAbstraction<XType>::S, EpsType>::value
                            && std::is_same<typename MatrixAbstraction<YType>::S, EpsType>::value,
                        bool>::type
dune_float_cmp_eq(const XType& xx, const YType& yy, const EpsType& eps)
{
  const auto rows = MatrixAbstraction<XType>::rows(xx);
  const auto cols = MatrixAbstraction<XType>::cols(xx);
  if (MatrixAbstraction<YType>::rows(yy) != rows || MatrixAbstraction<YType>::cols(yy) != cols)
    return false;
  for (size_t ii = 0; ii < rows; ++ii)
    for (size_t jj = 0; jj < cols; ++jj)
      if (!dune_float_cmp_eq<style>(
              MatrixAbstraction<XType>::get_entry(xx, ii, jj), MatrixAbstraction<YType>::get_entry(yy, ii, jj), eps))
        return false;
  return true;
} // ... float_cmp(...){

template <class T>
typename std::enable_if<is_arithmetic<T>::value, bool>::type cmp_gt(const T& xx, const T& yy)
{
  return std::greater<T>()(xx, yy);
}

template <class T>
bool cmp_gt(const std::complex<T>& xx, const std::complex<T>& yy)
{
  using namespace std;
  return cmp_gt(real(xx), real(yy)) && cmp_gt(imag(xx), imag(yy));
}

template <class XType, class YType>
typename std::enable_if<
    is_vector<XType>::value && is_vector<YType>::value
        && std::is_same<typename VectorAbstraction<XType>::S, typename VectorAbstraction<YType>::S>::value,
    bool>::type
cmp_gt(const XType& xx, const YType& yy)
{
  const auto sz = xx.size();
  if (yy.size() != sz)
    return false;
  for (size_t ii = 0; ii < sz; ++ii)
    if (!(cmp_gt(xx[ii], yy[ii])))
      return false;
  return true;
} // ... cmp_gt(...)

template <class XType, class YType>
typename std::enable_if<
    is_matrix<XType>::value && is_matrix<YType>::value
        && std::is_same<typename MatrixAbstraction<XType>::S, typename MatrixAbstraction<YType>::S>::value,
    bool>::type
cmp_gt(const XType& xx, const YType& yy)
{
  const auto rows = MatrixAbstraction<XType>::rows(xx);
  const auto cols = MatrixAbstraction<XType>::cols(xx);
  if (MatrixAbstraction<YType>::rows(yy) != rows || MatrixAbstraction<YType>::cols(yy) != cols)
    return false;
  for (size_t ii = 0; ii < rows; ++ii)
    for (size_t jj = 0; jj < cols; ++jj)
      if (!cmp_gt(MatrixAbstraction<XType>::get_entry(xx, ii, jj), MatrixAbstraction<YType>::get_entry(yy, ii, jj)))
        return false;
  return true;
} // ... cmp_gt(...)

template <class T>
typename std::enable_if<is_arithmetic<T>::value, bool>::type cmp_lt(const T& xx, const T& yy)
{
  return std::less<T>()(xx, yy);
}

template <class T>
typename std::enable_if<is_complex<T>::value, bool>::type cmp_lt(const T& xx, const T& yy)
{
  return cmp_lt(real(xx), real(yy)) && cmp_lt(imag(xx), imag(yy));
}

template <class XType, class YType>
typename std::enable_if<
    is_vector<XType>::value && is_vector<YType>::value
        && std::is_same<typename VectorAbstraction<XType>::S, typename VectorAbstraction<YType>::S>::value,
    bool>::type
cmp_lt(const XType& xx, const YType& yy)
{
  const auto sz = xx.size();
  if (yy.size() != sz)
    return false;
  for (size_t ii = 0; ii < sz; ++ii)
    if (!(cmp_lt(xx[ii], yy[ii])))
      return false;
  return true;
} // ... cmp_lt(...)

template <class XType, class YType>
typename std::enable_if<
    is_matrix<XType>::value && is_matrix<YType>::value
        && std::is_same<typename MatrixAbstraction<XType>::S, typename MatrixAbstraction<YType>::S>::value,
    bool>::type
cmp_lt(const XType& xx, const YType& yy)
{
  const auto rows = MatrixAbstraction<XType>::rows(xx);
  const auto cols = MatrixAbstraction<XType>::cols(xx);
  if (MatrixAbstraction<YType>::rows(yy) != rows || MatrixAbstraction<YType>::cols(yy) != cols)
    return false;
  for (size_t ii = 0; ii < rows; ++ii)
    for (size_t jj = 0; jj < cols; ++jj)
      if (!cmp_lt(MatrixAbstraction<XType>::get_entry(xx, ii, jj), MatrixAbstraction<YType>::get_entry(yy, ii, jj)))
        return false;
  return true;
} // ... cmp_gt(...)

template <Dune::FloatCmp::CmpStyle style, class T>
typename std::enable_if<is_arithmetic<T>::value, bool>::type dune_cmp_ge(const T& xx, const T& yy, const T& eps)
{
  return std::greater<T>()(xx, yy) || Dune::FloatCmp::eq<T, style>(xx, yy, eps);
}

template <Dune::FloatCmp::CmpStyle style, class T>
typename std::enable_if<is_complex<T>::value, bool>::type dune_cmp_ge(const T& xx, const T& yy, const T& eps)
{
  return dune_cmp_ge<style>(real(xx), real(yy), real(eps)) && dune_cmp_ge<style>(imag(xx), imag(yy), imag(eps));
}

template <Dune::FloatCmp::CmpStyle style, class XType, class YType, class EpsType>
typename std::enable_if<is_vector<XType>::value && is_vector<YType>::value
                            && std::is_same<typename VectorAbstraction<XType>::S, EpsType>::value
                            && std::is_same<typename VectorAbstraction<YType>::S, EpsType>::value,
                        bool>::type
dune_cmp_ge(const XType& xx, const YType& yy, const EpsType& eps)
{
  const auto sz = xx.size();
  if (yy.size() != sz)
    return false;
  for (size_t ii = 0; ii < sz; ++ii)
    if (!(dune_cmp_ge<style>(xx[ii], yy[ii], eps)))
      return false;
  return true;
} // ... dune_cmp_ge(...)

template <Dune::FloatCmp::CmpStyle style, class XType, class YType, class EpsType>
typename std::enable_if<is_matrix<XType>::value && is_matrix<YType>::value
                            && std::is_same<typename MatrixAbstraction<XType>::S, EpsType>::value
                            && std::is_same<typename MatrixAbstraction<YType>::S, EpsType>::value,
                        bool>::type
dune_cmp_ge(const XType& xx, const YType& yy, const EpsType& eps)
{
  const auto rows = MatrixAbstraction<XType>::rows(xx);
  const auto cols = MatrixAbstraction<XType>::cols(xx);
  if (MatrixAbstraction<YType>::rows(yy) != rows || MatrixAbstraction<YType>::cols(yy) != cols)
    return false;
  for (size_t ii = 0; ii < rows; ++ii)
    for (size_t jj = 0; jj < cols; ++jj)
      if (!dune_cmp_ge<style>(
              MatrixAbstraction<XType>::get_entry(xx, ii, jj), MatrixAbstraction<YType>::get_entry(yy, ii, jj), eps))
        return false;
  return true;
} // ... dune_cmp_ge(...)

template <Dune::FloatCmp::CmpStyle style, class T>
typename std::enable_if<is_arithmetic<T>::value, bool>::type dune_cmp_le(const T& xx, const T& yy, const T& eps)
{
  return std::less<T>()(xx, yy) || Dune::FloatCmp::eq<T, style>(xx, yy, eps);
}

template <Dune::FloatCmp::CmpStyle style, class T>
typename std::enable_if<is_complex<T>::value, bool>::type dune_cmp_le(const T& xx, const T& yy, const T& eps)
{
  return dune_cmp_le<style>(real(xx), real(yy), real(eps)) && dune_cmp_le<style>(imag(xx), imag(yy), imag(eps));
}

template <Dune::FloatCmp::CmpStyle style, class XType, class YType, class EpsType>
typename std::enable_if<is_vector<XType>::value && is_vector<YType>::value
                            && std::is_same<typename VectorAbstraction<XType>::S, EpsType>::value
                            && std::is_same<typename VectorAbstraction<YType>::S, EpsType>::value,
                        bool>::type
dune_cmp_le(const XType& xx, const YType& yy, const EpsType& eps)
{
  const auto sz = xx.size();
  if (yy.size() != sz)
    return false;
  for (size_t ii = 0; ii < sz; ++ii)
    if (!(dune_cmp_le<style>(xx[ii], yy[ii], eps)))
      return false;
  return true;
} // ... dune_cmp_le(...)

template <Dune::FloatCmp::CmpStyle style, class XType, class YType, class EpsType>
typename std::enable_if<is_matrix<XType>::value && is_matrix<YType>::value
                            && std::is_same<typename MatrixAbstraction<XType>::S, EpsType>::value
                            && std::is_same<typename MatrixAbstraction<YType>::S, EpsType>::value,
                        bool>::type
dune_cmp_le(const XType& xx, const YType& yy, const EpsType& eps)
{
  const auto rows = MatrixAbstraction<XType>::rows(xx);
  const auto cols = MatrixAbstraction<XType>::cols(xx);
  if (MatrixAbstraction<YType>::rows(yy) != rows || MatrixAbstraction<YType>::cols(yy) != cols)
    return false;
  for (size_t ii = 0; ii < rows; ++ii)
    for (size_t jj = 0; jj < cols; ++jj)
      if (!dune_cmp_le<style>(
              MatrixAbstraction<XType>::get_entry(xx, ii, jj), MatrixAbstraction<YType>::get_entry(yy, ii, jj), eps))
        return false;
  return true;
} // ... dune_cmp_le(...)


template <class T>
typename std::enable_if<is_arithmetic<T>::value, bool>::type
cmp_ge(const T& xx, const T& yy, const T& rtol, const T& atol)
{
  return std::greater<T>()(xx, yy) || float_cmp_eq<T>(xx, yy, rtol, atol);
}


template <class T>
typename std::enable_if<is_complex<T>::value, bool>::type cmp_ge(const T& xx, const T& yy, const T& rtol, const T& atol)
{
  return cmp_ge(real(xx), real(yy), real(rtol), real(atol)) && cmp_ge(imag(xx), imag(yy), imag(rtol), imag(atol));
}

template <class XType, class YType, class TolType>
typename std::enable_if<is_vector<XType>::value && is_vector<YType>::value
                            && std::is_same<typename VectorAbstraction<XType>::S, TolType>::value
                            && std::is_same<typename VectorAbstraction<YType>::S, TolType>::value,
                        bool>::type
cmp_ge(const XType& xx, const YType& yy, const TolType& rtol, const TolType& atol)
{
  const auto sz = xx.size();
  if (yy.size() != sz)
    return false;
  for (size_t ii = 0; ii < sz; ++ii)
    if (!(cmp_ge(xx[ii], yy[ii], rtol, atol)))
      return false;
  return true;
} // ... cmp_ge(...)

template <class XType, class YType, class TolType>
typename std::enable_if<is_matrix<XType>::value && is_matrix<YType>::value
                            && std::is_same<typename MatrixAbstraction<XType>::S, TolType>::value
                            && std::is_same<typename MatrixAbstraction<YType>::S, TolType>::value,
                        bool>::type
cmp_ge(const XType& xx, const YType& yy, const TolType& rtol, const TolType& atol)
{
  const auto rows = MatrixAbstraction<XType>::rows(xx);
  const auto cols = MatrixAbstraction<XType>::cols(xx);
  if (MatrixAbstraction<YType>::rows(yy) != rows || MatrixAbstraction<YType>::cols(yy) != cols)
    return false;
  for (size_t ii = 0; ii < rows; ++ii)
    for (size_t jj = 0; jj < cols; ++jj)
      if (!cmp_ge(MatrixAbstraction<XType>::get_entry(xx, ii, jj),
                  MatrixAbstraction<YType>::get_entry(yy, ii, jj),
                  rtol,
                  atol))
        return false;
  return true;
} // ... cmp_ge(...)

template <class T>
typename std::enable_if<is_arithmetic<T>::value, bool>::type
cmp_le(const T& xx, const T& yy, const T& rtol, const T& atol)
{
  return std::less<T>()(xx, yy) || float_cmp_eq<T>(xx, yy, rtol, atol);
}


template <class T>
typename std::enable_if<is_complex<T>::value, bool>::type cmp_le(const T& xx, const T& yy, const T& rtol, const T& atol)
{
  return cmp_le(real(xx), real(yy), real(rtol), real(atol)) && cmp_le(imag(xx), imag(yy), imag(rtol), imag(atol));
}

template <class XType, class YType, class TolType>
typename std::enable_if<is_vector<XType>::value && is_vector<YType>::value
                            && std::is_same<typename VectorAbstraction<XType>::S, TolType>::value
                            && std::is_same<typename VectorAbstraction<YType>::S, TolType>::value,
                        bool>::type
cmp_le(const XType& xx, const YType& yy, const TolType& rtol, const TolType& atol)
{
  const auto sz = xx.size();
  if (yy.size() != sz)
    return false;
  for (size_t ii = 0; ii < sz; ++ii)
    if (!(cmp_le(xx[ii], yy[ii], rtol, atol)))
      return false;
  return true;
} // ... cmp_le(...)

template <class XType, class YType, class TolType>
typename std::enable_if<is_matrix<XType>::value && is_matrix<YType>::value
                            && std::is_same<typename MatrixAbstraction<XType>::S, TolType>::value
                            && std::is_same<typename MatrixAbstraction<YType>::S, TolType>::value,
                        bool>::type
cmp_le(const XType& xx, const YType& yy, const TolType& rtol, const TolType& atol)
{
  const auto rows = MatrixAbstraction<XType>::rows(xx);
  const auto cols = MatrixAbstraction<XType>::cols(xx);
  if (MatrixAbstraction<YType>::rows(yy) != rows || MatrixAbstraction<YType>::cols(yy) != cols)
    return false;
  for (size_t ii = 0; ii < rows; ++ii)
    for (size_t jj = 0; jj < cols; ++jj)
      if (!cmp_le(MatrixAbstraction<XType>::get_entry(xx, ii, jj),
                  MatrixAbstraction<YType>::get_entry(yy, ii, jj),
                  rtol,
                  atol))
        return false;
  return true;
} // ... cmp_le(...)


template <class FirstType, class SecondType, class ToleranceType, Style style>
struct Call
{
  static bool eq(const FirstType& first, const SecondType& second, const ToleranceType& rtol, const ToleranceType& /**/)
  {
    return dune_float_cmp_eq<internal::ConvertStyle<style>::value>(first, second, rtol);
  }

  static bool ne(const FirstType& first, const SecondType& second, const ToleranceType& rtol, const ToleranceType& /**/)
  {
    return !dune_float_cmp_eq<internal::ConvertStyle<style>::value>(first, second, rtol);
  }

  static bool gt(const FirstType& first, const SecondType& second, const ToleranceType& rtol, const ToleranceType& /**/)
  {
    /**
     * todo: we need to check why dune_float_cmp_eq is required here. Normally cmp_gt should suffice.
     */
    return !dune_float_cmp_eq<internal::ConvertStyle<style>::value>(first, second, rtol) && cmp_gt(first, second);
  }

  static bool lt(const FirstType& first, const SecondType& second, const ToleranceType& rtol, const ToleranceType& /**/)
  {
    /**
     * todo: we need to check why dune_float_cmp_eq is required here. Normally cmp_lt should suffice.
     */
    return !dune_float_cmp_eq<internal::ConvertStyle<style>::value>(first, second, rtol) && cmp_lt(first, second);
  }

  static bool ge(const FirstType& first, const SecondType& second, const ToleranceType& rtol, const ToleranceType& /**/)
  {
    return dune_cmp_ge<internal::ConvertStyle<style>::value>(first, second, rtol);
  }

  static bool le(const FirstType& first, const SecondType& second, const ToleranceType& rtol, const ToleranceType& /**/)
  {
    return dune_cmp_le<internal::ConvertStyle<style>::value>(first, second, rtol);
  }
};

template <class FirstType, class SecondType, class ToleranceType>
struct Call<FirstType, SecondType, ToleranceType, Style::numpy>
{
  static bool eq(const FirstType& first, const SecondType& second, const ToleranceType& rtol, const ToleranceType& atol)
  {
    return float_cmp_eq(first, second, rtol, atol);
  }

  static bool ne(const FirstType& first, const SecondType& second, const ToleranceType& rtol, const ToleranceType& atol)
  {
    return !float_cmp_eq(first, second, rtol, atol);
  }

  static bool gt(const FirstType& first, const SecondType& second, const ToleranceType& rtol, const ToleranceType& atol)
  {
    /**
     * todo: we need to check why float_cmp_eq is required here. Normally cmp_gt should suffice.
     */
    return !float_cmp_eq(first, second, rtol, atol) && cmp_gt(first, second);
  }

  static bool lt(const FirstType& first, const SecondType& second, const ToleranceType& rtol, const ToleranceType& atol)
  {
    /**
     * todo: we need to check why float_cmp_eq is required here. Normally cmp_lt should suffice.
     */
    return !float_cmp_eq(first, second, rtol, atol) && cmp_lt(first, second);
  }

  static bool ge(const FirstType& first, const SecondType& second, const ToleranceType& rtol, const ToleranceType& atol)
  {
    return cmp_ge(first, second, rtol, atol);
  }

  static bool le(const FirstType& first, const SecondType& second, const ToleranceType& rtol, const ToleranceType& atol)
  {
    return cmp_le(first, second, rtol, atol);
  }
};

template <class FirstType, class SecondType, class ToleranceType>
struct cmp_type_check
{
  static constexpr bool is_ok_scalar =
      (is_arithmetic<FirstType>::value || is_complex<FirstType>::value) && std::is_same<FirstType, SecondType>::value;
  static constexpr bool is_ok_vector = is_vector<FirstType>::value && is_vector<SecondType>::value
                                       && std::is_same<ToleranceType, typename VectorAbstraction<FirstType>::S>::value
                                       && std::is_same<ToleranceType, typename VectorAbstraction<SecondType>::S>::value;
  static constexpr bool is_ok_matrix = is_matrix<FirstType>::value && is_matrix<SecondType>::value
                                       && std::is_same<ToleranceType, typename MatrixAbstraction<FirstType>::S>::value
                                       && std::is_same<ToleranceType, typename MatrixAbstraction<SecondType>::S>::value;
  static constexpr bool value = is_ok_scalar || is_ok_vector || is_ok_matrix;
};

} // namespace internal
} // namespace FloatCmp
} // namespace Common
} // namespace XT
} // namespace Dune

#include <dune/xt/common/vector.hh>

#endif // DUNE_XT_COMMON_FLOAT_CMP_INTERNAL_HH
