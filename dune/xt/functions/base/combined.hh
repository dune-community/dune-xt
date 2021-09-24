// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2018, 2020)
//   René Fritze     (2014 - 2018, 2020)
//   Tim Keil        (2018)
//   Tobias Leibner  (2014, 2017)

#ifndef DUNE_XT_FUNCTIONS_BASE_COMBINED_HH
#define DUNE_XT_FUNCTIONS_BASE_COMBINED_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/parameter.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/functions/exceptions.hh>
#include <dune/xt/functions/type_traits.hh>

namespace Dune::XT::Functions::internal {


template <typename comb, size_t L_r, size_t L_rC, size_t R_r, size_t R_rC>
struct CombinedDim
{
  static constexpr bool available()
  {
    return true;
  };
  static constexpr size_t r()
  {
    return L_r;
  }
  static constexpr size_t rC()
  {
    return L_rC;
  }
};

// general case: non-scalar matrix * non-scalar {matrix or vector}
template <size_t L_r, size_t L_rC, size_t R_r, size_t R_rC>
struct CombinedDim<CombinationType::product, L_r, L_rC, R_r, R_rC>
{
  static constexpr size_t r()
  {
    if constexpr (R_rC == 1 && L_rC == 1 && L_r != 1) {
      return 1;
    } else if constexpr (L_r == 1 && L_rC == 1 && R_r * R_rC != 1) {
      return R_r;
    } else if constexpr (R_r == 1 && R_rC == 1) {
      return L_r;
    } else if constexpr (L_rC * R_rC != 1) {
      return L_r;
    }
    return 0;
  }

  static constexpr size_t rC()
  {
    if constexpr (R_rC == 1 && L_rC == 1 && L_r != 1) {
      return 1;
    } else if constexpr (L_r == 1 && L_rC == 1 && R_r * R_rC != 1) {
      return R_rC;
    } else if constexpr (R_r == 1 && R_rC == 1) {
      return L_rC;
    } else if constexpr (L_rC * R_rC != 1) {
      return R_rC;
    }
    return 0;
  }

  static constexpr bool available()
  {
    return r() != 0 && rC() != 0;
  }
};

template <class L, class R>
static int
compute_combined_order(const L& left, const R& right, const Common::Parameter& param, CombinationType::difference)
{
  return std::max(left.order(param), right.order(param));
}

template <class L, class R>
static int
compute_combined_order(const L& left, const R& right, const Common::Parameter& param, CombinationType::fraction)
{
  return left.order(param) + right.order(param);
}
template <class L, class R>
static int
compute_combined_order(const L& left, const R& right, const Common::Parameter& param, CombinationType::product)
{
  return left.order(param) + right.order(param);
}
template <class L, class R>
static int compute_combined_order(const L& left, const R& right, const Common::Parameter& param, CombinationType::sum)
{
  return std::max(left.order(param), right.order(param));
}

template <class Left, class Right, class D>
static auto compute_combined_eval(
    const Left& left, const Right& right, const D& point, const Common::Parameter& param, CombinationType::difference)
{
  return left.evaluate(point, param) - right.evaluate(point, param);
}

template <class Left, class Right, class D>
static auto compute_combined_eval(
    const Left& left, const Right& right, const D& point, const Common::Parameter& param, CombinationType::fraction)
{
  auto left_value = left.evaluate(point, param);
  const auto right_value = right.evaluate(point, param)[0];
  left_value /= right_value;
  return left_value;
}

template <class Left, class Right, class D>
static auto compute_combined_eval(
    const Left& left, const Right& right, const D& point, const Common::Parameter& param, CombinationType::sum)
{
  return left.evaluate(point, param) + right.evaluate(point, param);
}

template <typename comb, class L_R, class R_R, size_t L_r, size_t L_rC, size_t R_r, size_t R_rC>
struct CombinedEval
{
  using R = typename std::conditional<std::is_same<CombinationType::fraction, comb>::value,
                                      typename Common::multiplication_promotion<L_R, R_R>::type,
                                      typename Common::plus_promotion<L_R, R_R>::type>::type;
  static constexpr size_t r = CombinedDim<comb, L_r, L_rC, R_r, R_rC>::r();
  static constexpr size_t rC = CombinedDim<comb, L_r, L_rC, R_r, R_rC>::rC();

  template <class Left, class Right, class D>
  static auto compute(const Left& left, const Right& right, const D& point, const Common::Parameter& param)
  {
    return compute_combined_eval(left, right, point, param, comb{});
  }
};

template <class L_R, class R_R, size_t L_r, size_t L_rC, size_t R_r, size_t R_rC>
struct CombinedEval<CombinationType::product, L_R, R_R, L_r, L_rC, R_r, R_rC>
{
  static constexpr size_t r = CombinedDim<CombinationType::product, L_r, L_rC, R_r, R_rC>::r();
  static constexpr size_t rC = CombinedDim<CombinationType::product, L_r, L_rC, R_r, R_rC>::rC();
  using R = typename Common::multiplication_promotion<L_R, R_R>::type;
  using RangeReturnType = typename RangeTypeSelector<R, r, rC>::return_type;
  template <class Left, class Right, class D>
  static RangeReturnType compute(const Left& left, const Right& right, const D& point, const Common::Parameter& param)
  {
    // general case: non-scalar matrix * non-scalar {matrix or vector}
    if constexpr (R_r == L_rC && L_rC * R_rC != 1) {
      if constexpr (L_rC == 1) {
        // need to wrap the left vector into a matrix, for * to do the right thing
        XT::Common::FieldMatrix<R, 1, L_r> left_value;
        left_value[0] = left.evaluate(point, param);
        return left_value.transpose() * right.evaluate(point, param);
      } else
        return left.evaluate(point, param) * right.evaluate(point, param);
    }

    // special case: non-scalar elementwise multiplication by scalar from the left
    else if constexpr (L_r == 1 && L_rC == 1 && R_r * R_rC != 1) {
      const auto left_value = left.evaluate(point, param)[0];
      RangeReturnType right_value = right.evaluate(point, param);
      right_value *= left_value;
      return right_value;
    }
    // special case: elementwise multiplication by scalar from the right AND all scalar case
    else if constexpr (std::is_convertible<decltype(left.evaluate(point, param)), RangeReturnType>::value && R_r == 1
                       && R_rC == 1) {
      RangeReturnType left_value = left.evaluate(point, param);
      const auto right_value = right.evaluate(point, param)[0];
      left_value *= right_value;
      return left_value;
    }
    // special case: non-scalar vectors
    else if constexpr (L_r == R_r && R_rC == 1 && L_rC == 1 && L_r != 1) {
      return left.evaluate(point, param) * right.evaluate(point, param);
    }
    DUNE_THROW(NotImplemented, "no product compute combination available");
  }
};

template <class comb,
          class L_R,
          class R_R,
          size_t d,
          size_t L_r,
          size_t L_rC,
          size_t R_r,
          size_t R_rC,
          typename anything = void>
struct CombinedJac
{
  static constexpr size_t r = CombinedDim<comb, L_r, L_rC, R_r, R_rC>::r();
  static constexpr size_t rC = CombinedDim<comb, L_r, L_rC, R_r, R_rC>::rC();

  using R = typename CombinedEval<comb, L_R, R_R, L_r, L_rC, R_r, R_rC>::R;
  using DerivativeRangeReturnType = typename DerivativeRangeTypeSelector<d, R, r, rC>::return_type;

  template <class Left, class Right, class D>
  static DerivativeRangeReturnType
  compute(const Left& /*left*/, const Right& /*right*/, const D& /*point*/, const Common::Parameter& /*param*/)
  {
    DUNE_THROW(Exceptions::combined_error,
               "Not available for a " << get_combination_name(comb{}) << "of " << L_r << "x" << L_rC << " and " << R_r
                                      << "x" << R_rC << "!");
    return DerivativeRangeReturnType();
  }
};

template <class L_R, class R_R, size_t d, size_t L_r, size_t L_rC, typename a>
struct CombinedJac<CombinationType::difference, L_R, R_R, d, L_r, L_rC, L_r, L_rC, a>
{
  static constexpr size_t r = CombinedDim<CombinationType::difference, L_r, L_rC, L_r, L_rC>::r();
  static constexpr size_t rC = CombinedDim<CombinationType::difference, L_r, L_rC, L_r, L_rC>::rC();

  using R = typename CombinedEval<CombinationType::difference, L_R, R_R, L_r, L_rC, L_r, L_rC>::R;
  using DerivativeRangeReturnType = typename DerivativeRangeTypeSelector<d, R, r, rC>::return_type;

  template <class Left, class Right, class D>
  static DerivativeRangeReturnType
  compute(const Left& left, const Right& right, const D& point, const Common::Parameter& param)
  {
    return left.jacobian(point, param) - right.jacobian(point, param);
  }
};

template <class L_R, class R_R, size_t d, size_t L_r, size_t L_rC, typename a>
struct CombinedJac<CombinationType::sum, L_R, R_R, d, L_r, L_rC, L_r, L_rC, a>
{
  static constexpr size_t r = CombinedDim<CombinationType::sum, L_r, L_rC, L_r, L_rC>::r();
  static constexpr size_t rC = CombinedDim<CombinationType::sum, L_r, L_rC, L_r, L_rC>::rC();

  using R = typename CombinedEval<CombinationType::sum, L_R, R_R, L_r, L_rC, L_r, L_rC>::R;
  using DerivativeRangeReturnType = typename DerivativeRangeTypeSelector<d, R, r, rC>::return_type;

  template <class Left, class Right, class D>
  static DerivativeRangeReturnType
  compute(const Left& left, const Right& right, const D& point, const Common::Parameter& param)
  {
    return left.jacobian(point, param) - right.jacobian(point, param);
  }
};


/**
 * \brief Helper class defining types of combined functions, if available.
 *
 * \note Most likely you do not want to use this class directly, but CombinedConstElementFunction or
 *       CombinedElementFunction.
 */
template <class Left, class Right, typename comb>
struct CombinedHelper
{
  static_assert(is_element_function<Left>::value || is_function<Left>::value || is_grid_function<Left>::value);
  static_assert(is_element_function<Right>::value || is_function<Right>::value || is_grid_function<Right>::value);
  static_assert(Left::d == Right::d);

  static constexpr size_t d = Left::d;
  static constexpr size_t L_r = Left::r;
  static constexpr size_t L_rC = Left::rC;
  static constexpr size_t R_r = Right::r;
  static constexpr size_t R_rC = Right::rC;

  using CombinedDimHelper = CombinedDim<comb, L_r, L_rC, R_r, R_rC>;
  static constexpr bool available = CombinedDimHelper::available();
  static constexpr size_t r = CombinedDimHelper::r();
  static constexpr size_t rC = CombinedDimHelper::rC();

  static int order(const Left& left, const Right& right, const Common::Parameter& param)
  {
    return compute_combined_order(left, right, param, comb{});
  }

  using L_R = typename Left::R;
  using R_R = typename Right::R;
  using CombinedEvalHelper = CombinedEval<comb, L_R, R_R, L_r, L_rC, R_r, R_rC>;
  using R = typename CombinedEvalHelper::R;

  template <class D>
  static auto evaluate(const Left& left, const Right& right, const D& point, const Common::Parameter& param)
  {
    return CombinedEvalHelper::compute(left, right, point, param);
  }

  using CombinedJacHelper = CombinedJac<comb, L_R, R_R, d, L_r, L_rC, R_r, R_rC>;
  using DerivativeRangeReturnType = typename CombinedJacHelper::DerivativeRangeReturnType;

  template <class D>
  static DerivativeRangeReturnType
  jacobian(const Left& left, const Right& right, const D& point, const Common::Parameter& param)
  {
    return CombinedJacHelper::compute(left, right, point, param);
  }
}; // struct CombinedHelper


template <class LeftType, class RightType>
struct CombinedStorageProvider
{
  using ThisType = CombinedStorageProvider;

  XT::Common::StorageProvider<LeftType> left;
  XT::Common::StorageProvider<RightType> right;

  CombinedStorageProvider(LeftType& lft, RightType& rght)
    : left(lft)
    , right(rght)
  {}

  CombinedStorageProvider(std::shared_ptr<LeftType> lft, std::shared_ptr<RightType> rght)
    : left(lft)
    , right(rght)
  {}

  CombinedStorageProvider(std::unique_ptr<LeftType>&& lft, std::unique_ptr<RightType>&& rght)
    : left(std::move(lft))
    , right(std::move(rght))
  {}

  CombinedStorageProvider(const ThisType&) = default;

  CombinedStorageProvider(ThisType&) = default;

  CombinedStorageProvider(ThisType&&) = default;
}; // struct CombinedStorageProvider


} // namespace Dune::XT::Functions::internal

#endif // DUNE_XT_FUNCTIONS_BASE_COMBINED_HH
