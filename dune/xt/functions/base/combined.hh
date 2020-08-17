// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef DUNE_XT_FUNCTIONS_BASE_COMBINED_HH
#define DUNE_XT_FUNCTIONS_BASE_COMBINED_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/parameter.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/functions/exceptions.hh>
#include <dune/xt/functions/type_traits.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace internal {


template <CombinationType comb, size_t L_r, size_t L_rC, size_t R_r, size_t R_rC, typename anything = void>
struct CombinedDim
{
  static const constexpr bool available = false;
};

template <size_t L_r, size_t L_rC, typename a>
struct CombinedDim<CombinationType::difference, L_r, L_rC, L_r, L_rC, a>
{
  static const constexpr bool available = true;
  static const constexpr size_t r = L_r;
  static const constexpr size_t rC = L_rC;
};

template <size_t L_r, size_t L_rC, typename a>
struct CombinedDim<CombinationType::fraction, L_r, L_rC, 1, 1, a>
{
  static const constexpr bool available = true;
  static const constexpr size_t r = L_r;
  static const constexpr size_t rC = L_rC;
};

// general case: non-scalar matrix * non-scalar {matrix or vector}
template <size_t L_r, size_t L_rC, size_t R_rC>
struct CombinedDim<CombinationType::product, L_r, L_rC, L_rC, R_rC, std::enable_if_t<L_rC * R_rC != 1, void>>
{
  static const constexpr bool available = true;
  static const constexpr size_t r = L_r;
  static const constexpr size_t rC = R_rC;
};

// special case: elementwise multiplication by scalar from the right AND all scalar case
template <size_t L_r, size_t L_rC, typename a>
struct CombinedDim<CombinationType::product, L_r, L_rC, 1, 1, a>
{
  static const constexpr bool available = true;
  static const constexpr size_t r = L_r;
  static const constexpr size_t rC = L_rC;
};

// special case: non-scalar elementwise multiplication by scalar from the left
template <size_t R_r, size_t R_rC>
struct CombinedDim<CombinationType::product, 1, 1, R_r, R_rC, std::enable_if_t<R_r * R_rC != 1, void>>
{
  static const constexpr bool available = true;
  static const constexpr size_t r = R_r;
  static const constexpr size_t rC = R_rC;
};

// special case: non-scalar vectors
template <size_t L_r>
struct CombinedDim<CombinationType::product, L_r, 1, L_r, 1, std::enable_if_t<L_r != 1, void>>
{
  static const constexpr bool available = true;
  static const constexpr size_t r = 1;
  static const constexpr size_t rC = 1;
};

template <size_t L_r, size_t L_rC, typename a>
struct CombinedDim<CombinationType::sum, L_r, L_rC, L_r, L_rC, a>
{
  static const constexpr bool available = true;
  static const constexpr size_t r = L_r;
  static const constexpr size_t rC = L_rC;
};


template <CombinationType comb, size_t L_r, size_t L_rC, size_t R_r, size_t R_rC, typename anything = void>
struct CombinedOrder
{
  static_assert(CombinedDim<comb, L_r, L_rC, R_r, R_rC>::available, "Not available for these dimensions!");
};

template <size_t L_r, size_t L_rC, size_t R_r, size_t R_rC, typename a>
struct CombinedOrder<CombinationType::difference, L_r, L_rC, R_r, R_rC, a>
{
  template <class L, class R>
  static int compute(const L& left, const R& right, const Common::Parameter& param)
  {
    return std::max(left.order(param), right.order(param));
  }
};

/// \todo Any better ideas?
template <size_t L_r, size_t L_rC, size_t R_r, size_t R_rC, typename a>
struct CombinedOrder<CombinationType::fraction, L_r, L_rC, R_r, R_rC, a>
{
  template <class L, class R>
  static int compute(const L& left, const R& right, const Common::Parameter& param)
  {
    return left.order(param) + right.order(param);
  }
};

template <size_t L_r, size_t L_rC, size_t R_r, size_t R_rC, typename a>
struct CombinedOrder<CombinationType::product, L_r, L_rC, R_r, R_rC, a>
{
  template <class L, class R>
  static int compute(const L& left, const R& right, const Common::Parameter& param)
  {
    return left.order(param) + right.order(param);
  }
};

template <size_t L_r, size_t L_rC, size_t R_r, size_t R_rC, typename a>
struct CombinedOrder<CombinationType::sum, L_r, L_rC, R_r, R_rC, a>
{
  template <class L, class R>
  static int compute(const L& left, const R& right, const Common::Parameter& param)
  {
    return std::max(left.order(param), right.order(param));
  }
};


template <CombinationType comb,
          class L_R,
          class R_R,
          size_t L_r,
          size_t L_rC,
          size_t R_r,
          size_t R_rC,
          typename anything = void>
struct CombinedEval
{};

template <class L_R, class R_R, size_t L_r, size_t L_rC, typename a>
struct CombinedEval<CombinationType::difference, L_R, R_R, L_r, L_rC, L_r, L_rC, a>
{
  static const constexpr size_t r = CombinedDim<CombinationType::difference, L_r, L_rC, L_r, L_rC>::r;
  static const constexpr size_t rC = CombinedDim<CombinationType::difference, L_r, L_rC, L_r, L_rC>::rC;

  using R = typename Common::plus_promotion<L_R, R_R>::type;
  using RangeReturnType = typename RangeTypeSelector<R, r, rC>::return_type;

  template <class Left, class Right, class D>
  static RangeReturnType compute(const Left& left, const Right& right, const D& point, const Common::Parameter& param)
  {
    return left.evaluate(point, param) - right.evaluate(point, param);
  }
};

template <class L_R, class R_R, size_t L_r, size_t L_rC, typename a>
struct CombinedEval<CombinationType::fraction, L_R, R_R, L_r, L_rC, 1, 1, a>
{
  static const constexpr size_t r = CombinedDim<CombinationType::fraction, L_r, L_rC, 1, 1>::r;
  static const constexpr size_t rC = CombinedDim<CombinationType::fraction, L_r, L_rC, 1, 1>::rC;

  using R = typename Common::multiplication_promotion<L_R, R_R>::type;
  using RangeReturnType = typename RangeTypeSelector<R, r, rC>::return_type;

  template <class Left, class Right, class D>
  static RangeReturnType compute(const Left& left, const Right& right, const D& point, const Common::Parameter& param)
  {
    RangeReturnType left_value = left.evaluate(point, param);
    auto right_value = right.evaluate(point, param)[0];
    left_value /= right_value;
    return left_value;
  }
};

// general case: non-scalar matrix * non-scalar {matrix or vector}
template <class L_R, class R_R, size_t L_r, size_t L_rC, size_t R_rC>
struct CombinedEval<CombinationType::product, L_R, R_R, L_r, L_rC, L_rC, R_rC, std::enable_if_t<L_rC * R_rC != 1, void>>
{
  static const constexpr size_t r = CombinedDim<CombinationType::product, L_r, L_rC, L_rC, R_rC>::r;
  static const constexpr size_t rC = CombinedDim<CombinationType::product, L_r, L_rC, L_rC, R_rC>::rC;

  using R = typename Common::multiplication_promotion<L_R, R_R>::type;
  using RangeReturnType = typename RangeTypeSelector<R, r, rC>::return_type;

  template <class Left, class Right, class D>
  static RangeReturnType compute(const Left& left, const Right& right, const D& point, const Common::Parameter& param)
  {
    if constexpr (L_rC == 1) {
      // need to wrap the left vector into a matrix, for * to do the right thing
      XT::Common::FieldMatrix<R, 1, L_r> left_value;
      left_value[0] = left.evaluate(point, param);
      return left_value.transpose() * right.evaluate(point, param);
    } else
      return left.evaluate(point, param) * right.evaluate(point, param);
  }
};

// special case: elementwise multiplication by scalar from the right AND all scalar case
template <class L_R, class R_R, size_t L_r, size_t L_rC, typename a>
struct CombinedEval<CombinationType::product, L_R, R_R, L_r, L_rC, 1, 1, a>
{
  static const constexpr size_t r = CombinedDim<CombinationType::product, L_r, L_rC, 1, 1>::r;
  static const constexpr size_t rC = CombinedDim<CombinationType::product, L_r, L_rC, 1, 1>::rC;

  using R = typename Common::multiplication_promotion<L_R, R_R>::type;
  using RangeReturnType = typename RangeTypeSelector<R, r, rC>::return_type;

  template <class Left, class Right, class D>
  static RangeReturnType compute(const Left& left, const Right& right, const D& point, const Common::Parameter& param)
  {
    RangeReturnType left_value = left.evaluate(point, param);
    auto right_value = right.evaluate(point, param)[0];
    left_value *= right_value;
    return left_value;
  }
};

// special case: non-scalar elementwise multiplication by scalar from the left
template <class L_R, class R_R, size_t R_r, size_t R_rC>
struct CombinedEval<CombinationType::product, L_R, R_R, 1, 1, R_r, R_rC, std::enable_if_t<R_r * R_rC != 1, void>>
{
  static const constexpr size_t r = CombinedDim<CombinationType::product, 1, 1, R_r, R_rC>::r;
  static const constexpr size_t rC = CombinedDim<CombinationType::product, 1, 1, R_r, R_rC>::rC;

  using R = typename Common::multiplication_promotion<L_R, R_R>::type;
  using RangeReturnType = typename RangeTypeSelector<R, r, rC>::return_type;

  template <class Left, class Right, class D>
  static RangeReturnType compute(const Left& left, const Right& right, const D& point, const Common::Parameter& param)
  {
    auto left_value = left.evaluate(point, param)[0];
    RangeReturnType right_value = right.evaluate(point, param);
    right_value *= left_value;
    return right_value;
  }
};

// special case: non-scalar vectors
template <class L_R, class R_R, size_t L_r>
struct CombinedEval<CombinationType::product, L_R, R_R, L_r, 1, L_r, 1, std::enable_if_t<L_r != 1, void>>
{
  static const constexpr size_t r = CombinedDim<CombinationType::product, L_r, 1, L_r, 1>::r;
  static const constexpr size_t rC = CombinedDim<CombinationType::product, L_r, 1, L_r, 1>::rC;

  using R = typename Common::multiplication_promotion<L_R, R_R>::type;
  using RangeReturnType = typename RangeTypeSelector<R, r, rC>::return_type;

  template <class Left, class Right, class D>
  static RangeReturnType compute(const Left& left, const Right& right, const D& point, const Common::Parameter& param)
  {
    return left.evaluate(point, param) * right.evaluate(point, param);
  }
};

template <class L_R, class R_R, size_t L_r, size_t L_rC, typename a>
struct CombinedEval<CombinationType::sum, L_R, R_R, L_r, L_rC, L_r, L_rC, a>
{
  static const constexpr size_t r = CombinedDim<CombinationType::sum, L_r, L_rC, L_r, L_rC>::r;
  static const constexpr size_t rC = CombinedDim<CombinationType::sum, L_r, L_rC, L_r, L_rC>::rC;

  using R = typename Common::plus_promotion<L_R, R_R>::type;
  using RangeReturnType = typename RangeTypeSelector<R, r, rC>::return_type;

  template <class Left, class Right, class D>
  static RangeReturnType compute(const Left& left, const Right& right, const D& point, const Common::Parameter& param)
  {
    return left.evaluate(point, param) + right.evaluate(point, param);
  }
};


template <CombinationType comb,
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
  static const constexpr size_t r = CombinedDim<comb, L_r, L_rC, R_r, R_rC>::r;
  static const constexpr size_t rC = CombinedDim<comb, L_r, L_rC, R_r, R_rC>::rC;

  using R = typename CombinedEval<comb, L_R, R_R, L_r, L_rC, R_r, R_rC>::R;
  using DerivativeRangeReturnType = typename DerivativeRangeTypeSelector<d, R, r, rC>::return_type;

  template <class Left, class Right, class D>
  static DerivativeRangeReturnType
  compute(const Left& /*left*/, const Right& /*right*/, const D& /*point*/, const Common::Parameter& /*param*/)
  {
    DUNE_THROW(Exceptions::combined_error,
               "Not available for a " << GetCombination<comb>::name() << "of " << L_r << "x" << L_rC << " and " << R_r
                                      << "x" << R_rC << "!");
    return DerivativeRangeReturnType();
  }
};

template <class L_R, class R_R, size_t d, size_t L_r, size_t L_rC, typename a>
struct CombinedJac<CombinationType::difference, L_R, R_R, d, L_r, L_rC, L_r, L_rC, a>
{
  static const constexpr size_t r = CombinedDim<CombinationType::difference, L_r, L_rC, L_r, L_rC>::r;
  static const constexpr size_t rC = CombinedDim<CombinationType::difference, L_r, L_rC, L_r, L_rC>::rC;

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
  static const constexpr size_t r = CombinedDim<CombinationType::sum, L_r, L_rC, L_r, L_rC>::r;
  static const constexpr size_t rC = CombinedDim<CombinationType::sum, L_r, L_rC, L_r, L_rC>::rC;

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
template <class Left, class Right, CombinationType comb>
struct CombinedHelper
{
  static_assert(is_element_function<Left>::value || is_function<Left>::value || is_grid_function<Left>::value, "");
  static_assert(is_element_function<Right>::value || is_function<Right>::value || is_grid_function<Right>::value, "");
  static_assert(Left::d == Right::d, "");

  static const constexpr size_t d = Left::d;
  static const constexpr size_t L_r = Left::r;
  static const constexpr size_t L_rC = Left::rC;
  static const constexpr size_t R_r = Right::r;
  static const constexpr size_t R_rC = Right::rC;

  using CombinedDimHelper = CombinedDim<comb, L_r, L_rC, R_r, R_rC>;
  static const constexpr bool available = CombinedDimHelper::available;
  static const constexpr size_t r = CombinedDimHelper::r;
  static const constexpr size_t rC = CombinedDimHelper::rC;

  static int order(const Left& left, const Right& right, const Common::Parameter& param)
  {
    return CombinedOrder<comb, L_r, L_rC, R_r, R_rC>::compute(left, right, param);
  }

  using L_R = typename Left::R;
  using R_R = typename Right::R;
  using CombinedEvalHelper = CombinedEval<comb, L_R, R_R, L_r, L_rC, R_r, R_rC>;
  using R = typename CombinedEvalHelper::R;
  using RangeReturnType = typename CombinedEvalHelper::RangeReturnType;

  template <class D>
  static RangeReturnType evaluate(const Left& left, const Right& right, const D& point, const Common::Parameter& param)
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


} // namespace internal
} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_BASE_COMBINED_HH
