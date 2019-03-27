// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_XT_FUNCTIONS_BASE_COMBINED_ELEMENT_FUNCTIONS_HH
#define DUNE_XT_FUNCTIONS_BASE_COMBINED_ELEMENT_FUNCTIONS_HH

#include <dune/xt/common/memory.hh>

#include <dune/xt/functions/interfaces/element-functions.hh>
#include <dune/xt/functions/type_traits.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace internal {


/**
 * \brief Helper class defining types of combined functions, if available.
 *
 * \note Most likely you do not want to use this class directly, but CombinedConstElementFunction or
 *       CombinedElementFunction.
 */
template <class LeftType, class RightType, CombinationType comb>
class CombinedElementFunctionHelper
{
  static_assert(is_element_function<LeftType>::value, "");
  static_assert(is_element_function<RightType>::value, "");

public:
  using E = typename LeftType::E;
  using R = typename LeftType::R;

private:
  using D = typename LeftType::D;
  static const constexpr size_t d = LeftType::d;

private:
  static_assert(std::is_same<typename RightType::E, E>::value, "");
  static_assert(std::is_same<typename RightType::D, D>::value, "");
  static_assert(RightType::d == d, "");
  static_assert(std::is_same<typename RightType::R, R>::value, "");

  template <class L, class R>
  class dim_switch
  {
    template <CombinationType cc = comb,
              size_t rL = L::r,
              size_t rCL = L::rC,
              size_t rR = R::r,
              size_t rCR = R::rC,
              bool anything = true>
    class Dimension
    {
    public:
      static const bool available = false;
    };

    template <size_t r_, size_t rC_, bool anything>
    class Dimension<CombinationType::difference, r_, rC_, r_, rC_, anything>
    {
    public:
      static const bool available = true;
      static const size_t r = r_;
      static const size_t rC = rC_;
    };

    template <size_t r_, size_t rC_, bool anything>
    class Dimension<CombinationType::sum, r_, rC_, r_, rC_, anything>
    {
    public:
      static const bool available = true;
      static const size_t r = r_;
      static const size_t rC = rC_;
    };

    template <bool anything>
    class Dimension<CombinationType::product, 1, 1, 1, 1, anything>
    {
    public:
      static const bool available = true;
      static const size_t r = 1;
      static const size_t rC = 1;
    };

    template <size_t r_, size_t rC_, bool anything>
    class Dimension<CombinationType::product, r_, rC_, 1, 1, anything>
    {
    public:
      static const bool available = true;
      static const size_t r = r_;
      static const size_t rC = rC_;
    };

    template <size_t r_, size_t rC_, bool anything>
    class Dimension<CombinationType::product, 1, 1, r_, rC_, anything>
    {
    public:
      static const bool available = true;
      static const size_t r = r_;
      static const size_t rC = rC_;
    };

    template <size_t rL, size_t c, size_t rR, bool anything>
    class Dimension<CombinationType::product, rL, c, c, rR, anything>
    {
    public:
      static const bool available = true;
      static const size_t r = rL;
      static const size_t rC = rR;
    };

  public:
    static const bool available = Dimension<>::available;
    static const size_t r = Dimension<>::r;
    static const size_t rC = Dimension<>::rC;
  }; // class dim_switch

public:
  static const size_t r = dim_switch<LeftType, RightType>::r;
  static const size_t rC = dim_switch<LeftType, RightType>::rC;
  static const bool available = dim_switch<LeftType, RightType>::available;

private:
  using DomainType = typename ElementFunctionInterface<E, r, rC, R>::DomainType;
  using RangeReturnType = typename RangeTypeSelector<R, r, rC>::return_type;
  using DerivativeRangeReturnType = typename DerivativeRangeTypeSelector<d, R, r, rC>::return_type;

  template <CombinationType cc = comb,
            size_t rL = LeftType::r,
            size_t rCL = LeftType::rC,
            size_t rR = RightType::r,
            size_t rCR = RightType::rC,
            bool anything = true>
  class CombinationTypeSwitch
  {
    static_assert(!anything, "Nothing available for these CombinationType!");
  };

  template <size_t r_, size_t rC_, bool anything>
  class CombinationTypeSwitch<CombinationType::difference, r_, rC_, r_, rC_, anything>
  {
  public:
    static int order(const int left_order, const int right_order)
    {
      return std::max(left_order, right_order);
    }

    static RangeReturnType evaluate(const LeftType& left_local,
                                    const RightType& right_local,
                                    const DomainType& point_in_reference_element,
                                    const Common::Parameter& param)
    {
      return left_local.evaluate(point_in_reference_element, param)
             - right_local.evaluate(point_in_reference_element, param);
    }

    static DerivativeRangeReturnType jacobian(const LeftType& left_local,
                                              const RightType& right_local,
                                              const DomainType& point_in_reference_element,
                                              const Common::Parameter& param)
    {
      return left_local.jacobian(point_in_reference_element, param)
             - right_local.jacobian(point_in_reference_element, param);
    }
  }; // class CombinationTypeSwitch< ..., difference >

  template <size_t r_, size_t rC_, bool anything>
  class CombinationTypeSwitch<CombinationType::sum, r_, rC_, r_, rC_, anything>
  {
  public:
    static int order(const int left_order, const int right_order)
    {
      return std::max(left_order, right_order);
    }

    static RangeReturnType evaluate(const LeftType& left_local,
                                    const RightType& right_local,
                                    const DomainType& point_in_reference_element,
                                    const Common::Parameter& param)
    {
      return left_local.evaluate(point_in_reference_element, param)
             + right_local.evaluate(point_in_reference_element, param);
    }

    static DerivativeRangeReturnType jacobian(const LeftType& left_local,
                                              const RightType& right_local,
                                              const DomainType& point_in_reference_element,
                                              const Common::Parameter& param)
    {
      return left_local.jacobian(point_in_reference_element, param)
             + right_local.jacobian(point_in_reference_element, param);
    }
  }; // class CombinationTypeSwitch< ..., sum >

  template <bool anything>
  class CombinationTypeSwitch<CombinationType::product, 1, 1, 1, 1, anything>
  {
  public:
    static int order(const int left_order, const int right_order)
    {
      return left_order + right_order;
    }

    static RangeReturnType evaluate(const LeftType& left_local,
                                    const RightType& right_local,
                                    const DomainType& point_in_reference_element,
                                    const Common::Parameter& param)
    {
      return left_local.evaluate(point_in_reference_element, param)
             * right_local.evaluate(point_in_reference_element, param);
    }

    static DerivativeRangeReturnType jacobian(const LeftType& /*left_local*/,
                                              const RightType& /*right_local*/,
                                              const DomainType& /*point_in_reference_element*/,
                                              const Common::Parameter& /*param*/)
    {
      DUNE_THROW(NotImplemented, "If you need this, implement it!");
      return DerivativeRangeReturnType();
    }
  }; // class CombinationTypeSwitch< ..., product >

  template <size_t r_, size_t rC_, bool anything>
  class CombinationTypeSwitch<CombinationType::product, r_, rC_, 1, 1, anything>
  {
  public:
    static int order(const int left_order, const int right_order)
    {
      return left_order + right_order;
    }

    static RangeReturnType evaluate(const LeftType& left_local,
                                    const RightType& right_local,
                                    const DomainType& point_in_reference_element,
                                    const Common::Parameter& param)
    {
      auto result = left_local.evaluate(point_in_reference_element, param);
      result *= right_local.evaluate(point_in_reference_element, param);
      return result;
    }

    static DerivativeRangeReturnType jacobian(const LeftType& /*left_local*/,
                                              const RightType& /*right_local*/,
                                              const DomainType& /*point_in_reference_element*/,
                                              const Common::Parameter& /*param*/)
    {
      DUNE_THROW(NotImplemented, "If you need this, implement it!");
      return DerivativeRangeReturnType();
    }
  }; // class CombinationTypeSwitch< ..., product >

  template <size_t r_, size_t rC_, bool anything>
  class CombinationTypeSwitch<CombinationType::product, 1, 1, r_, rC_, anything>
  {
  public:
    static int order(const int left_order, const int right_order)
    {
      return left_order + right_order;
    }

    static RangeReturnType evaluate(const LeftType& left_local,
                                    const RightType& right_local,
                                    const DomainType& point_in_reference_element,
                                    const Common::Parameter& param)
    {
      auto result = right_local.evaluate(point_in_reference_element, param);
      result *= left_local.evaluate(point_in_reference_element, param);
      return result;
    }

    static DerivativeRangeReturnType jacobian(const LeftType& /*left_local*/,
                                              const RightType& /*right_local*/,
                                              const DomainType& /*point_in_reference_element*/,
                                              const Common::Parameter& /*param*/)
    {
      DUNE_THROW(NotImplemented, "If you need this, implement it!");
      return DerivativeRangeReturnType();
    }
  }; // class CombinationTypeSwitch< ..., product >

  template <size_t rL, size_t c_, size_t rR_, bool anything>
  class CombinationTypeSwitch<CombinationType::product, rL, c_, c_, rR_, anything>
  {
  public:
    static int order(const int left_order, const int right_order)
    {
      return left_order + right_order;
    }

    static RangeReturnType evaluate(const LeftType& left_local,
                                    const RightType& right_local,
                                    const DomainType& point_in_reference_element,
                                    const Common::Parameter& param)
    {
      return left_local.evaluate(point_in_reference_element, param)
             * right_local.evaluate(point_in_reference_element, param);
    }

    static DerivativeRangeReturnType jacobian(const LeftType& /*left*/,
                                              const RightType& /*right*/,
                                              const DomainType& /*xx*/,
                                              const Common::Parameter& /*param*/)
    {
      DUNE_THROW(NotImplemented, "If you need this, implement it!");
      return DerivativeRangeReturnType();
    }
  }; // class CombinationTypeSwitch< ..., product >

public:
  static int order(const LeftType& left, const RightType& right, const Common::Parameter& param)
  {
    return CombinationTypeSwitch<>::order(left.order(param), right.order(param));
  }

  static RangeReturnType
  evaluate(const LeftType& left, const RightType& right, const DomainType& xx, const Common::Parameter& param)
  {
    return CombinationTypeSwitch<>::evaluate(left, right, xx, param);
  }

  static DerivativeRangeReturnType
  jacobian(const LeftType& left, const RightType& right, const DomainType& xx, const Common::Parameter& param)
  {
    return CombinationTypeSwitch<>::jacobian(left, right, xx, param);
  }
}; // class CombinedElementFunctionHelper


template <class LeftType, class RightType>
struct DualStorageProvider
{
  XT::Common::StorageProvider<LeftType> left;
  XT::Common::StorageProvider<RightType> right;

  DualStorageProvider(LeftType& lft, RightType& rght)
    : left(lft)
    , right(rght)
  {}

  DualStorageProvider(std::shared_ptr<LeftType> lft, std::shared_ptr<RightType> rght)
    : left(lft)
    , right(rght)
  {}

  DualStorageProvider(std::unique_ptr<LeftType>&& lft, std::unique_ptr<RightType>&& rght)
    : left(std::move(lft))
    , right(std::move(rght))
  {}
}; // struct DualStorageProvider


} // namespace internal


template <class LeftType, class RightType, CombinationType combination>
class CombinedConstElementFunction
  : public ElementFunctionInterface<
        typename internal::CombinedElementFunctionHelper<LeftType, RightType, combination>::E,
        internal::CombinedElementFunctionHelper<LeftType, RightType, combination>::r,
        internal::CombinedElementFunctionHelper<LeftType, RightType, combination>::rC,
        typename internal::CombinedElementFunctionHelper<LeftType, RightType, combination>::R>
{
  using BaseType =
      ElementFunctionInterface<typename internal::CombinedElementFunctionHelper<LeftType, RightType, combination>::E,
                               internal::CombinedElementFunctionHelper<LeftType, RightType, combination>::r,
                               internal::CombinedElementFunctionHelper<LeftType, RightType, combination>::rC,
                               typename internal::CombinedElementFunctionHelper<LeftType, RightType, combination>::R>;

  using Select = internal::CombinedElementFunctionHelper<LeftType, RightType, combination>;

public:
  using typename BaseType::DerivativeRangeReturnType;
  using typename BaseType::DerivativeRangeType;
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::RangeReturnType;
  using typename BaseType::RangeType;

  CombinedConstElementFunction(const LeftType& left, const RightType& right)
    : BaseType()
    , left_(left)
    , right_(right)
    , bind_is_temporarily_ok_(false)
  {
    bind_if_arguments_were_bound();
  }

  CombinedConstElementFunction(std::shared_ptr<const LeftType> left, std::shared_ptr<const RightType> right)
    : BaseType()
    , left_(left)
    , right_(right)
    , bind_is_temporarily_ok_(false)
  {
    bind_if_arguments_were_bound();
  }

  CombinedConstElementFunction(std::unique_ptr<const LeftType>&& left, std::unique_ptr<const RightType>&& right)
    : BaseType()
    , left_(std::move(left))
    , right_(std::move(right))
    , bind_is_temporarily_ok_(false)
  {
    bind_if_arguments_were_bound();
  }

protected:
  void post_bind(const ElementType& /*element*/) override
  {
    DUNE_THROW_IF(!bind_is_temporarily_ok_, XT::Common::Exceptions::you_are_using_this_wrong, "");
  }

public:
  int order(const XT::Common::Parameter& param = {}) const override final
  {
    return Select::order(left_.access(), right_.access(), param);
  }

  RangeReturnType evaluate(const DomainType& point_in_reference_element,
                           const Common::Parameter& param = {}) const override final
  {
    return Select::evaluate(left_.access(), right_.access(), point_in_reference_element, param);
  }

  DerivativeRangeReturnType jacobian(const DomainType& point_in_reference_element,
                                     const Common::Parameter& param = {}) const override final
  {
    return Select::jacobian(left_.access(), right_.access(), point_in_reference_element, param);
  }

private:
  void bind_if_arguments_were_bound()
  {
    bool left_bound = false;
    bool right_bound = false;
    try {
      left_.access().element();
      left_bound = true;
    } catch (const XT::Grid::Exceptions::not_bound_to_an_element_yet&) {
      left_bound = false;
    }
    try {
      right_.access().element();
      right_bound = true;
    } catch (const XT::Grid::Exceptions::not_bound_to_an_element_yet&) {
      right_bound = false;
    }
    if (left_bound && right_bound) {
      DUNE_THROW_IF(left_.access().element() != right_.access().element(), Exceptions::wrong_input_given, "");
      bind_is_temporarily_ok_ = true;
      this->bind(left_.access().element());
      bind_is_temporarily_ok_ = false;
    }
  } // ... bind_if_arguments_were_bound(...)

  Common::ConstStorageProvider<LeftType> left_;
  Common::ConstStorageProvider<RightType> right_;
  bool bind_is_temporarily_ok_;
}; // class CombinedConstElementFunction


template <class LeftType, class RightType, CombinationType combination>
class CombinedElementFunction
  : internal::DualStorageProvider<LeftType, RightType>
  , public CombinedConstElementFunction<LeftType, RightType, combination>
{
  using Storage = internal::DualStorageProvider<LeftType, RightType>;
  using BaseType = CombinedConstElementFunction<LeftType, RightType, combination>;

public:
  using typename BaseType::ElementType;

  CombinedElementFunction(LeftType& left, RightType& right)
    : Storage(left, right)
    , BaseType(Storage::left.access(), Storage::right.access())
  {}

  CombinedElementFunction(std::shared_ptr<LeftType> left, std::shared_ptr<RightType> right)
    : Storage(left, right)
    , BaseType(Storage::left.access(), Storage::right.access())
  {}

  CombinedElementFunction(std::unique_ptr<LeftType>&& left, std::unique_ptr<RightType>&& right)
    : Storage(std::move(left), std::move(right))
    , BaseType(Storage::left.access(), Storage::right.access())
  {}

protected:
  void post_bind(const ElementType& element) override final
  {
    Storage::left.access().bind(element);
    Storage::right.access().bind(element);
  }
}; // class CombinedElementFunction


/**
 * \brief Element function representing the difference between two element functions.
 *
 * \sa CombinedConstElementFunction
 */
template <class MinuendType, class SubtrahendType>
class ConstDifferenceElementFunction
  : public CombinedConstElementFunction<MinuendType, SubtrahendType, CombinationType::difference>
{
  using BaseType = CombinedConstElementFunction<MinuendType, SubtrahendType, CombinationType::difference>;

public:
  template <class... Args>
  explicit ConstDifferenceElementFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
};


/**
 * \brief Element function representing the difference between two element functions.
 *
 * \sa CombinedElementFunction
 */
template <class MinuendType, class SubtrahendType>
class DifferenceElementFunction
  : public CombinedElementFunction<MinuendType, SubtrahendType, CombinationType::difference>
{
  using BaseType = CombinedElementFunction<MinuendType, SubtrahendType, CombinationType::difference>;

public:
  template <class... Args>
  explicit DifferenceElementFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
};


/**
 * \brief Element function representing the sum of two element functions.
 *
 * \sa CombinedConstElementFunction
 */
template <class LeftSummandType, class RightSummandType>
class ConstSumElementFunction
  : public CombinedConstElementFunction<LeftSummandType, RightSummandType, CombinationType::sum>
{
  using BaseType = CombinedConstElementFunction<LeftSummandType, RightSummandType, CombinationType::sum>;

public:
  template <class... Args>
  explicit ConstSumElementFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
};


/**
 * \brief Element function representing the sum of two element functions.
 *
 * \sa CombinedElementFunction
 */
template <class LeftSummandType, class RightSummandType>
class SumElementFunction : public CombinedElementFunction<LeftSummandType, RightSummandType, CombinationType::sum>
{
  using BaseType = CombinedElementFunction<LeftSummandType, RightSummandType, CombinationType::sum>;

public:
  template <class... Args>
  explicit SumElementFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
};


/**
 * \brief Element function representing the product of two element functions.
 *
 * \sa CombinedElementFunction
 */
template <class LeftFactorType, class RightFactorType>
class ProductElementFunction : public CombinedElementFunction<LeftFactorType, RightFactorType, CombinationType::product>
{
  using BaseType = CombinedElementFunction<LeftFactorType, RightFactorType, CombinationType::product>;

public:
  template <class... Args>
  explicit ProductElementFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
};


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_BASE_COMBINED_HH
