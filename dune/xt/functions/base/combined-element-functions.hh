// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2019)
//   Tobias Leibner  (2019 - 2020)

#ifndef DUNE_XT_FUNCTIONS_BASE_COMBINED_ELEMENT_FUNCTIONS_HH
#define DUNE_XT_FUNCTIONS_BASE_COMBINED_ELEMENT_FUNCTIONS_HH

#include <dune/xt/common/type_traits.hh>

#include <dune/xt/functions/exceptions.hh>
#include <dune/xt/functions/interfaces/element-functions.hh>
#include <dune/xt/functions/type_traits.hh>

#include "combined.hh"

namespace Dune::XT::Functions {


template <class LeftType, class RightType, typename comb>
class CombinedConstElementFunction
  : public ElementFunctionInterface<typename LeftType::E,
                                    internal::CombinedHelper<LeftType, RightType, comb>::r,
                                    internal::CombinedHelper<LeftType, RightType, comb>::rC,
                                    typename internal::CombinedHelper<LeftType, RightType, comb>::R>
{
  static_assert(is_element_function<LeftType>::value);
  static_assert(is_element_function<RightType>::value);
  static_assert(std::is_same<typename LeftType::E, typename RightType::E>::value);

  using ThisType = CombinedConstElementFunction;
  using BaseType = ElementFunctionInterface<typename LeftType::E,
                                            internal::CombinedHelper<LeftType, RightType, comb>::r,
                                            internal::CombinedHelper<LeftType, RightType, comb>::rC,
                                            typename internal::CombinedHelper<LeftType, RightType, comb>::R>;
  using Helper = internal::CombinedHelper<LeftType, RightType, comb>;

public:
  using typename BaseType::DerivativeRangeReturnType;
  using typename BaseType::DerivativeRangeType;
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::RangeReturnType;
  using typename BaseType::RangeType;

  CombinedConstElementFunction(const LeftType& left, const RightType& right)
    : BaseType(left.parameter_type() + right.parameter_type())
    , left_(left)
    , right_(right)
    , bind_is_temporarily_ok_(false)
  {
    bind_if_arguments_were_bound();
  }

  CombinedConstElementFunction(const LeftType& left, RightType&& right)
    : BaseType(left.parameter_type() + right.parameter_type())
    , left_(left)
    , right_(std::move(right))
    , bind_is_temporarily_ok_(false)
  {
    bind_if_arguments_were_bound();
  }

  CombinedConstElementFunction(std::shared_ptr<const LeftType> left, std::shared_ptr<const RightType> right)
    : BaseType(left->parameter_type() + right->parameter_type())
    , left_(left)
    , right_(right)
    , bind_is_temporarily_ok_(false)
  {
    bind_if_arguments_were_bound();
  }

  CombinedConstElementFunction(std::unique_ptr<const LeftType>&& left, std::unique_ptr<const RightType>&& right)
    : BaseType(left->parameter_type() + right->parameter_type())
    , left_(std::move(left))
    , right_(std::move(right))
    , bind_is_temporarily_ok_(false)
  {
    bind_if_arguments_were_bound();
  }

  CombinedConstElementFunction(const ThisType&) = default;

  CombinedConstElementFunction(ThisType&&) = default;

protected:
  void post_bind(const ElementType& /*element*/) override
  {
    DUNE_THROW_IF(!bind_is_temporarily_ok_, XT::Common::Exceptions::you_are_using_this_wrong, "");
  }

public:
  int order(const XT::Common::Parameter& param = {}) const override final
  {
    return Helper::order(left_.access(), right_.access(), param);
  }

  RangeReturnType evaluate(const DomainType& point_in_reference_element,
                           const Common::Parameter& param = {}) const override final
  {
    return Helper::evaluate(left_.access(), right_.access(), point_in_reference_element, param);
  }

  DerivativeRangeReturnType jacobian(const DomainType& point_in_reference_element,
                                     const Common::Parameter& param = {}) const override final
  {
    return Helper::jacobian(left_.access(), right_.access(), point_in_reference_element, param);
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


template <class LeftType, class RightType, class combination>
class CombinedElementFunction
  : internal::CombinedStorageProvider<LeftType, RightType>
  , public CombinedConstElementFunction<LeftType, RightType, combination>
{
  using ThisType = CombinedElementFunction;
  using Storage = internal::CombinedStorageProvider<LeftType, RightType>;
  using BaseType = CombinedConstElementFunction<LeftType, RightType, combination>;

public:
  using typename BaseType::ElementType;

  CombinedElementFunction(LeftType& lft, RightType& rght)
    : Storage(lft, rght)
    , BaseType(Storage::left.access(), Storage::right.access())
  {}

  CombinedElementFunction(LeftType& lft, RightType&& rght)
    : Storage(lft, std::move(rght))
    , BaseType(Storage::left.access(), Storage::right.access())
  {}

  CombinedElementFunction(std::shared_ptr<LeftType> lft, std::shared_ptr<RightType> rght)
    : Storage(lft, rght)
    , BaseType(Storage::left.access(), Storage::right.access())
  {}

  CombinedElementFunction(std::unique_ptr<LeftType>&& lft, std::unique_ptr<RightType>&& rght)
    : Storage(std::move(lft), std::move(rght))
    , BaseType(Storage::left.access(), Storage::right.access())
  {}

  CombinedElementFunction(const ThisType&) = default;

  CombinedElementFunction(ThisType&&) = default;

protected:
  void post_bind(const ElementType& element) override final
  {
    Storage::left.access().bind(element);
    Storage::right.access().bind(element);
  }
}; // class CombinedElementFunction


/**
 * \brief Element function representing the difference between two element functions (const variant).
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
 * \brief Element function representing the sum of two element functions (const variant).
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
 * \brief Element function representing the product of two element functions (const variant).
 *
 * \sa CombinedElementFunction
 */
template <class LeftFactorType, class RightFactorType>
class ConstProductElementFunction
  : public CombinedConstElementFunction<LeftFactorType, RightFactorType, CombinationType::product>
{
  using BaseType = CombinedConstElementFunction<LeftFactorType, RightFactorType, CombinationType::product>;

public:
  template <class... Args>
  explicit ConstProductElementFunction(Args&&... args)
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


/**
 * \brief Element function representing the fraction of two element functions (const variant).
 *
 * \sa CombinedElementFunction
 */
template <class LeftFactorType, class RightFactorType>
class ConstFractionElementFunction
  : public CombinedConstElementFunction<LeftFactorType, RightFactorType, CombinationType::fraction>
{
  using BaseType = CombinedConstElementFunction<LeftFactorType, RightFactorType, CombinationType::fraction>;

public:
  template <class... Args>
  explicit ConstFractionElementFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
};


/**
 * \brief Element function representing the fraction of two element functions.
 *
 * \sa CombinedElementFunction
 */
template <class LeftFactorType, class RightFactorType>
class FractionElementFunction
  : public CombinedElementFunction<LeftFactorType, RightFactorType, CombinationType::fraction>
{
  using BaseType = CombinedElementFunction<LeftFactorType, RightFactorType, CombinationType::fraction>;

public:
  template <class... Args>
  explicit FractionElementFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
};


} // namespace Dune::XT::Functions

#endif // DUNE_XT_FUNCTIONS_BASE_COMBINED_HH
