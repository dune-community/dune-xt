// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018 - 2020)
//   René Fritze     (2018 - 2020)
//   Tim Keil        (2018)
//   Tobias Leibner  (2018 - 2020)

#ifndef DUNE_XT_FUNCTIONS_BASE_COMBINED_FUNCTIONS_HH
#define DUNE_XT_FUNCTIONS_BASE_COMBINED_FUNCTIONS_HH

#include <dune/xt/functions/interfaces/function.hh>
#include <dune/xt/functions/type_traits.hh>

#include "combined.hh"


namespace Dune::XT::Functions {


/**
 * \brief Combined function.
 *
 *        This class combines two given functions of type LeftType and RightType
using the given combination
 *        Combination. This class (and any derived class, like Difference, Sum
or Product) can be used in two ways:
 *        - You can pass references of the left and right operand to this class.
This is done for instance when calling
 *          operator+, operator- or operator* on any function deriving from
FunctionInterface:
\code
using ConstantType = Functions::ConstantFunction< ..., double>;
ConstantType one( ... );
ConstantType two( ... );
// the following code
auto difference = one - two;
// is equivalent to
Difference< ConstantType, ConstantType > difference(one, two);
// and
internal::Combined< ConstantType, ConstantType, CombinationType::difference >
difference(one, tow);
\endcode
 *          In this situation you are responsible to ensure that the arguments
given are valid throughout the lifetime
 *          of this class. The following will lead to a segfault:
\code
using ConstantType = Functions::ConstantFunction< ..., double >;

Difference< ConstantType, ConstantType > stupid_difference()
{
  ConstantType one( ... );
  ConstantType two( ... );
  return one - two;
}
\endcode
 *        - You can pass shared_ptr of the left and right operands to this
class. In this case the following is valid:
\code
using ConstantType = Functions::ConstantFunction< ..., double >;

Difference< ConstantType, ConstantType > stupid_difference()
{
  auto one = std::make_shared< ConstantType >(1);
  auto two = std::make_shared< ConstantType >(2);
  return Difference< ConstantType, ConstantType >(one, two)
}
\endcode
 *
 * \note  Most likely you do not want to use this class diretly, but one of
Difference, Sum or Product.
 */
template <class LeftType, class RightType, class comb>
class CombinedFunction
  : public FunctionInterface<LeftType::d,
                             internal::CombinedHelper<LeftType, RightType, comb>::r,
                             internal::CombinedHelper<LeftType, RightType, comb>::rC,
                             typename internal::CombinedHelper<LeftType, RightType, comb>::R>
{
  static_assert(is_function<LeftType>::value);
  static_assert(is_function<RightType>::value);
  static_assert(LeftType::d == RightType::d);

  using BaseType = FunctionInterface<LeftType::d,
                                     internal::CombinedHelper<LeftType, RightType, comb>::r,
                                     internal::CombinedHelper<LeftType, RightType, comb>::rC,
                                     typename internal::CombinedHelper<LeftType, RightType, comb>::R>;
  using ThisType = CombinedFunction;

  using Helper = internal::CombinedHelper<LeftType, RightType, comb>;

public:
  using typename BaseType::DerivativeRangeReturnType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeReturnType;

  CombinedFunction(const LeftType& left, const RightType& right, const std::string& nm = "")
    : left_(std::move(left.copy_as_function()))
    , right_(std::move(right.copy_as_function()))
    , name_(nm.empty() ? "(" + left_->name() + " " + GetCombination<comb>::symbol() + " " + right_->name() + ")" : nm)
  {}

  CombinedFunction(const ThisType& other)
    : BaseType(other)
    , left_(other.left_->copy_as_function())
    , right_(other.right_->copy_as_function())
    , name_(other.name_)
  {}

  CombinedFunction(ThisType&& source) = default;

private:
  ThisType* copy_as_function_impl() const override
  {
    return new ThisType(*this);
  }

public:
  std::unique_ptr<ThisType> copy_as_function() const
  {
    return std::unique_ptr<ThisType>(this->copy_as_function_impl());
  }

  std::string name() const override final
  {
    return name_;
  }

  int order(const XT::Common::Parameter& param = {}) const override final
  {
    return Helper::order(*left_, *right_, param);
  }

  RangeReturnType evaluate(const DomainType& point_in_global_coordinates,
                           const Common::Parameter& param = {}) const override final
  {
    return Helper::evaluate(*left_, *right_, point_in_global_coordinates, param);
  }

  DerivativeRangeReturnType jacobian(const DomainType& point_in_global_coordinates,
                                     const Common::Parameter& param = {}) const override final
  {
    return Helper::jacobian(*left_, *right_, point_in_global_coordinates, param);
  }

private:
  std::unique_ptr<LeftType> left_;
  std::unique_ptr<RightType> right_;
  const std::string name_;
}; // class Combined


/**
 * \brief Function representing the difference between two functions.
 *
 * \see internal::Combined
 */
template <class MinuendType, class SubtrahendType>
class DifferenceFunction : public CombinedFunction<MinuendType, SubtrahendType, CombinationType::difference>
{
  using BaseType = CombinedFunction<MinuendType, SubtrahendType, CombinationType::difference>;

public:
  template <class... Args>
  explicit DifferenceFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
}; // class DifferenceFunction


/**
 * \brief Function representing the sum of two functions.
 *
 * \see internal::Combined
 */
template <class LeftSummandType, class RightSummandType>
class SumFunction : public CombinedFunction<LeftSummandType, RightSummandType, CombinationType::sum>
{
  using BaseType = CombinedFunction<LeftSummandType, RightSummandType, CombinationType::sum>;

public:
  template <class... Args>
  explicit SumFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
}; // class SumFunction


/**
 * \brief Function representing the fraction of two functions.
 *
 * \see CombinedFunction
 */
template <class NominatorType, class DenominatorType>
class FractionFunction : public CombinedFunction<NominatorType, DenominatorType, CombinationType::fraction>
{
  using BaseType = CombinedFunction<NominatorType, DenominatorType, CombinationType::fraction>;

public:
  template <class... Args>
  explicit FractionFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
}; // class FractionFunction


/**
 * \brief Function representing the product of two functions.
 *
 * \see internal::Combined
 */
template <class LeftFactorType, class RightFactorType>
class ProductFunction : public CombinedFunction<LeftFactorType, RightFactorType, CombinationType::product>
{
  using BaseType = CombinedFunction<LeftFactorType, RightFactorType, CombinationType::product>;

public:
  template <class... Args>
  explicit ProductFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
}; // class ProductFunction


} // namespace Dune::XT::Functions

#endif // DUNE_XT_FUNCTIONS_BASE_COMBINED_HH
