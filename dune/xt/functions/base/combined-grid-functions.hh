// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018 - 2019)
//   René Fritze     (2018 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2018 - 2020)

#ifndef DUNE_XT_FUNCTIONS_BASE_COMBINED_GRID_FUNCTIONS_HH
#define DUNE_XT_FUNCTIONS_BASE_COMBINED_GRID_FUNCTIONS_HH

#include <dune/xt/functions/interfaces/grid-function.hh>

#include "combined.hh"

namespace Dune::XT::Functions {


/**
 * \brief Base combined grid function.
 *
\code
using IndicatorType = Functions::IndicatorFunction< ..., double>;
IndicatorType one( ... );
IndicatorType two( ... );
// the following code
auto difference = one - two;
// is equivalent to
Difference< IndicatorType, IndicatorType > difference(one, two);
// and
internal::Combined< IndicatorType, IndicatorType, CombinationType::difference > difference(one, tow);
\endcode
 *          In this situation you are responsible to ensure that the arguments given are valid throughout the lifetime
 *          of this class. The following will lead to a segfault:
\code
using IndicatorType = Functions::IndicatorFunction< ..., double >;

Difference< IndicatorType, IndicatorType > stupid_difference()
{
  IndicatorType one( ... );
  IndicatorType two( ... );
  return one - two;
}
\endcode
 *        - You can pass shared_ptr of the left and right operands to this class. In this case the following is valid:
\code
using IndicatorType = Functions::IndicatorFunction< ..., double >;

Difference<IndicatorType, IndicatorType> stupid_difference()
{
  auto one = std::make_shared<IndicatorType>(1);
  auto two = std::make_shared<IndicatorType>(2);
  return Difference<IndicatorType, IndicatorType>(one, two)
}
\endcode
 *
 * \note  Most likely you do not want to use this class diretly, but one of Difference, Fraction, Sum or Product.
 *
 * \todo Implement custom local function to hold a copy of this!
 */
template <class LeftType, class RightType, class comb>
class CombinedGridFunction
  : public GridFunctionInterface<typename LeftType::E,
                                 internal::CombinedHelper<LeftType, RightType, comb>::r,
                                 internal::CombinedHelper<LeftType, RightType, comb>::rC,
                                 typename internal::CombinedHelper<LeftType, RightType, comb>::R>
{
  static_assert(is_grid_function<LeftType>::value);
  static_assert(is_grid_function<RightType>::value);
  static_assert(std::is_same<typename LeftType::E, typename RightType::E>::value);

  using ThisType = CombinedGridFunction;
  using BaseType = GridFunctionInterface<typename LeftType::E,
                                         internal::CombinedHelper<LeftType, RightType, comb>::r,
                                         internal::CombinedHelper<LeftType, RightType, comb>::rC,
                                         typename internal::CombinedHelper<LeftType, RightType, comb>::R>;

public:
  using ElementType = typename BaseType::ElementType;
  using LocalFunctionType = typename BaseType::LocalFunctionType;

  CombinedGridFunction(const LeftType& left,
                       const RightType& right,
                       const std::string nm = "",
                       const std::string& logging_prefix = "")
    : BaseType(left.parameter_type() + right.parameter_type(),
               logging_prefix.empty() ? Common::to_camel_case(get_combination_name(comb{}) + "GridFunction")
                                      : logging_prefix,
               logging_prefix.empty())
    , left_(left.copy_as_grid_function())
    , right_(right.copy_as_grid_function())
    , name_(nm.empty() ? "(" + left_->name() + GetCombination<comb>::symbol() + right_->name() + ")" : nm)
  {
    LOG_(debug) << Common::to_camel_case(get_combination_name(comb{}) + "GridFunction") << "(left=" << &left
                << ", right=" << &right << ", nm=\"" << nm << "\")" << std::endl;
  }

  CombinedGridFunction(LeftType*&& left,
                       RightType*&& right,
                       const std::string nm = "",
                       const std::string& logging_prefix = "")
    : BaseType(left->parameter_type() + right->parameter_type(),
               logging_prefix.empty() ? Common::to_camel_case(get_combination_name(comb{}) + "GridFunction")
                                      : logging_prefix,
               logging_prefix.empty())
    , left_(std::move(left))
    , right_(std::move(right))
    , name_(nm.empty() ? "(" + left_->name() + GetCombination<comb>::symbol() + right_->name() + ")" : nm)
  {
    LOG_(debug) << Common::to_camel_case(get_combination_name(comb{}) + "GridFunction") << "(left=" << left
                << ", right=" << right << ", nm=\"" << nm << "\")" << std::endl;
  }

  CombinedGridFunction(const ThisType& other)
    : BaseType(other)
    , left_(other.left_->copy_as_grid_function())
    , right_(other.right_->copy_as_grid_function())
    , name_(other.name_)
  {}

  CombinedGridFunction(ThisType&& source) noexcept = default;

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    LOG_(debug) << Common::to_camel_case(get_combination_name(comb{}) + "GridFunction") + "::local_function()"
                << std::endl;
    using LeftLF = typename LeftType::LocalFunctionType;
    using RightLF = typename RightType::LocalFunctionType;
    return std::make_unique<CombinedElementFunction<LeftLF, RightLF, comb>>(std::move(left_->local_function()),
                                                                            std::move(right_->local_function()));
  } // ... local_function(...)


private:
  ThisType* copy_as_grid_function_impl() const override
  {
    return new ThisType(*this);
  }

public:
  std::unique_ptr<ThisType> copy_as_grid_function() const
  {
    return std::unique_ptr<ThisType>(this->copy_as_grid_function_impl());
  }

  std::string name() const override final
  {
    return name_;
  }

private:
  std::unique_ptr<GridFunctionInterface<typename LeftType::E, LeftType::r, LeftType::rC, typename LeftType::R>> left_;
  std::unique_ptr<GridFunctionInterface<typename RightType::E, RightType::r, RightType::rC, typename RightType::R>>
      right_;
  const std::string name_;
}; // class CombinedGridFunction


/**
 * \brief Function representing the difference between two functions.
 *
 * \see CombinedGridFunction
 */
template <class MinuendType, class SubtrahendType>
class DifferenceGridFunction : public CombinedGridFunction<MinuendType, SubtrahendType, CombinationType::difference>
{
  using BaseType = CombinedGridFunction<MinuendType, SubtrahendType, CombinationType::difference>;

public:
  template <class... Args>
  explicit DifferenceGridFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
}; // class DifferenceGridFunction


/**
 * \brief Function representing the sum of two functions.
 *
 * \see CombinedGridFunction
 */
template <class LeftSummandType, class RightSummandType>
class SumGridFunction : public CombinedGridFunction<LeftSummandType, RightSummandType, CombinationType::sum>
{
  using BaseType = CombinedGridFunction<LeftSummandType, RightSummandType, CombinationType::sum>;

public:
  template <class... Args>
  explicit SumGridFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
}; // class SumGridFunction


/**
 * \brief Function representing the fraction of two functions.
 *
 * \see CombinedGridFunction
 */
template <class NominatorType, class DenominatorType>
class FractionGridFunction : public CombinedGridFunction<NominatorType, DenominatorType, CombinationType::fraction>
{
  using BaseType = CombinedGridFunction<NominatorType, DenominatorType, CombinationType::fraction>;

public:
  template <class... Args>
  explicit FractionGridFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
}; // class FractionGridFunction


/**
 * \brief Grid function representing the product of two grid functions.
 *
 * \see CombinedGridFunction
 */
template <class LeftFactorType, class RightFactorType>
class ProductGridFunction : public CombinedGridFunction<LeftFactorType, RightFactorType, CombinationType::product>
{
  using BaseType = CombinedGridFunction<LeftFactorType, RightFactorType, CombinationType::product>;

public:
  template <class... Args>
  explicit ProductGridFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
}; // class ProductGridFunction


} // namespace Dune::XT::Functions

#endif // DUNE_XT_FUNCTIONS_BASE_COMBINED_HH
