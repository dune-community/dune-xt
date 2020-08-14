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
#include <dune/xt/functions/type_traits.hh>

#include "combined-element-functions.hh"

namespace Dune {
namespace XT {
namespace Functions {


/**
 * \brief Base combined grid function.
 *
 *        This class combines two given grid functions of type LeftType and RightType using the given combination
 *        Combination. This class (and any derived class, like Difference, Sum or Product) can be used in two ways:
 *        - You can pass references of the left and right operand to this class. This is done for instance when calling
 *          operator+, operator- or operator* on any function deriving from GridFunctionInterface:
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
 */
template <class LeftType, class RightType, CombinationType combination>
class CombinedGridFunction
  : public GridFunctionInterface<typename internal::CombinedElementFunctionHelper<LeftType, RightType, combination>::E,
                                 internal::CombinedElementFunctionHelper<LeftType, RightType, combination>::r,
                                 internal::CombinedElementFunctionHelper<LeftType, RightType, combination>::rC,
                                 typename internal::CombinedElementFunctionHelper<LeftType, RightType, combination>::R>
{
  using ThisType = CombinedGridFunction;
  using BaseType =
      GridFunctionInterface<typename internal::CombinedElementFunctionHelper<LeftType, RightType, combination>::E,
                            internal::CombinedElementFunctionHelper<LeftType, RightType, combination>::r,
                            internal::CombinedElementFunctionHelper<LeftType, RightType, combination>::rC,
                            typename internal::CombinedElementFunctionHelper<LeftType, RightType, combination>::R>;

public:
  using ElementType = typename BaseType::ElementType;
  using LocalFunctionType = typename BaseType::LocalFunctionType;

  CombinedGridFunction(const LeftType& left,
                       const RightType& right,
                       const std::string nm = "",
                       const std::string& logging_prefix = "")
    : BaseType(left.parameter_type() + right.parameter_type(),
               logging_prefix.empty() ? Common::to_camel_case(GetCombination<combination>::name() + "GridFunction")
                                      : logging_prefix,
               logging_prefix.empty())
    , left_(left)
    , right_(right)
    , name_(nm.empty()
                ? "(" + left_.access().name() + GetCombination<combination>::symbol() + right_.access().name() + ")"
                : nm)
  {
    LOG_(debug) << Common::to_camel_case(GetCombination<combination>::name() + "GridFunction") << "(left=" << &left
                << ", right=" << &right << ", nm=\"" << nm << "\")" << std::endl;
  }

  CombinedGridFunction(const LeftType& left,
                       RightType&& right,
                       const std::string nm = "",
                       const std::string& logging_prefix = "")
    : BaseType(left.parameter_type() + right.parameter_type(),
               logging_prefix.empty() ? Common::to_camel_case(GetCombination<combination>::name() + "GridFunction")
                                      : logging_prefix,
               logging_prefix.empty())
    , left_(left)
    , right_(std::move(right))
    , name_(nm.empty()
                ? "(" + left_.access().name() + GetCombination<combination>::symbol() + right_.access().name() + ")"
                : nm)
  {
    LOG_(debug) << Common::to_camel_case(GetCombination<combination>::name() + "GridFunction") << "(left=" << &left
                << ", right=" << &right << ", nm=\"" << nm << "\")" << std::endl;
  }

  CombinedGridFunction(std::shared_ptr<const LeftType> left,
                       std::shared_ptr<const RightType> right,
                       const std::string nm = "",
                       const std::string& logging_prefix = "")
    : BaseType(left->parameter_type() + right->parameter_type(),
               logging_prefix.empty() ? Common::to_camel_case(GetCombination<combination>::name() + "GridFunction")
                                      : logging_prefix,
               logging_prefix.empty())
    , left_(left)
    , right_(right)
    , name_(nm.empty()
                ? "(" + left_.access().name() + GetCombination<combination>::symbol() + right_.access().name() + ")"
                : nm)
  {
    LOG_(debug) << Common::to_camel_case(GetCombination<combination>::name() + "GridFunction")
                << "(left_shrd_ptr=" << left << ", right_shrd_ptr=" << right << ", nm=\"" << nm << "\")" << std::endl;
  }

  CombinedGridFunction(LeftType*&& left,
                       RightType*&& right,
                       const std::string nm = "",
                       const std::string& logging_prefix = "")
    : BaseType(left->parameter_type() + right->parameter_type(),
               logging_prefix.empty() ? Common::to_camel_case(GetCombination<combination>::name() + "GridFunction")
                                      : logging_prefix,
               logging_prefix.empty())
    , left_(std::move(left))
    , right_(std::move(right))
    , name_(nm.empty()
                ? "(" + left_.access().name() + GetCombination<combination>::symbol() + right_.access().name() + ")"
                : nm)
  {
    LOG_(debug) << Common::to_camel_case(GetCombination<combination>::name() + "GridFunction")
                << "(left_raw_ptr=" << left << ", right_raw_ptr=" << right << ", nm=\"" << nm << "\")" << std::endl;
  }

  CombinedGridFunction(const ThisType& other) = default;

  CombinedGridFunction(ThisType&& source) = default;

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    LOG_(debug) << Common::to_camel_case(GetCombination<combination>::name() + "GridFunction") + "::local_function()"
                << std::endl;
    using LeftLF = typename LeftType::LocalFunctionType;
    using RightLF = typename RightType::LocalFunctionType;
    return std::make_unique<CombinedElementFunction<LeftLF, RightLF, combination>>(
        std::move(left_.access().local_function()), std::move(right_.access().local_function()));
  } // ... local_function(...)

  std::string name() const override final
  {
    return name_;
  }

private:
  Common::ConstStorageProvider<LeftType> left_;
  Common::ConstStorageProvider<RightType> right_;
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
template <class LeftSummandType, class RightSummandType>
class ProductGridFunction : public CombinedGridFunction<LeftSummandType, RightSummandType, CombinationType::product>
{
  using BaseType = CombinedGridFunction<LeftSummandType, RightSummandType, CombinationType::product>;

public:
  template <class... Args>
  explicit ProductGridFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
}; // class ProductGridFunction


template <class T1, class T2, class... Args>
std::shared_ptr<DifferenceGridFunction<T1, T2>> make_difference(const T1& left, const T2& right, Args&&... args)
{
  return std::make_shared<DifferenceGridFunction<T1, T2>>(left, right, std::forward<Args>(args)...);
}

template <class T1, class T2, class... Args>
std::shared_ptr<DifferenceGridFunction<T1, T2>>
make_difference(std::shared_ptr<T1> left, std::shared_ptr<T2> right, Args&&... args)
{
  return std::make_shared<DifferenceGridFunction<T1, T2>>(left, right, std::forward<Args>(args)...);
}

template <class T1, class T2, class... Args>
std::shared_ptr<DifferenceGridFunction<T1, T2>> make_difference(T1*&& left, T2*&& right, Args&&... args)
{
  return std::make_shared<DifferenceGridFunction<T1, T2>>(
      std::move(left), std::move(right), std::forward<Args>(args)...);
}


template <class T1, class T2, class... Args>
std::shared_ptr<SumGridFunction<T1, T2>> make_sum(const T1& left, const T2& right, Args&&... args)
{
  return std::make_shared<SumGridFunction<T1, T2>>(left, right, std::forward<Args>(args)...);
}

template <class T1, class T2, class... Args>
std::shared_ptr<SumGridFunction<T1, T2>> make_sum(std::shared_ptr<T1> left, std::shared_ptr<T2> right, Args&&... args)
{
  return std::make_shared<SumGridFunction<T1, T2>>(left, right, std::forward<Args>(args)...);
}

template <class T1, class T2, class... Args>
std::shared_ptr<SumGridFunction<T1, T2>> make_sum(T1*&& left, T2*&& right, Args&&... args)
{
  return std::make_shared<SumGridFunction<T1, T2>>(std::move(left), std::move(right), std::forward<Args>(args)...);
}


template <class T1, class T2, class... Args>
std::shared_ptr<FractionGridFunction<T1, T2>> make_fraction(const T1& left, const T2& right, Args&&... args)
{
  return std::make_shared<FractionGridFunction<T1, T2>>(left, right, std::forward<Args>(args)...);
}

template <class T1, class T2, class... Args>
std::shared_ptr<FractionGridFunction<T1, T2>>
make_fraction(std::shared_ptr<T1> left, std::shared_ptr<T2> right, Args&&... args)
{
  return std::make_shared<FractionGridFunction<T1, T2>>(left, right, std::forward<Args>(args)...);
}

template <class T1, class T2, class... Args>
std::shared_ptr<FractionGridFunction<T1, T2>> make_fraction(T1*&& left, T2*&& right, Args&&... args)
{
  return std::make_shared<FractionGridFunction<T1, T2>>(std::move(left), std::move(right), std::forward<Args>(args)...);
}


template <class T1, class T2, class... Args>
std::shared_ptr<ProductGridFunction<T1, T2>> make_product(const T1& left, const T2& right, Args&&... args)
{
  return std::make_shared<ProductGridFunction<T1, T2>>(left, right, std::forward<Args>(args)...);
}

template <class T1, class T2, class... Args>
std::shared_ptr<ProductGridFunction<T1, T2>>
make_product(std::shared_ptr<T1> left, std::shared_ptr<T2> right, Args&&... args)
{
  return std::make_shared<ProductGridFunction<T1, T2>>(left, right, std::forward<Args>(args)...);
}

template <class T1, class T2, class... Args>
std::shared_ptr<ProductGridFunction<T1, T2>> make_product(T1*&& left, T2*&& right, Args&&... args)
{
  return std::make_shared<ProductGridFunction<T1, T2>>(std::move(left), std::move(right), std::forward<Args>(args)...);
}


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_BASE_COMBINED_HH
