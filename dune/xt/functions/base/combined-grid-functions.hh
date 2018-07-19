// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights
// reserved.
// License: Dual licensed as BSD 2-Clause License
// (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014 - 2018)
//   Tobias Leibner  (2014, 2017)

#ifndef DUNE_XT_FUNCTIONS_BASE_COMBINED_GRID_FUNCTIONS_HH
#define DUNE_XT_FUNCTIONS_BASE_COMBINED_GRID_FUNCTIONS_HH

#include <dune/xt/functions/base/combined-functions.hh>
#include <dune/xt/functions/interfaces/function.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace internal {


/**
 * \brief Helper class defining types of combined functions, if available.
 *
 * \note Most likely you do not want to use this class directly, but Combined.
 */
template <class LeftType, class RightType, Combination comb>
class SelectCombinedGridFunction
{
  static_assert(is_grid_function<LeftType>::value, "");
  static_assert(is_grid_function<RightType>::value, "");

public:
  using E = typename LeftType::ElementType;
  using D = typename LeftType::DomainFieldType;
  static const size_t d = LeftType::domain_dim;
  using R = typename LeftType::RangeFieldType;

private:
  static_assert(std::is_same<typename RightType::ElementType, E>::value, "Types do not match!");
  static_assert(std::is_same<typename RightType::DomainFieldType, D>::value, "Types do not match!");
  static_assert(RightType::domain_dim == d, "Dimensions do not match!");
  static_assert(std::is_same<typename RightType::RangeFieldType, R>::value, "Types do not match!");

  template <class L, class R>
  class Choose
  {
    template <size_t rL, size_t rR, size_t rCL, size_t rcR, Combination cc, bool anything = true>
    class Dimension
    {
      static_assert(!anything, "No combination for these dimensions available!");
    };

    template <size_t r_in, size_t rC_in, bool anything>
    class Dimension<r_in, r_in, rC_in, rC_in, Combination::difference, anything>
    {
    public:
      static const size_t r = r_in;
      static const size_t rC = rC_in;
    };

    template <size_t r_in, size_t rC_in, bool anything>
    class Dimension<r_in, r_in, rC_in, rC_in, Combination::sum, anything>
    {
    public:
      static const size_t r = r_in;
      static const size_t rC = rC_in;
    };

    template <size_t r_in, size_t rC_in, bool anything>
    class Dimension<1, r_in, 1, rC_in, Combination::product, anything>
    {
    public:
      static const size_t r = r_in;
      static const size_t rC = rC_in;
    };

  public:
    static const size_t r = Dimension<L::range_dim, R::range_dim, L::range_dim_cols, R::range_dim_cols, comb>::r;
    static const size_t rC = Dimension<L::range_dim, R::range_dim, L::range_dim_cols, R::range_dim_cols, comb>::rC;
  }; // class Choose

public:
  static const size_t r = Choose<LeftType, RightType>::r;
  static const size_t rC = Choose<LeftType, RightType>::rC;

  using LeftLocalFunctionType = typename LeftType::LocalFunctionType;
  using RightLocalFunctionType = typename RightType::LocalFunctionType;
  using DomainType = typename ElementFunctionInterface<E, r, rC, R>::DomainType;
  using RangeType = typename RightType::LocalFunctionType::RangeType;
  using ScalarRangeType = typename LeftType::LocalFunctionType::RangeType;
  using DerivativeRangeType = typename ElementFunctionInterface<E, r, rC, R>::DerivativeRangeType;
  using DerivativeRangeReturnType = typename ElementFunctionInterface<E, r, rC, R>::DerivativeRangeReturnType;

private:
  template <Combination cc, bool anything = true>
  class Call
  {
    static_assert(!anything, "Nothing available for these Combination!");
  }; // class Call

  template <bool anything>
  class Call<Combination::difference, anything>
  {
  public:
    static std::string type()
    {
      return "difference";
    }

    static size_t order(const size_t left_order, const size_t right_order)
    {
      return std::max(left_order, right_order);
    }

    static RangeType evaluate(const LeftLocalFunctionType& left_local,
                              const RightLocalFunctionType& right_local,
                              const DomainType& point_in_reference_element,
                              const Common::Parameter& param)
    {
      return left_local.evaluate(point_in_reference_element, param)
             - right_local.evaluate(point_in_reference_element, param);
    }

    static DerivativeRangeReturnType jacobian(const LeftLocalFunctionType& left_local,
                                              const RightLocalFunctionType& right_local,
                                              const DomainType& point_in_reference_element,
                                              const Common::Parameter& param)
    {
      return left_local.jacobian(point_in_reference_element, param)
             - right_local.jacobian(point_in_reference_element, param);
    } // ... jacobian(...)
  }; // class Call< ..., difference >

  template <bool anything>
  class Call<Combination::sum, anything>
  {
  public:
    static std::string type()
    {
      return "sum";
    }

    static size_t order(const size_t left_order, const size_t right_order)
    {
      return std::max(left_order, right_order);
    }

    static RangeType evaluate(const LeftLocalFunctionType& left_local,
                              const RightLocalFunctionType& right_local,
                              const DomainType& point_in_reference_element,
                              const Common::Parameter& param)
    {
      return left_local.evaluate(point_in_reference_element, param)
             + right_local.evaluate(point_in_reference_element, param);
    } // ... evaluate(...)

    static DerivativeRangeReturnType jacobian(const LeftLocalFunctionType& left_local,
                                              const RightLocalFunctionType& right_local,
                                              const DomainType& point_in_reference_element,
                                              const Common::Parameter& param)
    {
      return left_local.jacobian(point_in_reference_element, param)
             + right_local.jacobian(point_in_reference_element, param);
    } // ... jacobian(...)
  }; // class Call< ..., sum >

  // left only scalar atm
  template <bool anything>
  class Call<Combination::product, anything>
  {
  public:
    static std::string type()
    {
      return "product";
    }

    static size_t order(const size_t left_order, const size_t right_order)
    {
      return left_order + right_order;
    }

    static RangeType evaluate(const LeftLocalFunctionType& left_local,
                              const RightLocalFunctionType& right_local,
                              const DomainType& point_in_reference_element,
                              const Common::Parameter& param)
    {
      ScalarRangeType left_eval = left_local.evaluate(point_in_reference_element, param);
      RangeType right_eval = right_local.evaluate(point_in_reference_element, param);
      if (left_eval.size() != 1)
        DUNE_THROW(NotImplemented, "Only available for scalar left type!");
      right_eval *= left_eval[0];
      return right_eval;
    } // ... evaluate(...)

    static DerivativeRangeReturnType jacobian(const LeftLocalFunctionType& /*left_local*/,
                                              const RightLocalFunctionType& /*right_local*/,
                                              const DomainType& /*point_in_reference_element*/,
                                              const Common::Parameter& /*param*/)
    {
      DUNE_THROW(NotImplemented, "If you need this, implement it!");
      return DerivativeRangeReturnType();
    }
  }; // class Call< ..., product >

public:
  static std::string type()
  {
    return Call<comb>::type();
  }

  static size_t order(const size_t left_order, const size_t right_order)
  {
    return Call<comb>::order(left_order, right_order);
  }

  static RangeType evaluate(const LeftLocalFunctionType& left_local,
                            const RightLocalFunctionType& right_local,
                            const DomainType& point_in_reference_element,
                            const Common::Parameter& param)
  {
    return Call<comb>::evaluate(left_local, right_local, point_in_reference_element, param);
  }

  static DerivativeRangeReturnType jacobian(const LeftLocalFunctionType& left_local,
                                            const RightLocalFunctionType& right_local,
                                            const DomainType& point_in_reference_element,
                                            const Common::Parameter& param)
  {
    return Call<comb>::jacobian(left_local, right_local, point_in_reference_element, param);
  }
}; // class SelectCombinedGridFunction

/**
 * \brief Generic combined local function.
 *
 * \note Most likely you do not want to use this class directly, but Combined.
 */
template <class LeftType, class RightType, Combination type>
class CombinedLocalFunction
    : public ElementFunctionInterface<typename SelectCombinedGridFunction<LeftType, RightType, type>::E,
                                      SelectCombinedGridFunction<LeftType, RightType, type>::r,
                                      SelectCombinedGridFunction<LeftType, RightType, type>::rC,
                                      typename SelectCombinedGridFunction<LeftType, RightType, type>::R>
{
  using BaseType = ElementFunctionInterface<typename SelectCombinedGridFunction<LeftType, RightType, type>::E,
                                            SelectCombinedGridFunction<LeftType, RightType, type>::r,
                                            SelectCombinedGridFunction<LeftType, RightType, type>::rC,
                                            typename SelectCombinedGridFunction<LeftType, RightType, type>::R>;

  using Select = SelectCombinedGridFunction<LeftType, RightType, type>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::DerivativeRangeType;
  using typename BaseType::RangeReturnType;
  using typename BaseType::DerivativeRangeReturnType;

  CombinedLocalFunction(const LeftType& left, const RightType& right)
    : BaseType()
    , left_local_(left.local_function())
    , right_local_(right.local_function())
  {
  }

protected:
  void post_bind(const ElementType& element) override final
  {
    left_local_->bind(element);
    right_local_->bind(element);
  }

public:
  int order(const XT::Common::Parameter& param = {}) const override final
  {
    return Select::order(left_local_->order(param), right_local_->order(param));
  }

  RangeReturnType evaluate(const DomainType& point_in_reference_element,
                           const Common::Parameter& param = {}) const override final
  {
    return Select::evaluate(*left_local_, *right_local_, point_in_reference_element, param);
  }

  DerivativeRangeReturnType jacobian(const DomainType& point_in_reference_element,
                                     const Common::Parameter& param = {}) const override final
  {
    return Select::jacobian(*left_local_, *right_local_, point_in_reference_element, param);
  }

private:
  std::unique_ptr<typename LeftType::LocalFunctionType> left_local_;
  std::unique_ptr<typename RightType::LocalFunctionType> right_local_;
}; // class CombinedLocalFunction

/**
 * \brief Generic combined function.
 *
 *        This class combines two given functions of type LeftType and RightType
using the given combination
 *        Combination. This class (and any derived class, like Difference, Sum
or Product) can be used in two ways:
 *        - You can pass references of the left and right operand to this class.
This is done for instance when calling
 *          operator+, operator- or operator* on any function deriving from
GridFunctionInterface:
\code
using IndicatorType = Functions::IndicatorFunction< ..., double>;
IndicatorType one( ... );
IndicatorType two( ... );
// the following code
auto difference = one - two;
// is equivalent to
Difference< IndicatorType, IndicatorType > difference(one, two);
// and
internal::Combined< IndicatorType, IndicatorType, Combination::difference >
difference(one, tow);
\endcode
 *          In this situation you are responsible to ensure that the arguments
given are valid throughout the lifetime
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
 *        - You can pass shared_ptr of the left and right operands to this
class. In this case the following is valid:
\code
using IndicatorType = Functions::IndicatorFunction< ..., double >;

Difference< IndicatorType, IndicatorType > stupid_difference()
{
  auto one = std::make_shared< IndicatorType >(1);
  auto two = std::make_shared< IndicatorType >(2);
  return Difference< IndicatorType, IndicatorType >(one, two)
}
\endcode
 *
 * \note  Most likely you do not want to use this class diretly, but one of
Difference, Sum or Product.
 */
template <class LeftType, class RightType, Combination comb>
class CombinedGridFunction
    : public GridFunctionInterface<typename SelectCombinedGridFunction<LeftType, RightType, comb>::E,
                                   SelectCombinedGridFunction<LeftType, RightType, comb>::r,
                                   SelectCombinedGridFunction<LeftType, RightType, comb>::rC,
                                   typename SelectCombinedGridFunction<LeftType, RightType, comb>::R>
{
  using BaseType = GridFunctionInterface<typename SelectCombinedGridFunction<LeftType, RightType, comb>::E,
                                         SelectCombinedGridFunction<LeftType, RightType, comb>::r,
                                         SelectCombinedGridFunction<LeftType, RightType, comb>::rC,
                                         typename SelectCombinedGridFunction<LeftType, RightType, comb>::R>;

  using LeftStorageType = Common::ConstStorageProvider<LeftType>;
  using RightStorageType = Common::ConstStorageProvider<RightType>;
  using ThisType = CombinedGridFunction<LeftType, RightType, comb>;

public:
  using ElementType = typename BaseType::ElementType;
  using LocalFunctionType = typename BaseType::LocalFunctionType;

  CombinedGridFunction(const LeftType& left, const RightType& right, const std::string nm = "")
    : left_(Common::make_unique<LeftStorageType>(left))
    , right_(Common::make_unique<RightStorageType>(right))
    , name_(nm.empty()
                ? SelectCombinedGridFunction<LeftType, RightType, comb>::type() + " of '" + left.name() + "' and '"
                      + right.name()
                      + "'"
                : nm)
  {
  }

  CombinedGridFunction(const std::shared_ptr<const LeftType> left,
                       const std::shared_ptr<const RightType> right,
                       const std::string nm = "")
    : left_(Common::make_unique<LeftStorageType>(left))
    , right_(Common::make_unique<RightStorageType>(right))
    , name_(nm.empty()
                ? SelectCombinedGridFunction<LeftType, RightType, comb>::type() + " of '" + left_->access().name()
                      + "' and '"
                      + right_->access().name()
                      + "'"
                : nm)
  {
  }

  CombinedGridFunction(ThisType&& source) = default;

  CombinedGridFunction(const ThisType& other) = delete;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& other) = delete;

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    using RealLocalFunctionType = CombinedLocalFunction<LeftType, RightType, comb>;
    assert(left_);
    assert(right_);
    return Common::make_unique<RealLocalFunctionType>(left_->access(), right_->access());
  } // ... local_function(...)

  std::string type() const override final
  {
    return SelectCombinedGridFunction<LeftType, RightType, comb>::type() + " of '" + left_->access().type() + "' and '"
           + right_->access().type() + "'";
  } // ... type(...)

  std::string name() const override final
  {
    return name_;
  }

private:
  std::unique_ptr<const LeftStorageType> left_;
  std::unique_ptr<const RightStorageType> right_;
  const std::string name_;
}; // class CombinedGridFunction

} // namespace internal

/**
 * \brief Function representing the difference between two functions.
 *
 * \see internal::CombinedGridFunction
 */
template <class MinuendType, class SubtrahendType>
class DifferenceGridFunction
    : public internal::CombinedGridFunction<MinuendType, SubtrahendType, internal::Combination::difference>
{
  using BaseType = internal::CombinedGridFunction<MinuendType, SubtrahendType, internal::Combination::difference>;

public:
  template <class... Args>
  explicit DifferenceGridFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
}; // class DifferenceGridFunction

/**
 * \brief Function representing the sum of two functions.
 *
 * \see internal::CombinedGridFunction
 */
template <class LeftSummandType, class RightSummandType>
class SumGridFunction
    : public internal::CombinedGridFunction<LeftSummandType, RightSummandType, internal::Combination::sum>
{
  using BaseType = internal::CombinedGridFunction<LeftSummandType, RightSummandType, internal::Combination::sum>;

public:
  template <class... Args>
  explicit SumGridFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
}; // class SumGridFunction

/**
 * \brief Function representing the product of two functions.
 *
 * \see internal::CombinedGridFunction
 */
template <class LeftSummandType, class RightSummandType>
class ProductGridFunction
    : public internal::CombinedGridFunction<LeftSummandType, RightSummandType, internal::Combination::product>
{
  using BaseType = internal::CombinedGridFunction<LeftSummandType, RightSummandType, internal::Combination::product>;

public:
  template <class... Args>
  explicit ProductGridFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
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

} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_BASE_COMBINED_HH
