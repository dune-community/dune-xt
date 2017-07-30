// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014 - 2016)
//   Tobias Leibner  (2014, 2017)

#ifndef DUNE_XT_FUNCTIONS_COMBINED_HH
#define DUNE_XT_FUNCTIONS_COMBINED_HH

#include <memory>
#include <type_traits>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/memory.hh>

#include "interfaces.hh"

namespace Dune {
namespace XT {
namespace Functions {
namespace internal {


enum class Combination
{
  difference,
  sum,
  product
}; // enum class Combination


/**
 * \brief Helper class defining types of combined functions, if available.
 *
 * \note Most likely you do not want to use this class directly, but Combined.
 */
template <class LeftType, class RightType, Combination comb>
class SelectCombined
{
  static_assert(is_localizable_function<LeftType>::value, "LeftType has to be a LocalizableFunction!");
  static_assert(is_localizable_function<RightType>::value, "RightType has to be a LocalizableFunction!");

public:
  typedef typename LeftType::EntityType E;
  typedef typename LeftType::DomainFieldType D;
  static const size_t d = LeftType::dimDomain;
  typedef typename LeftType::RangeFieldType R;

private:
  static_assert(std::is_same<typename RightType::EntityType, E>::value, "Types do not match!");
  static_assert(std::is_same<typename RightType::DomainFieldType, D>::value, "Types do not match!");
  static_assert(RightType::dimDomain == d, "Dimensions do not match!");
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
    static const size_t r = Dimension<L::dimRange, R::dimRange, L::dimRangeCols, R::dimRangeCols, comb>::r;
    static const size_t rC = Dimension<L::dimRange, R::dimRange, L::dimRangeCols, R::dimRangeCols, comb>::rC;
  }; // class Choose

public:
  static const size_t r = Choose<LeftType, RightType>::r;
  static const size_t rC = Choose<LeftType, RightType>::rC;

  typedef typename LeftType::LocalfunctionType LeftLocalfunctionType;
  typedef typename RightType::LocalfunctionType RightLocalfunctionType;
  typedef typename LocalfunctionInterface<E, D, d, R, r, rC>::DomainType DomainType;
  typedef typename LocalfunctionInterface<E, D, d, R, r, rC>::RangeType RangeType;
  typedef typename LocalfunctionInterface<E, D, d, R, r, rC>::JacobianRangeType JacobianRangeType;

private:
  template <Combination cc, bool anything = true>
  class Call
  {
    static_assert(!anything, "Nothing available for these combinations!");
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

    static void evaluate(const LeftLocalfunctionType& left_local,
                         const RightLocalfunctionType& right_local,
                         const DomainType& xx,
                         RangeType& ret,
                         const Common::Parameter& mu,
                         RangeType& tmp_ret)
    {
      left_local.evaluate(xx, ret, mu);
      right_local.evaluate(xx, tmp_ret, mu);
      ret -= tmp_ret;
    } // ... evaluate(...)

    static void jacobian(const LeftLocalfunctionType& left_local,
                         const RightLocalfunctionType& right_local,
                         const DomainType& xx,
                         JacobianRangeType& ret,
                         const Common::Parameter& mu,
                         JacobianRangeType& tmp_ret)
    {
      left_local.jacobian(xx, ret, mu);
      right_local.jacobian(xx, tmp_ret, mu);
      ret -= tmp_ret;
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

    static void evaluate(const LeftLocalfunctionType& left_local,
                         const RightLocalfunctionType& right_local,
                         const DomainType& xx,
                         RangeType& ret,
                         const Common::Parameter& mu,
                         RangeType& tmp_ret)
    {
      left_local.evaluate(xx, ret, mu);
      right_local.evaluate(xx, tmp_ret, mu);
      ret += tmp_ret;
    } // ... evaluate(...)

    static void jacobian(const LeftLocalfunctionType& left_local,
                         const RightLocalfunctionType& right_local,
                         const DomainType& xx,
                         JacobianRangeType& ret,
                         const Common::Parameter& mu,
                         JacobianRangeType& tmp_ret)
    {
      left_local.jacobian(xx, ret, mu);
      right_local.jacobian(xx, tmp_ret, mu);
      ret += tmp_ret;
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

    static void evaluate(const LeftLocalfunctionType& left_local,
                         const RightLocalfunctionType& right_local,
                         const DomainType& xx,
                         RangeType& ret,
                         const Common::Parameter& mu,
                         RangeType& /*tmp_ret*/)
    {
      auto left_value = left_local.evaluate(xx, mu);
      right_local.evaluate(xx, ret, mu);
      ret *= left_value;
    } // ... evaluate(...)

    static void jacobian(const LeftLocalfunctionType& /*left_local*/,
                         const RightLocalfunctionType& /*right_local*/,
                         const DomainType& /*xx*/,
                         JacobianRangeType& /*ret*/,
                         const Common::Parameter& /*mu*/,
                         JacobianRangeType& /*tmp_ret*/)
    {
      DUNE_THROW(NotImplemented, "If you need this, implement it!");
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

  static void evaluate(const LeftLocalfunctionType& left_local,
                       const RightLocalfunctionType& right_local,
                       const DomainType& xx,
                       RangeType& ret,
                       const Common::Parameter& mu,
                       RangeType& tmp_ret)
  {
    Call<comb>::evaluate(left_local, right_local, xx, ret, mu, tmp_ret);
  }

  static void jacobian(const LeftLocalfunctionType& left_local,
                       const RightLocalfunctionType& right_local,
                       const DomainType& xx,
                       JacobianRangeType& ret,
                       const Common::Parameter& mu,
                       JacobianRangeType& tmp_ret)
  {
    Call<comb>::jacobian(left_local, right_local, xx, ret, mu, tmp_ret);
  }
}; // class SelectCombined

/**
 * \brief Generic combined local function.
 *
 * \note Most likely you do not want to use this class directly, but Combined.
 */
template <class LeftType, class RightType, Combination type>
class CombinedLocalFunction : public LocalfunctionInterface<typename SelectCombined<LeftType, RightType, type>::E,
                                                            typename SelectCombined<LeftType, RightType, type>::D,
                                                            SelectCombined<LeftType, RightType, type>::d,
                                                            typename SelectCombined<LeftType, RightType, type>::R,
                                                            SelectCombined<LeftType, RightType, type>::r,
                                                            SelectCombined<LeftType, RightType, type>::rC>
{
  typedef LocalfunctionInterface<typename SelectCombined<LeftType, RightType, type>::E,
                                 typename SelectCombined<LeftType, RightType, type>::D,
                                 SelectCombined<LeftType, RightType, type>::d,
                                 typename SelectCombined<LeftType, RightType, type>::R,
                                 SelectCombined<LeftType, RightType, type>::r,
                                 SelectCombined<LeftType, RightType, type>::rC>
      BaseType;

  typedef SelectCombined<LeftType, RightType, type> Select;

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  CombinedLocalFunction(const LeftType& left, const RightType& right, const EntityType& ent)
    : BaseType(ent)
    , left_local_(left.local_function(this->entity()))
    , right_local_(right.local_function(this->entity()))
    , tmp_range_(0.0)
    , tmp_jacobian_(0.0)
  {
  }

  virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    return Select::order(left_local_->order(), right_local_->order());
  }

  virtual void evaluate(const DomainType& xx, RangeType& ret, const Common::Parameter& mu = {}) const override final
  {
    Select::evaluate(*left_local_, *right_local_, xx, ret, mu, tmp_range_);
  }

  virtual void
  jacobian(const DomainType& xx, JacobianRangeType& ret, const Common::Parameter& mu = {}) const override final
  {
    Select::jacobian(*left_local_, *right_local_, xx, ret, mu, tmp_jacobian_);
  }

private:
  const std::unique_ptr<const typename LeftType::LocalfunctionType> left_local_;
  const std::unique_ptr<const typename RightType::LocalfunctionType> right_local_;
  mutable RangeType tmp_range_;
  mutable JacobianRangeType tmp_jacobian_;
}; // class CombinedLocalFunction

/**
 * \brief Generic combined function.
 *
 *        This class combines two given functions of type LeftType and RightType using the given combination
 *        Combination. This class (and any derived class, like Difference, Sum or Product) can be used in two was:
 *        - You can pass references of the left and right operand to this class. This is done for instance when calling
 *          operator+, operator- or operator* on any function deriving from LocalizableFunctionInterface:
\code
typedef Functions::Constant< ..., double, 2, double 1 > ConstantType;
ConstantType one(1);
ConstantType two(2);
// the following code
auto difference = one - two;
// is equivalent to
Difference< ConstantType, ConstantType > difference(one, two);
// and
internal::Combined< ConstantType, ConstantType, Combination::difference > difference(one, tow);
\endcode
 *          In this situation you are responsible to ensure that the arguments given are valid throughout the lifetime
 *          of this class. The following will lead to a segfault:
\code
typedef Functions::Constant< ..., double, 2, double 1 > ConstantType;

Difference< ConstantType, ConstantType > stupid_difference()
{
  ConstantType one(1);
  ConstantType two(2);
  return one - two;
}
\endcode
 *        - You can pass shared_ptr of the left and right operands to this class. In this case the following is valid:
\code
typedef Functions::Constant< ..., double, 2, double 1 > ConstantType;

Difference< ConstantType, ConstantType > stupid_difference()
{
  auto one = std::make_shared< ConstantType >(1);
  auto two = std::make_shared< ConstantType >(2);
  return Difference< ConstantType, ConstantType >(one, two)
}
\endcode
 *
 * \note  Most likely you do not want to use this class diretly, but one of Difference, Sum or Product.
 */
template <class LeftType, class RightType, Combination comb>
class Combined : public LocalizableFunctionInterface<typename SelectCombined<LeftType, RightType, comb>::E,
                                                     typename SelectCombined<LeftType, RightType, comb>::D,
                                                     SelectCombined<LeftType, RightType, comb>::d,
                                                     typename SelectCombined<LeftType, RightType, comb>::R,
                                                     SelectCombined<LeftType, RightType, comb>::r,
                                                     SelectCombined<LeftType, RightType, comb>::rC>
{
  typedef LocalizableFunctionInterface<typename SelectCombined<LeftType, RightType, comb>::E,
                                       typename SelectCombined<LeftType, RightType, comb>::D,
                                       SelectCombined<LeftType, RightType, comb>::d,
                                       typename SelectCombined<LeftType, RightType, comb>::R,
                                       SelectCombined<LeftType, RightType, comb>::r,
                                       SelectCombined<LeftType, RightType, comb>::rC>
      BaseType;
  typedef Common::ConstStorageProvider<LeftType> LeftStorageType;
  typedef Common::ConstStorageProvider<RightType> RightStorageType;
  typedef Combined<LeftType, RightType, comb> ThisType;

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  Combined(const LeftType& left, const RightType& right, const std::string nm = "")
    : left_(Common::make_unique<LeftStorageType>(left))
    , right_(Common::make_unique<RightStorageType>(right))
    , name_(nm.empty()
                ? SelectCombined<LeftType, RightType, comb>::type() + " of '" + left.name() + "' and '" + right.name()
                      + "'"
                : nm)
  {
  }

  Combined(const std::shared_ptr<const LeftType> left,
           const std::shared_ptr<const RightType> right,
           const std::string nm = "")
    : left_(Common::make_unique<LeftStorageType>(left))
    , right_(Common::make_unique<RightStorageType>(right))
    , name_(nm.empty()
                ? SelectCombined<LeftType, RightType, comb>::type() + " of '" + left_->access().name() + "' and '"
                      + right_->access().name()
                      + "'"
                : nm)
  {
  }

  Combined(ThisType&& source) = default;

  Combined(const ThisType& other) = delete;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& other) = delete;

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override final
  {
    typedef CombinedLocalFunction<LeftType, RightType, comb> RealLocalFunctionType;
    assert(left_);
    assert(right_);
    return Common::make_unique<RealLocalFunctionType>(left_->access(), right_->access(), entity);
  } // ... local_function(...)

  virtual ThisType* copy() const
  {
    DUNE_THROW(NotImplemented, "Are you kidding me?");
  }

  virtual std::string type() const override final
  {
    return SelectCombined<LeftType, RightType, comb>::type() + " of '" + left_->access().type() + "' and '"
           + right_->access().type() + "'";
  } // ... type(...)

  virtual std::string name() const override final
  {
    return name_;
  }

private:
  std::unique_ptr<const LeftStorageType> left_;
  std::unique_ptr<const RightStorageType> right_;
  const std::string name_;
}; // class Combined

} // namespace internal

/**
 * \brief Function representing the difference between two functions.
 *
 * \see internal::Combined
 */
template <class MinuendType, class SubtrahendType>
class DifferenceFunction : public internal::Combined<MinuendType, SubtrahendType, internal::Combination::difference>
{
  typedef internal::Combined<MinuendType, SubtrahendType, internal::Combination::difference> BaseType;

public:
  template <class... Args>
  DifferenceFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
}; // class DifferenceFunction

/**
 * \brief Function representing the sum of two functions.
 *
 * \see internal::Combined
 */
template <class LeftSummandType, class RightSummandType>
class SumFunction : public internal::Combined<LeftSummandType, RightSummandType, internal::Combination::sum>
{
  typedef internal::Combined<LeftSummandType, RightSummandType, internal::Combination::sum> BaseType;

public:
  template <class... Args>
  SumFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
}; // class SumFunction

/**
 * \brief Function representing the product of two functions.
 *
 * \see internal::Combined
 */
template <class LeftSummandType, class RightSummandType>
class ProductFunction : public internal::Combined<LeftSummandType, RightSummandType, internal::Combination::product>
{
  typedef internal::Combined<LeftSummandType, RightSummandType, internal::Combination::product> BaseType;

public:
  template <class... Args>
  ProductFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
}; // class ProductFunction

template <class T1, class T2, class... Args>
std::shared_ptr<DifferenceFunction<T1, T2>> make_difference(const T1& left, const T2& right, Args&&... args)
{
  return std::make_shared<DifferenceFunction<T1, T2>>(left, right, std::forward<Args>(args)...);
}

template <class T1, class T2, class... Args>
std::shared_ptr<DifferenceFunction<T1, T2>>
make_difference(std::shared_ptr<T1> left, std::shared_ptr<T2> right, Args&&... args)
{
  return std::make_shared<DifferenceFunction<T1, T2>>(left, right, std::forward<Args>(args)...);
}

template <class T1, class T2, class... Args>
std::shared_ptr<SumFunction<T1, T2>> make_sum(const T1& left, const T2& right, Args&&... args)
{
  return std::make_shared<SumFunction<T1, T2>>(left, right, std::forward<Args>(args)...);
}

template <class T1, class T2, class... Args>
std::shared_ptr<SumFunction<T1, T2>> make_sum(std::shared_ptr<T1> left, std::shared_ptr<T2> right, Args&&... args)
{
  return std::make_shared<SumFunction<T1, T2>>(left, right, std::forward<Args>(args)...);
}

template <class T1, class T2, class... Args>
std::shared_ptr<ProductFunction<T1, T2>> make_product(const T1& left, const T2& right, Args&&... args)
{
  return std::make_shared<ProductFunction<T1, T2>>(left, right, std::forward<Args>(args)...);
}

template <class T1, class T2, class... Args>
std::shared_ptr<ProductFunction<T1, T2>>
make_product(std::shared_ptr<T1> left, std::shared_ptr<T2> right, Args&&... args)
{
  return std::make_shared<ProductFunction<T1, T2>>(left, right, std::forward<Args>(args)...);
}


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_COMBINED_HH
