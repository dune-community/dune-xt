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

#ifndef DUNE_XT_FUNCTIONS_BASE_COMBINED_FUNCTIONS_HH
#define DUNE_XT_FUNCTIONS_BASE_COMBINED_FUNCTIONS_HH

#include <dune/xt/common/memory.hh>

#include <dune/xt/functions/interfaces/function.hh>
#include <dune/xt/functions/type_traits.hh>


namespace Dune {
namespace XT {
namespace Functions {
namespace internal {


/**
 * \brief Helper class defining types of combined functions, if available.
 *
 * \note Most likely you do not want to use this class directly, but Combined.
 *
 * \todo Update product handling as in CombinedElementFunctionHelper to allow for more combinations!
 */
template <class LeftType, class RightType, CombinationType comb>
class SelectCombined
{

public:
  using D = typename LeftType::DomainFieldType;
  static const size_t d = LeftType::domain_dim;
  using R = typename LeftType::RangeFieldType;

private:
  static_assert(std::is_same<typename RightType::DomainFieldType, D>::value, "Types do not match!");
  static_assert(RightType::domain_dim == d, "Dimensions do not match!");
  static_assert(std::is_same<typename RightType::RangeFieldType, R>::value, "Types do not match!");

  template <class L, class R>
  class Choose
  {
    //! last tpl arg cannot be dropped due to gcc bug https://gcc.gnu.org/bugzilla/show_bug.cgi?id=85282
    template <size_t rL, size_t rR, size_t rCL, size_t rcR, CombinationType cc, bool anything = true>
    class Dimension
    {
      static_assert(!anything, "No combination for these dimensions available!");
    };

    template <size_t r_in, size_t rC_in, bool anything>
    class Dimension<r_in, r_in, rC_in, rC_in, CombinationType::difference, anything>
    {
    public:
      static const size_t r = r_in;
      static const size_t rC = rC_in;
    };

    template <size_t r_in, size_t rC_in, bool anything>
    class Dimension<r_in, r_in, rC_in, rC_in, CombinationType::sum, anything>
    {
    public:
      static const size_t r = r_in;
      static const size_t rC = rC_in;
    };

    template <size_t r_in, size_t rC_in, bool anything>
    class Dimension<1, r_in, 1, rC_in, CombinationType::product, anything>
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

  using DomainType = typename FunctionInterface<d, r, rC, R>::DomainType;
  using RangeReturnType = typename RightType::RangeReturnType;
  using ScalarRangeReturnType = typename LeftType::RangeReturnType;
  using DerivativeRangeReturnType = typename FunctionInterface<d, r, rC, R>::DerivativeRangeReturnType;

private:
  template <CombinationType cc, bool anything = true>
  class Call
  {
    static_assert(!anything, "Nothing available for these combinations!");
  }; // class Call

  template <bool anything>
  class Call<CombinationType::difference, anything>
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

    static RangeReturnType evaluate(const LeftType& left_,
                                    const RightType& right_,
                                    const DomainType& point_in_global_coordinates,
                                    const Common::Parameter& param)
    {
      return left_.evaluate(point_in_global_coordinates, param) - right_.evaluate(point_in_global_coordinates, param);
    }

    static DerivativeRangeReturnType jacobian(const LeftType& left_,
                                              const RightType& right_,
                                              const DomainType& point_in_global_coordinates,
                                              const Common::Parameter& param)
    {
      return left_.jacobian(point_in_global_coordinates, param) - right_.jacobian(point_in_global_coordinates, param);
    } // ... jacobian(...)
  }; // class Call< ..., difference >

  template <bool anything>
  class Call<CombinationType::sum, anything>
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

    static RangeReturnType evaluate(const LeftType& left_,
                                    const RightType& right_,
                                    const DomainType& point_in_global_coordinates,
                                    const Common::Parameter& param)
    {
      return left_.evaluate(point_in_global_coordinates, param) + right_.evaluate(point_in_global_coordinates, param);
    } // ... evaluate(...)

    static DerivativeRangeReturnType jacobian(const LeftType& left_,
                                              const RightType& right_,
                                              const DomainType& point_in_global_coordinates,
                                              const Common::Parameter& param)
    {
      return left_.jacobian(point_in_global_coordinates, param) + right_.jacobian(point_in_global_coordinates, param);
    } // ... jacobian(...)
  }; // class Call< ..., sum >

  // left only scalar atm
  template <bool anything>
  class Call<CombinationType::product, anything>
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

    static RangeReturnType evaluate(const LeftType& left_,
                                    const RightType& right_,
                                    const DomainType& point_in_global_coordinates,
                                    const Common::Parameter& param)
    {
      ScalarRangeReturnType left_eval = left_.evaluate(point_in_global_coordinates, param);
      RangeReturnType right_eval = right_.evaluate(point_in_global_coordinates, param);
      if (left_eval.size() != 1)
        DUNE_THROW(NotImplemented, "Only available for scalar left type!");
      right_eval *= left_eval[0];
      return right_eval;
    } // ... evaluate(...)

    static DerivativeRangeReturnType jacobian(const LeftType& /*left_*/,
                                              const RightType& /*right_*/,
                                              const DomainType& /*point_in_global_coordinates*/,
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

  static RangeReturnType evaluate(const LeftType& left_,
                                  const RightType& right_,
                                  const DomainType& point_in_global_coordinates,
                                  const Common::Parameter& param)
  {
    return Call<comb>::evaluate(left_, right_, point_in_global_coordinates, param);
  }

  static DerivativeRangeReturnType jacobian(const LeftType& left_,
                                            const RightType& right_,
                                            const DomainType& point_in_global_coordinates,
                                            const Common::Parameter& param)
  {
    return Call<comb>::jacobian(left_, right_, point_in_global_coordinates, param);
  }
}; // class SelectCombined


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
template <class LeftType, class RightType, CombinationType comb>
class CombinedFunction
  : public FunctionInterface<LeftType::domain_dim,
                             SelectCombined<LeftType, RightType, comb>::r,
                             SelectCombined<LeftType, RightType, comb>::rC,
                             typename SelectCombined<LeftType, RightType, comb>::R>
{
  using BaseType = FunctionInterface<LeftType::domain_dim,
                                     SelectCombined<LeftType, RightType, comb>::r,
                                     SelectCombined<LeftType, RightType, comb>::rC,
                                     typename SelectCombined<LeftType, RightType, comb>::R>;

  using ThisType = CombinedFunction;

  using Select = SelectCombined<LeftType, RightType, comb>;

  using LeftStorageType = Common::ConstStorageProvider<LeftType>;
  using RightStorageType = Common::ConstStorageProvider<RightType>;

public:
  CombinedFunction(const LeftType& left, const RightType& right, const std::string nm = "")
    : left_(std::make_unique<LeftStorageType>(left))
    , right_(std::make_unique<RightStorageType>(right))
    , name_(get_name(left, right, nm))
  {}

  CombinedFunction(const std::shared_ptr<const LeftType> left,
                   const std::shared_ptr<const RightType> right,
                   const std::string nm = "")
    : left_(std::make_unique<LeftStorageType>(left))
    , right_(std::make_unique<RightStorageType>(right))
    , name_(get_name(*left, *right, nm))
  {}

  CombinedFunction(const LeftType& left, const std::shared_ptr<const RightType> right, const std::string nm = "")
    : left_(std::make_unique<LeftStorageType>(left))
    , right_(std::make_unique<RightStorageType>(right))
    , name_(get_name(left, *right, nm))
  {}

  CombinedFunction(const std::shared_ptr<const LeftType> left, const RightType& right, const std::string nm = "")
    : left_(std::make_unique<LeftStorageType>(left))
    , right_(std::make_unique<RightStorageType>(right))
    , name_(get_name(*left, right, nm))
  {}

  CombinedFunction(LeftType*&& left, RightType*&& right, const std::string nm = "")
    : left_(std::make_unique<LeftStorageType>(std::move(left)))
    , right_(std::make_unique<RightStorageType>(std::move(right)))
    , name_(nm.empty() ? SelectCombined<LeftType, RightType, comb>::type() + " of '" + left_->access().name()
                             + "' and '" + right_->access().name() + "'"
                       : nm)
  {}

  CombinedFunction(ThisType&& source) = default;

  CombinedFunction(const ThisType& other) = delete;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& other) = delete;

  std::string name() const override final
  {
    return name_;
  }

  using typename BaseType::DerivativeRangeReturnType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeReturnType;

  int order(const XT::Common::Parameter& param = {}) const override final
  {
    auto ret = Select::order(left_->access().order(param), right_->access().order(param));
    assert(ret < std::numeric_limits<int>::max());
    return static_cast<int>(ret);
  }

  RangeReturnType evaluate(const DomainType& point_in_global_coordinates,
                           const Common::Parameter& param = {}) const override final
  {
    return Select::evaluate(left_->access(), right_->access(), point_in_global_coordinates, param);
  }

  DerivativeRangeReturnType jacobian(const DomainType& point_in_global_coordinates,
                                     const Common::Parameter& param = {}) const override final
  {
    return Select::jacobian(left_->access(), right_->access(), point_in_global_coordinates, param);
  }

private:
  static std::string get_name(const LeftType& left, const RightType& right, const std::string& nm)
  {
    return nm.empty() ? SelectCombined<LeftType, RightType, comb>::type() + " of '" + left.name() + "' and '"
                            + right.name() + "'"
                      : nm;
  }

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
class DifferenceFunction : public internal::CombinedFunction<MinuendType, SubtrahendType, CombinationType::difference>
{
  using BaseType = internal::CombinedFunction<MinuendType, SubtrahendType, CombinationType::difference>;

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
class SumFunction : public internal::CombinedFunction<LeftSummandType, RightSummandType, CombinationType::sum>
{
  using BaseType = internal::CombinedFunction<LeftSummandType, RightSummandType, CombinationType::sum>;

public:
  template <class... Args>
  explicit SumFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
}; // class SumFunction


/**
 * \brief Function representing the product of two functions.
 *
 * \see internal::Combined
 */
template <class LeftSummandType, class RightSummandType>
class ProductFunction : public internal::CombinedFunction<LeftSummandType, RightSummandType, CombinationType::product>
{
  using BaseType = internal::CombinedFunction<LeftSummandType, RightSummandType, CombinationType::product>;

public:
  template <class... Args>
  explicit ProductFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
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

#endif // DUNE_XT_FUNCTIONS_BASE_COMBINED_HH
