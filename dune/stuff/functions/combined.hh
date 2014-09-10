// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTIONS_COMBINED_HH
#define DUNE_STUFF_FUNCTIONS_COMBINED_HH

#include <type_traits>

#include <dune/common/typetraits.hh>

#include <dune/stuff/common/memory.hh>

#include "interfaces.hh"

namespace Dune {
namespace Stuff {
namespace Functions {
namespace internal {


enum class Combination
{
  difference,
  sum,
  product
}; // enum class Combination


template <class LeftType, class RightType, Combination comb>
class SelectCombined
{
  static_assert(std::is_base_of<Tags::LocalizableFunction, LeftType>::value,
                "LeftType has to be a LocalizableFunction!");
  static_assert(std::is_base_of<Tags::LocalizableFunction, RightType>::value,
                "RightType has to be a LocalizableFunction!");

public:
  typedef typename LeftType::EntityType E;
  typedef typename LeftType::DomainFieldType D;
  static const unsigned int d = LeftType::dimDomain;
  typedef typename LeftType::RangeFieldType R;

private:
  static_assert(std::is_same<typename RightType::EntityType, E>::value, "Types do not match!");
  static_assert(std::is_same<typename RightType::DomainFieldType, D>::value, "Types do not match!");
  static_assert(RightType::dimDomain == d, "Dimensions do not match!");
  static_assert(std::is_same<typename RightType::RangeFieldType, R>::value, "Types do not match!");

  template <class L, class R>
  class Choose
  {
    template <int rL, int rR, int rCL, int rcR, Combination cc, bool anything = true>
    class Dimension
    {
      static_assert(!anything, "No combination for these dimensions available!");
    };

    template <int r_in, int rC_in, bool anything>
    class Dimension<r_in, r_in, rC_in, rC_in, Combination::difference, anything>
    {
    public:
      static const unsigned int r  = r_in;
      static const unsigned int rC = rC_in;
    };

    template <int r_in, int rC_in, bool anything>
    class Dimension<r_in, r_in, rC_in, rC_in, Combination::sum, anything>
    {
    public:
      static const unsigned int r  = r_in;
      static const unsigned int rC = rC_in;
    };

    template <int r_in, int rC_in, bool anything>
    class Dimension<1, r_in, 1, rC_in, Combination::product, anything>
    {
    public:
      static const unsigned int r  = r_in;
      static const unsigned int rC = rC_in;
    };

  public:
    static const unsigned int r  = Dimension<L::dimRange, R::dimRange, L::dimRangeCols, R::dimRangeCols, comb>::r;
    static const unsigned int rC = Dimension<L::dimRange, R::dimRange, L::dimRangeCols, R::dimRangeCols, comb>::rC;
  }; // class Choose

public:
  static const unsigned int r  = Choose<LeftType, RightType>::r;
  static const unsigned int rC = Choose<LeftType, RightType>::r;

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

    static void evaluate(const LeftLocalfunctionType& left_local, const RightLocalfunctionType& right_local,
                         const DomainType& xx, RangeType& ret, RangeType& tmp_ret)
    {
      left_local.evaluate(xx, ret);
      right_local.evaluate(xx, tmp_ret);
      ret -= tmp_ret;
    } // ... evaluate(...)

    static void jacobian(const LeftLocalfunctionType& left_local, const RightLocalfunctionType& right_local,
                         const DomainType& xx, JacobianRangeType& ret, JacobianRangeType& tmp_ret)
    {
      left_local.jacobian(xx, ret);
      right_local.jacobian(xx, tmp_ret);
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

    static void evaluate(const LeftLocalfunctionType& left_local, const RightLocalfunctionType& right_local,
                         const DomainType& xx, RangeType& ret, RangeType& tmp_ret)
    {
      left_local.evaluate(xx, ret);
      right_local.evaluate(xx, tmp_ret);
      ret += tmp_ret;
    } // ... evaluate(...)

    static void jacobian(const LeftLocalfunctionType& left_local, const RightLocalfunctionType& right_local,
                         const DomainType& xx, JacobianRangeType& ret, JacobianRangeType& tmp_ret)
    {
      left_local.jacobian(xx, ret);
      right_local.jacobian(xx, tmp_ret);
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

    static void evaluate(const LeftLocalfunctionType& left_local, const RightLocalfunctionType& right_local,
                         const DomainType& xx, RangeType& ret, RangeType& /*tmp_ret*/)
    {
      auto left_value = left_local.evaluate(xx);
      right_local.evaluate(xx, ret);
      ret *= left_value;
    } // ... evaluate(...)

    static void jacobian(const LeftLocalfunctionType& /*left_local*/, const RightLocalfunctionType& /*right_local*/,
                         const DomainType& /*xx*/, JacobianRangeType& /*ret*/, JacobianRangeType& /*tmp_ret*/)
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

  static void evaluate(const LeftLocalfunctionType& left_local, const RightLocalfunctionType& right_local,
                       const DomainType& xx, RangeType& ret, RangeType& tmp_ret)
  {
    Call<comb>::evaluate(left_local, right_local, xx, ret, tmp_ret);
  }

  static void jacobian(const LeftLocalfunctionType& left_local, const RightLocalfunctionType& right_local,
                       const DomainType& xx, JacobianRangeType& ret, JacobianRangeType& tmp_ret)
  {
    Call<comb>::jacobian(left_local, right_local, xx, ret, tmp_ret);
  }
}; // class SelectCombined


template <class LeftType, class RightType, Combination type>
class CombinedLocalFunction
    : public LocalfunctionInterface<
          typename SelectCombined<LeftType, RightType, type>::E, typename SelectCombined<LeftType, RightType, type>::D,
          SelectCombined<LeftType, RightType, type>::d, typename SelectCombined<LeftType, RightType, type>::R,
          SelectCombined<LeftType, RightType, type>::r, SelectCombined<LeftType, RightType, type>::rC>
{
  typedef LocalfunctionInterface<
      typename SelectCombined<LeftType, RightType, type>::E, typename SelectCombined<LeftType, RightType, type>::D,
      SelectCombined<LeftType, RightType, type>::d, typename SelectCombined<LeftType, RightType, type>::R,
      SelectCombined<LeftType, RightType, type>::r, SelectCombined<LeftType, RightType, type>::rC> BaseType;

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

  virtual size_t order() const DS_OVERRIDE DS_FINAL
  {
    return Select::order(left_local_->order(), right_local_->order());
  }

  virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_OVERRIDE DS_FINAL
  {
    Select::evaluate(*left_local_, *right_local_, xx, ret, tmp_range_);
  }

  virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_OVERRIDE DS_FINAL
  {
    Select::jacobian(*left_local_, *right_local_, xx, ret, tmp_jacobian_);
  }

private:
  const std::unique_ptr<const typename LeftType::LocalfunctionType> left_local_;
  const std::unique_ptr<const typename RightType::LocalfunctionType> right_local_;
  mutable RangeType tmp_range_;
  mutable JacobianRangeType tmp_jacobian_;
}; // class CombinedLocalFunction


template <class LeftType, class RightType, Combination comb>
class Combined
    : public LocalizableFunctionInterface<
          typename SelectCombined<LeftType, RightType, comb>::E, typename SelectCombined<LeftType, RightType, comb>::D,
          SelectCombined<LeftType, RightType, comb>::d, typename SelectCombined<LeftType, RightType, comb>::R,
          SelectCombined<LeftType, RightType, comb>::r, SelectCombined<LeftType, RightType, comb>::rC>
{
  typedef LocalizableFunctionInterface<
      typename SelectCombined<LeftType, RightType, comb>::E, typename SelectCombined<LeftType, RightType, comb>::D,
      SelectCombined<LeftType, RightType, comb>::d, typename SelectCombined<LeftType, RightType, comb>::R,
      SelectCombined<LeftType, RightType, comb>::r, SelectCombined<LeftType, RightType, comb>::rC> BaseType;
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

  Combined(const std::shared_ptr<const LeftType> left, const std::shared_ptr<const RightType> right,
           const std::string nm = "")
    : left_(Common::make_unique<LeftStorageType>(left))
    , right_(Common::make_unique<RightStorageType>(right))
    , name_(nm.empty()
                ? SelectCombined<LeftType, RightType, comb>::type() + " of '" + left_->storage_access().name()
                      + "' and '" + right_->storage_access().name() + "'"
                : nm)
  {
  }

  Combined(ThisType&& source)
    : left_(std::move(source.left_))
    , right_(std::move(source.right_))
    , name_(source.name_)
  {
  }

  Combined(const ThisType& other) = delete;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& other) = delete;

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const /*DS_OVERRIDE DS_FINAL*/
  {
    typedef CombinedLocalFunction<LeftType, RightType, comb> RealLocalFunctionType;
    return std::unique_ptr<RealLocalFunctionType>(
        new RealLocalFunctionType(left_->storage_access(), right_->storage_access(), entity));
  } // ... local_function(...)

  virtual ThisType* copy() const
  {
    DUNE_THROW(NotImplemented, "Are you kidding me?");
  }

  virtual std::string type() const DS_OVERRIDE DS_FINAL
  {
    return SelectCombined<LeftType, RightType, comb>::type() + " of '" + left_->storage_access().type() + "' and '"
           + right_->storage_access().type() + "'";
  } // ... type(...)

  virtual std::string name() const DS_OVERRIDE DS_FINAL
  {
    return name_;
  }

private:
  std::unique_ptr<const LeftStorageType> left_;
  std::unique_ptr<const RightStorageType> right_;
  const std::string name_;
}; // class Combined


} // namespace internal


template <class MinuendType, class SubtrahendType>
class Difference : public internal::Combined<MinuendType, SubtrahendType, internal::Combination::difference>
{
  typedef internal::Combined<MinuendType, SubtrahendType, internal::Combination::difference> BaseType;

public:
  template <class... Args>
  Difference(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
}; // class Difference


template <class LeftSummandType, class RightSummandType>
class Sum : public internal::Combined<LeftSummandType, RightSummandType, internal::Combination::sum>
{
  typedef internal::Combined<LeftSummandType, RightSummandType, internal::Combination::sum> BaseType;

public:
  template <class... Args>
  Sum(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
}; // class Sum


template <class LeftSummandType, class RightSummandType>
class Product : public internal::Combined<LeftSummandType, RightSummandType, internal::Combination::product>
{
  typedef internal::Combined<LeftSummandType, RightSummandType, internal::Combination::product> BaseType;

public:
  template <class... Args>
  Product(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
}; // class Product


} // namespace Functions
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTIONS_COMBINED_HH
