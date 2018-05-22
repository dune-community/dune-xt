// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2015 - 2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_FUNCTIONS_DERIVED_HH
#define DUNE_XT_FUNCTIONS_DERIVED_HH
#if 0
#include <memory>
#include <type_traits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/memory.hh>

#include <dune/xt/functions/interfaces.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace internal {


enum class Derivative
{
  divergence
};


template <class FunctionType, Derivative derivative>
class SelectDerived
{
  static_assert(is_localizable_function<FunctionType>::value, "FunctionType has to be a LocalizableFunction!");

public:
  typedef typename FunctionType::EntityType E;
  typedef typename FunctionType::DomainFieldType D;
  static const size_t d = FunctionType::dimDomain;
  typedef typename FunctionType::RangeFieldType R;

private:
  template <class F>
  class Choose
  {
    template <size_t d, size_t r, size_t rC, Derivative der, bool anything = true>
    class Dimension
    {
      static_assert(!anything, "No derivative for these dimensions available!");
    };

    template <size_t d, bool anything>
    class Dimension<d, d, 1, Derivative::divergence, anything>
    {
    public:
      static const size_t r = 1;
      static const size_t rC = 1;
    };

  public:
    static const size_t r = Dimension<d, F::dimRange, F::dimRangeCols, derivative>::r;
    static const size_t rC = Dimension<d, F::dimRange, F::dimRangeCols, derivative>::rC;
  }; // class SelectDerived

public:
  static const size_t r = Choose<FunctionType>::r;
  static const size_t rC = Choose<FunctionType>::rC;

  typedef typename FunctionType::LocalfunctionType FunctionLocalfunctionType;
  typedef typename LocalfunctionInterface<E, D, d, R, r, rC>::DomainType DomainType;
  typedef typename LocalfunctionInterface<E, D, d, R, r, rC>::RangeType RangeType;
  typedef typename LocalfunctionInterface<E, D, d, R, r, rC>::JacobianRangeType JacobianRangeType;

private:
  template <Derivative der, bool anything = true>
  class Call
  {
    static_assert(!anything, "Nothing available for this derivative!");
  }; // class Call

  template <bool anything>
  class Call<Derivative::divergence, anything>
  {
  public:
    static std::string type()
    {
      return "divergence";
    }

    static size_t order(const size_t ord)
    {
      return boost::numeric_cast<size_t>(std::max(boost::numeric_cast<ssize_t>(ord) - 1, ssize_t(0)));
    }

    static void evaluate(const FunctionLocalfunctionType& func_local,
                         const DomainType& xx,
                         RangeType& ret,
                         const Common::Parameter& mu)
    {
      typename FunctionLocalfunctionType::JacobianRangeType tmp_jac(0.0);
      func_local.jacobian(xx, tmp_jac, mu);
      ret *= 0.0;
      for (size_t dd = 0; dd < d; ++dd)
        ret[0] += tmp_jac[dd][dd];
    } // ... evaluate(...)

    static void jacobian(const FunctionLocalfunctionType& /*func_local*/,
                         const DomainType& /*xx*/,
                         JacobianRangeType& /*ret*/,
                         const Common::Parameter& /*mu*/)
    {
      DUNE_THROW(NotImplemented, "for divergence!");
    }
  }; // class Call< ..., divergence >

public:
  static std::string type()
  {
    return Call<derivative>::type();
  }

  static size_t order(const size_t ord)
  {
    return Call<derivative>::order(ord);
  }

  static void evaluate(const FunctionLocalfunctionType& func_local,
                       const DomainType& xx,
                       RangeType& ret,
                       const Common::Parameter& mu)
  {
    Call<derivative>::evaluate(func_local, xx, ret, mu);
  }

  static void jacobian(const FunctionLocalfunctionType& func_local,
                       const DomainType& xx,
                       JacobianRangeType& ret,
                       const Common::Parameter& mu)
  {
    Call<derivative>::jacobian(func_local, xx, ret, mu);
  }
}; // class SelectDerived


template <class FunctionType, Derivative derivative>
class DerivedLocalFunction : public LocalfunctionInterface<typename SelectDerived<FunctionType, derivative>::E,
                                                           typename SelectDerived<FunctionType, derivative>::D,
                                                           SelectDerived<FunctionType, derivative>::d,
                                                           typename SelectDerived<FunctionType, derivative>::R,
                                                           SelectDerived<FunctionType, derivative>::r,
                                                           SelectDerived<FunctionType, derivative>::rC>
{
  typedef LocalfunctionInterface<typename SelectDerived<FunctionType, derivative>::E,
                                 typename SelectDerived<FunctionType, derivative>::D,
                                 SelectDerived<FunctionType, derivative>::d,
                                 typename SelectDerived<FunctionType, derivative>::R,
                                 SelectDerived<FunctionType, derivative>::r,
                                 SelectDerived<FunctionType, derivative>::rC>
      BaseType;

  typedef SelectDerived<FunctionType, derivative> Select;

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  DerivedLocalFunction(const FunctionType& func, const EntityType& ent)
    : BaseType(ent)
    , func_local_(func.local_function(this->entity()))
  {
  }

  virtual size_t order(const XT::Common::Parameter& mu = {}) const override final
  {
    return Select::order(func_local_->order(mu));
  }

  virtual void evaluate(const DomainType& xx, RangeType& ret, const Common::Parameter& mu = {}) const override final
  {
    Select::evaluate(*func_local_, xx, ret, mu);
  }

  virtual void
  jacobian(const DomainType& xx, JacobianRangeType& ret, const Common::Parameter& mu = {}) const override final
  {
    Select::jacobian(*func_local_, xx, ret, mu);
  }

private:
  const std::unique_ptr<const typename FunctionType::LocalfunctionType> func_local_;
}; // class DerivedLocalFunction


template <class FunctionType, Derivative derivative>
class Derived : public LocalizableFunctionInterface<typename SelectDerived<FunctionType, derivative>::E,
                                                    typename SelectDerived<FunctionType, derivative>::D,
                                                    SelectDerived<FunctionType, derivative>::d,
                                                    typename SelectDerived<FunctionType, derivative>::R,
                                                    SelectDerived<FunctionType, derivative>::r,
                                                    SelectDerived<FunctionType, derivative>::rC>
{
  typedef LocalizableFunctionInterface<typename SelectDerived<FunctionType, derivative>::E,
                                       typename SelectDerived<FunctionType, derivative>::D,
                                       SelectDerived<FunctionType, derivative>::d,
                                       typename SelectDerived<FunctionType, derivative>::R,
                                       SelectDerived<FunctionType, derivative>::r,
                                       SelectDerived<FunctionType, derivative>::rC>
      BaseType;
  typedef Common::ConstStorageProvider<FunctionType> FunctionStorageType;
  typedef Derived<FunctionType, derivative> ThisType;

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  Derived(const FunctionType& func, const std::string nm = "")
    : func_(Common::make_unique<FunctionStorageType>(func))
    , name_(nm.empty() ? SelectDerived<FunctionType, derivative>::type() + " of '" + func.name() + "'" : nm)
  {
  }

  Derived(const std::shared_ptr<const FunctionType> func, const std::string nm = "")
    : func_(Common::make_unique<FunctionStorageType>(func))
    , name_(nm.empty() ? SelectDerived<FunctionType, derivative>::type() + " of '" + func_->access().name() + "'" : nm)
  {
  }

  Derived(ThisType&& source) = default;
  Derived(const ThisType& other) = delete;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& other) = delete;

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override final
  {
    typedef DerivedLocalFunction<FunctionType, derivative> RealLocalFunctionType;
    assert(func_);
    return Common::make_unique<RealLocalFunctionType>(func_->access(), entity);
  } // ... local_function(...)

  virtual ThisType* copy() const
  {
    DUNE_THROW(NotImplemented, "Are you kidding me?");
  }

  virtual std::string type() const override final
  {
    return SelectDerived<FunctionType, derivative>::type() + " of '" + func_->access().type() + "'";
  }

  virtual std::string name() const override final
  {
    return name_;
  }

private:
  std::unique_ptr<const FunctionStorageType> func_;
  const std::string name_;
}; // class Derived


} // namespace internal


template <class FunctionType>
class DivergenceFunction : public internal::Derived<FunctionType, internal::Derivative::divergence>
{
  typedef internal::Derived<FunctionType, internal::Derivative::divergence> BaseType;

public:
  template <class... Args>
  DivergenceFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
}; // class DivergenceFunction


template <class T, class... Args>
std::shared_ptr<DivergenceFunction<T>> make_divergence(const T& func, Args&&... args)
{
  return std::make_shared<DivergenceFunction<T>>(func, std::forward<Args>(args)...);
}

template <class T, class... Args>
std::shared_ptr<DivergenceFunction<T>> make_divergence(std::shared_ptr<T> func, Args&&... args)
{
  return std::make_shared<DivergenceFunction<T>>(func, std::forward<Args>(args)...);
}


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif
#endif // DUNE_XT_FUNCTIONS_DERIVED_HH
