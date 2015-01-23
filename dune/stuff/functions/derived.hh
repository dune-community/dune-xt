// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTIONS_DERIVED_HH
#define DUNE_STUFF_FUNCTIONS_DERIVED_HH

#include <type_traits>
#include <memory>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/typetraits.hh>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/memory.hh>

#include "interfaces.hh"

namespace Dune {
namespace Stuff {
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
  static const unsigned int d = FunctionType::dimDomain;
  typedef typename FunctionType::RangeFieldType R;

private:
  template <class F>
  class Choose
  {
    template <int d, int r, int rC, Derivative der, bool anything = true>
    class Dimension
    {
      static_assert(!anything, "No derivative for these dimensions available!");
    };

    template <int d, bool anything>
    class Dimension<d, d, 1, Derivative::divergence, anything>
    {
    public:
      static const unsigned int r  = 1;
      static const unsigned int rC = 1;
    };

  public:
    static const unsigned int r  = Dimension<d, F::dimRange, F::dimRangeCols, derivative>::r;
    static const unsigned int rC = Dimension<d, F::dimRange, F::dimRangeCols, derivative>::rC;
  }; // class SelectDerived

public:
  static const unsigned int r  = Choose<FunctionType>::r;
  static const unsigned int rC = Choose<FunctionType>::rC;

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

    static void evaluate(const FunctionLocalfunctionType& func_local, const DomainType& xx, RangeType& ret)
    {
      typename FunctionLocalfunctionType::JacobianRangeType tmp_jac(0.0);
      func_local.jacobian(xx, tmp_jac);
      ret *= 0.0;
      for (size_t dd = 0; dd < d; ++dd)
        ret[0] += tmp_jac[dd][dd];
    } // ... evaluate(...)

    static void jacobian(const FunctionLocalfunctionType& /*func_local*/, const DomainType& /*xx*/,
                         JacobianRangeType& /*ret*/)
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

  static void evaluate(const FunctionLocalfunctionType& func_local, const DomainType& xx, RangeType& ret)
  {
    Call<derivative>::evaluate(func_local, xx, ret);
  }

  static void jacobian(const FunctionLocalfunctionType& func_local, const DomainType& xx, JacobianRangeType& ret)
  {
    Call<derivative>::jacobian(func_local, xx, ret);
  }
}; // class SelectDerived


template <class FunctionType, Derivative derivative>
class DerivedLocalFunction
    : public LocalfunctionInterface<
          typename SelectDerived<FunctionType, derivative>::E, typename SelectDerived<FunctionType, derivative>::D,
          SelectDerived<FunctionType, derivative>::d, typename SelectDerived<FunctionType, derivative>::R,
          SelectDerived<FunctionType, derivative>::r, SelectDerived<FunctionType, derivative>::rC>
{
  typedef LocalfunctionInterface<
      typename SelectDerived<FunctionType, derivative>::E, typename SelectDerived<FunctionType, derivative>::D,
      SelectDerived<FunctionType, derivative>::d, typename SelectDerived<FunctionType, derivative>::R,
      SelectDerived<FunctionType, derivative>::r, SelectDerived<FunctionType, derivative>::rC> BaseType;

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

  virtual size_t order() const override final
  {
    return Select::order(func_local_->order());
  }

  virtual void evaluate(const DomainType& xx, RangeType& ret) const override final
  {
    Select::evaluate(*func_local_, xx, ret);
  }

  virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const override final
  {
    Select::jacobian(*func_local_, xx, ret);
  }

private:
  const std::unique_ptr<const typename FunctionType::LocalfunctionType> func_local_;
}; // class DerivedLocalFunction


template <class FunctionType, Derivative derivative>
class Derived
    : public LocalizableFunctionInterface<
          typename SelectDerived<FunctionType, derivative>::E, typename SelectDerived<FunctionType, derivative>::D,
          SelectDerived<FunctionType, derivative>::d, typename SelectDerived<FunctionType, derivative>::R,
          SelectDerived<FunctionType, derivative>::r, SelectDerived<FunctionType, derivative>::rC>
{
  typedef LocalizableFunctionInterface<
      typename SelectDerived<FunctionType, derivative>::E, typename SelectDerived<FunctionType, derivative>::D,
      SelectDerived<FunctionType, derivative>::d, typename SelectDerived<FunctionType, derivative>::R,
      SelectDerived<FunctionType, derivative>::r, SelectDerived<FunctionType, derivative>::rC> BaseType;
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
    , name_(nm.empty()
                ? SelectDerived<FunctionType, derivative>::type() + " of '" + func_->storage_access().name() + "'"
                : nm)
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
    return DSC::make_unique<RealLocalFunctionType>(func_->storage_access(), entity);
  } // ... local_function(...)

  virtual ThisType* copy() const
  {
    DUNE_THROW(NotImplemented, "Are you kidding me?");
  }

  virtual std::string type() const override final
  {
    return SelectDerived<FunctionType, derivative>::type() + " of '" + func_->storage_access().type() + "'";
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
class Divergence : public internal::Derived<FunctionType, internal::Derivative::divergence>
{
  typedef internal::Derived<FunctionType, internal::Derivative::divergence> BaseType;

public:
  template <class... Args>
  Divergence(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
}; // class Divergence


template <class T, class... Args>
std::shared_ptr<Divergence<T>> make_divergence(const T& func, Args&&... args)
{
  return std::make_shared<Divergence<T>>(func, std::forward<Args>(args)...);
}

template <class T, class... Args>
std::shared_ptr<Divergence<T>> make_divergence(std::shared_ptr<T> func, Args&&... args)
{
  return std::make_shared<Divergence<T>>(func, std::forward<Args>(args)...);
}


} // namespace Functions
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTIONS_DERIVED_HH
