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

#ifndef DUNE_XT_FUNCTIONS_BASE_DERIVED_HH
#define DUNE_XT_FUNCTIONS_BASE_DERIVED_HH

#include <memory>
#include <type_traits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/memory.hh>

#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/interfaces/function.hh>


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
  static_assert(is_grid_function<FunctionType>::value, "FunctionType has to be a GridFunction!");

public:
  typedef typename FunctionType::ElementType E;
  typedef typename FunctionType::DomainFieldType D;
  static const size_t d = FunctionType::domain_dim;
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
    static const size_t r = Dimension<d, F::range_dim, F::range_dim_cols, derivative>::r;
    static const size_t rC = Dimension<d, F::range_dim, F::range_dim_cols, derivative>::rC;
  }; // class SelectDerived

public:
  static const size_t r = Choose<FunctionType>::r;
  static const size_t rC = Choose<FunctionType>::rC;

  using FunctionLocalFunctionType = typename FunctionType::LocalFunctionType;
  using DomainType = typename ElementFunctionInterface<E, r, rC, R>::DomainType;
  using RangeType = typename ElementFunctionInterface<E, r, rC, R>::RangeType;
  using DerivativeRangeType = typename ElementFunctionInterface<E, r, rC, R>::DerivativeRangeType;

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

    static int order(const size_t ord)
    {
      return boost::numeric_cast<size_t>(std::max(boost::numeric_cast<ssize_t>(ord) - 1, ssize_t(0)));
    }

    static RangeType
    evaluate(const FunctionLocalFunctionType& func_local, const DomainType& xx, const Common::Parameter& param)
    {
      return func_local.jacobian(xx, param);
    } // ... evaluate(...)

    static DerivativeRangeType jacobian(const FunctionLocalFunctionType& /*func_local*/,
                                        const DomainType& /*xx*/,
                                        const Common::Parameter& /*param*/)
    {
      DUNE_THROW(NotImplemented, "for divergence!");
    }
  }; // class Call< ..., divergence >

public:
  static std::string type()
  {
    return Call<derivative>::type();
  }

  static int order(const size_t ord)
  {
    return Call<derivative>::order(ord);
  }

  static RangeType
  evaluate(const FunctionLocalFunctionType& func_local, const DomainType& xx, const Common::Parameter& param)
  {
    return Call<derivative>::evaluate(func_local, xx, param);
  }

  static DerivativeRangeType
  jacobian(const FunctionLocalFunctionType& func_local, const DomainType& xx, const Common::Parameter& param)
  {
    return Call<derivative>::jacobian(func_local, xx, param);
  }
}; // class SelectDerived


template <class FunctionType, Derivative derivative>
class DerivedLocalFunction : public ElementFunctionInterface<typename SelectDerived<FunctionType, derivative>::E,
                                                             SelectDerived<FunctionType, derivative>::r,
                                                             SelectDerived<FunctionType, derivative>::rC,
                                                             typename SelectDerived<FunctionType, derivative>::R>
{
  using BaseType = ElementFunctionInterface<typename SelectDerived<FunctionType, derivative>::E,
                                            SelectDerived<FunctionType, derivative>::r,
                                            SelectDerived<FunctionType, derivative>::rC,
                                            typename SelectDerived<FunctionType, derivative>::R>;

  using Select = SelectDerived<FunctionType, derivative>;

public:
  using ElementType = typename BaseType::ElementType;
  using DomainType = typename BaseType::DomainType;
  using RangeType = typename BaseType::RangeType;
  using DerivativeRangeType = typename BaseType::DerivativeRangeType;

  DerivedLocalFunction(const FunctionType& func)
    : BaseType()
    , func_local_(func.local_function())
  {
  }

  int order(const XT::Common::Parameter& param = {}) const override final
  {
    return Select::order(func_local_->order(param));
  }

  RangeType evaluate(const DomainType& xx, const Common::Parameter& param = {}) const override final
  {
    return Select::evaluate(*func_local_, xx, param);
  }

  DerivativeRangeType jacobian(const DomainType& xx, const Common::Parameter& param = {}) const override final
  {
    return Select::jacobian(*func_local_, xx, param);
  }

protected:
  void post_bind(const ElementType& element)
  {
    func_local_->bind(element);
  }

private:
  std::unique_ptr<typename FunctionType::LocalFunctionType> func_local_;
}; // class DerivedLocalFunction


template <class FunctionType, Derivative derivative>
class Derived : public GridFunctionInterface<typename SelectDerived<FunctionType, derivative>::E,
                                             SelectDerived<FunctionType, derivative>::r,
                                             SelectDerived<FunctionType, derivative>::rC,
                                             typename SelectDerived<FunctionType, derivative>::R>
{
  using BaseType = GridFunctionInterface<typename SelectDerived<FunctionType, derivative>::E,
                                         SelectDerived<FunctionType, derivative>::r,
                                         SelectDerived<FunctionType, derivative>::rC,
                                         typename SelectDerived<FunctionType, derivative>::R>;
  using FunctionStorageType = Common::ConstStorageProvider<FunctionType>;
  using ThisType = Derived<FunctionType, derivative>;

public:
  using ElementType = typename BaseType::ElementType;
  using LocalFunctionType = typename BaseType::LocalFunctionType;

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

  virtual std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    typedef DerivedLocalFunction<FunctionType, derivative> RealLocalFunctionType;
    assert(func_);
    return Common::make_unique<RealLocalFunctionType>(func_->access());
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


#endif // DUNE_XT_FUNCTIONS_DERIVED_HH
