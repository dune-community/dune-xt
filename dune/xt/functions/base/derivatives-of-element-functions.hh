// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2019)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_FUNCTIONS_BASE_DERIVATIVES_OF_ELEMENT_FUNCTIONS_HH
#define DUNE_XT_FUNCTIONS_BASE_DERIVATIVES_OF_ELEMENT_FUNCTIONS_HH

#include <dune/xt/common/memory.hh>

#include <dune/xt/functions/interfaces/element-functions.hh>
#include <dune/xt/functions/type_traits.hh>


namespace Dune::XT::Functions {
namespace internal {


template <class ElementFunctionType, DerivativeType derivative>
class DerivativeElementFunctionHelper
{
  static_assert(is_element_function<ElementFunctionType>::value);

public:
  using E = typename ElementFunctionType::E;
  using R = typename ElementFunctionType::R;

private:
  using D = typename ElementFunctionType::D;
  static const constexpr size_t d = ElementFunctionType::d;

private:
  template <class F>
  class DimSwitch
  {
    template <size_t d, size_t r, size_t rC, DerivativeType der, bool anything = true>
    class Dimension
    {
    public:
      static const bool available = false;
    };

    template <size_t d, bool anything>
    class Dimension<d, d, 1, DerivativeType::divergence, anything>
    {
    public:
      static const bool available = true;
      static const size_t r = 1;
      static const size_t rC = 1;
    };

    template <size_t d, bool anything>
    class Dimension<d, 1, 1, DerivativeType::gradient, anything>
    {
    public:
      static const bool available = true;
      static const size_t r = d;
      static const size_t rC = 1;
    };

  public:
    static const bool available = Dimension<d, F::range_dim, F::range_dim_cols, derivative>::available;
    static const size_t r = Dimension<d, F::range_dim, F::range_dim_cols, derivative>::r;
    static const size_t rC = Dimension<d, F::range_dim, F::range_dim_cols, derivative>::rC;
  }; // class DerivativeElementFunctionHelper

public:
  static const bool available = DimSwitch<ElementFunctionType>::available;
  static const size_t r = DimSwitch<ElementFunctionType>::r;
  static const size_t rC = DimSwitch<ElementFunctionType>::rC;

  using DomainType = typename ElementFunctionInterface<E, r, rC, R>::DomainType;
  using RangeReturnType = typename RangeTypeSelector<R, r, rC>::return_type;
  using DerivativeRangeReturnType = typename DerivativeRangeTypeSelector<d, R, r, rC>::return_type;

private:
  template <DerivativeType der, bool anything = true>
  class DerivativeTypeSwitch
  {
    static_assert(!anything, "Nothing available for this derivative!");
  };

  template <bool anything>
  class DerivativeTypeSwitch<DerivativeType::divergence, anything>
  {
  public:
    static int order(const int ord)
    {
      return std::max(ord - 1, 0);
    }

    static RangeReturnType
    evaluate(const ElementFunctionType& func, const DomainType& xx, const Common::Parameter& param)
    {
      auto grads = func.jacobian(xx, param);
      RangeReturnType ret = 0;
      for (size_t dd = 0; dd < d; ++dd)
        ret += grads[dd][dd];
      return ret;
    }

    static DerivativeRangeReturnType
    jacobian(const ElementFunctionType& /*func*/, const DomainType& /*xx*/, const Common::Parameter& /*param*/)
    {
      DUNE_THROW(NotImplemented, "for divergence!");
      return DerivativeRangeReturnType();
    }
  }; // class DerivativeTypeSwitch< ..., divergence >

  template <bool anything>
  class DerivativeTypeSwitch<DerivativeType::gradient, anything>
  {
  public:
    static int order(const int ord)
    {
      return std::max(ord - 1, 0);
    }

    static RangeReturnType
    evaluate(const ElementFunctionType& func, const DomainType& xx, const Common::Parameter& param)
    {
      return func.jacobian(xx, param)[0];
    }

    static DerivativeRangeReturnType
    jacobian(const ElementFunctionType& /*func*/, const DomainType& /*xx*/, const Common::Parameter& /*param*/)
    {
      DUNE_THROW(NotImplemented, "for gradient!");
      return DerivativeRangeReturnType();
    }
  }; // class DerivativeTypeSwitch< ..., divergence >

public:
  static int order(const int ord)
  {
    return DerivativeTypeSwitch<derivative>::order(ord);
  }

  static RangeReturnType evaluate(const ElementFunctionType& func, const DomainType& xx, const Common::Parameter& param)
  {
    return DerivativeTypeSwitch<derivative>::evaluate(func, xx, param);
  }

  static DerivativeRangeReturnType
  jacobian(const ElementFunctionType& func, const DomainType& xx, const Common::Parameter& param)
  {
    return DerivativeTypeSwitch<derivative>::jacobian(func, xx, param);
  }
}; // class DerivativeElementFunctionHelper


} // namespace internal


template <class ElementFunctionType, DerivativeType derivative>
class DerivativeElementFunction
  : public ElementFunctionInterface<
        typename internal::DerivativeElementFunctionHelper<ElementFunctionType, derivative>::E,
        internal::DerivativeElementFunctionHelper<ElementFunctionType, derivative>::r,
        internal::DerivativeElementFunctionHelper<ElementFunctionType, derivative>::rC,
        typename internal::DerivativeElementFunctionHelper<ElementFunctionType, derivative>::R>
{
  using BaseType =
      ElementFunctionInterface<typename internal::DerivativeElementFunctionHelper<ElementFunctionType, derivative>::E,
                               internal::DerivativeElementFunctionHelper<ElementFunctionType, derivative>::r,
                               internal::DerivativeElementFunctionHelper<ElementFunctionType, derivative>::rC,
                               typename internal::DerivativeElementFunctionHelper<ElementFunctionType, derivative>::R>;

  using Select = internal::DerivativeElementFunctionHelper<ElementFunctionType, derivative>;

public:
  using typename BaseType::DerivativeRangeReturnType;
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::RangeReturnType;

  DerivativeElementFunction(ElementFunctionType& func)
    : BaseType()
    , func_(func)
    , do_post_bind_(true)
  {}

  DerivativeElementFunction(std::shared_ptr<ElementFunctionType> func)
    : BaseType()
    , func_(func)
    , do_post_bind_(true)
  {}

  DerivativeElementFunction(std::unique_ptr<ElementFunctionType>&& func)
    : BaseType()
    , func_(std::move(func))
    , do_post_bind_(true)
  {}

protected:
  void post_bind(const ElementType& element)
  {
    if (do_post_bind_)
      func_.access().bind(element);
  }

public:
  int order(const XT::Common::Parameter& param = {}) const final
  {
    return Select::order(func_.access().order(param));
  }

  RangeReturnType evaluate(const DomainType& xx, const Common::Parameter& param = {}) const final
  {
    return Select::evaluate(func_.access(), xx, param);
  }

  DerivativeRangeReturnType jacobian(const DomainType& xx, const Common::Parameter& param = {}) const final
  {
    return Select::jacobian(func_.access(), xx, param);
  }

private:
  void bind_if_argument_was_bound()
  {
    bool was_bound = false;
    try {
      func_.access().element();
      was_bound = true;
    } catch (const XT::Grid::Exceptions::not_bound_to_an_element_yet&) {
      was_bound = false;
    }
    if (was_bound) {
      do_post_bind_ = false;
      this->bind(func_.access().element());
      do_post_bind_ = true;
    }
  } // ... bind_if_argument_was_bound(...)

  XT::Common::StorageProvider<ElementFunctionType> func_;
  bool do_post_bind_;
}; // class DerivativeElementFunction


/**
 * \brief Element function representing the divergence of an element function.
 *
 * \sa DerivativeElementFunction
 */
template <class ElementFunctionType>
class DivergenceElementFunction : public DerivativeElementFunction<ElementFunctionType, DerivativeType::divergence>
{
  using BaseType = DerivativeElementFunction<ElementFunctionType, DerivativeType::divergence>;

public:
  template <class... Args>
  explicit DivergenceElementFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
}; // class DivergenceElementFunction


/**
 * \brief Element function representing the gradient of an element function.
 *
 * \sa DerivativeElementFunction
 */
template <class ElementFunctionType>
class GradientElementFunction : public DerivativeElementFunction<ElementFunctionType, DerivativeType::gradient>
{
  using BaseType = DerivativeElementFunction<ElementFunctionType, DerivativeType::gradient>;

public:
  template <class... Args>
  explicit GradientElementFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
}; // class GradientElementFunction


} // namespace Dune::XT::Functions


#endif // DUNE_XT_FUNCTIONS_DERIVED_HH
