// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019 - 2020)
//   René Fritze     (2019 - 2020)
//   Tobias Leibner  (2019 - 2020)

#ifndef DUNE_XT_FUNCTIONS_INVERSE_HH
#define DUNE_XT_FUNCTIONS_INVERSE_HH

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/la/matrix-inverter.hh>
#include <dune/xt/functions/grid-function.hh>
#include <dune/xt/functions/interfaces/element-functions.hh>
#include <dune/xt/functions/interfaces/function.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/exceptions.hh>
#include <dune/xt/functions/type_traits.hh>

namespace Dune::XT::Functions {
namespace internal {


template <class FunctionType>
class InverseFunctionHelper
{
  static_assert(is_element_function<FunctionType>::value || is_function<FunctionType>::value
                || is_grid_function<FunctionType>::value);

public:
  static constexpr size_t d = FunctionType::d;
  static constexpr size_t r = FunctionType::r;
  static constexpr size_t rC = FunctionType::rC;
  using R = typename FunctionType::R;
  using DomainType = Dune::FieldVector<double, d>;
  using RangeReturnType = typename RangeTypeSelector<R, r, rC>::return_type;


  static constexpr bool available = (FunctionType::rC == FunctionType::r);

  static RangeReturnType compute(const RangeReturnType& value)
  {
    if constexpr (FunctionType::rC == 1 && FunctionType::r == 1) {
      auto value_to_invert = value[0];
      DUNE_THROW_IF(XT::Common::FloatCmp::eq(value_to_invert, 0.),
                    Exceptions::wrong_input_given,
                    "Scalar function value was not invertible!\n\nvalue_to_invert = " << value_to_invert);
      return 1. / value_to_invert;
    } else if constexpr (available) {
      RangeReturnType inverse_matrix;
      try {
        inverse_matrix = XT::LA::invert_matrix(value);
      } catch (const XT::LA::Exceptions::matrix_invert_failed& ee) {
        DUNE_THROW(Exceptions::wrong_input_given,
                   "Matrix-valued function value was not invertible!\n\nmatrix_to_invert = "
                       << value << "\n\nThis was the original error: " << ee.what());
      }
      return inverse_matrix;
    } else {
      static_assert(AlwaysFalse<FunctionType>::value, "Not available for these dimensions!");
    }
  }
}; // class InverseFunctionHelper


} // namespace internal


template <class ElementFunctionType>
class InverseElementFunction
  : public ElementFunctionInterface<typename ElementFunctionType::E,
                                    internal::InverseFunctionHelper<ElementFunctionType>::r,
                                    internal::InverseFunctionHelper<ElementFunctionType>::rC,
                                    typename internal::InverseFunctionHelper<ElementFunctionType>::R>
{
  static_assert(is_element_function<ElementFunctionType>::value);

  using BaseType = ElementFunctionInterface<typename ElementFunctionType::E,
                                            internal::InverseFunctionHelper<ElementFunctionType>::r,
                                            internal::InverseFunctionHelper<ElementFunctionType>::rC,
                                            typename internal::InverseFunctionHelper<ElementFunctionType>::R>;

  using Helper = internal::InverseFunctionHelper<ElementFunctionType>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::ElementType;
  using typename BaseType::RangeReturnType;

  InverseElementFunction(ElementFunctionType& func, const int ord)
    : func_(func)
    , order_(ord)
  {}

  InverseElementFunction(std::shared_ptr<ElementFunctionType> func, const int ord)
    : func_(func)
    , order_(ord)
  {}

  InverseElementFunction(std::unique_ptr<ElementFunctionType>&& func, const int ord)
    : func_(std::move(func))
    , order_(ord)
  {}

protected:
  void post_bind(const ElementType& element) override final
  {
    func_.access().bind(element);
  }

public:
  int order(const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    return order_;
  }

  RangeReturnType evaluate(const DomainType& xx, const Common::Parameter& param = {}) const override final
  {
    return Helper::compute(func_.access().evaluate(xx, param));
  }

private:
  XT::Common::StorageProvider<ElementFunctionType> func_;
  const int order_;
}; // class InverseElementFunction


template <class FunctionType>
class InverseFunction
  : public FunctionInterface<FunctionType::d,
                             internal::InverseFunctionHelper<FunctionType>::r,
                             internal::InverseFunctionHelper<FunctionType>::rC,
                             typename internal::InverseFunctionHelper<FunctionType>::R>
{
  static_assert(is_function<FunctionType>::value);

  using ThisType = InverseFunction;
  using BaseType = FunctionInterface<FunctionType::d,
                                     internal::InverseFunctionHelper<FunctionType>::r,
                                     internal::InverseFunctionHelper<FunctionType>::rC,
                                     typename internal::InverseFunctionHelper<FunctionType>::R>;

  using Helper = internal::InverseFunctionHelper<FunctionType>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeReturnType;

  InverseFunction(const FunctionType& func, const int ord)
    : BaseType(func.parameter_type())
    , func_(func.copy_as_function())
    , order_(ord)
  {}

  InverseFunction(FunctionType&& func, const int ord)
    : BaseType(func->parameter_type())
    , func_(std::move(func))
    , order_(ord)
  {}

  InverseFunction(const ThisType& other)
    : BaseType(other)
    , func_(other.func_->copy_as_function())
    , order_(order)
  {}

  InverseFunction(ThisType&&) = default;

private:
  ThisType* copy_as_function_impl() const override
  {
    return new ThisType(*this);
  }

public:
  std::unique_ptr<ThisType> copy_as_function() const
  {
    return std::unique_ptr<ThisType>(this->copy_as_function_impl());
  }
  int order(const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    return order_;
  }

  RangeReturnType evaluate(const DomainType& xx, const Common::Parameter& param = {}) const override final
  {
    return Helper::compute(func_.evaluate(xx, param));
  }

private:
  const std::unique_ptr<FunctionType> func_;
  const int order_;
}; // class InverseFunction


/**
 * \todo Write custom local function to hold a copy af this!
 */
template <class GridFunctionType>
class InverseGridFunction
  : public GridFunctionInterface<typename GridFunctionType::E,
                                 internal::InverseFunctionHelper<GridFunctionType>::r,
                                 internal::InverseFunctionHelper<GridFunctionType>::rC,
                                 typename internal::InverseFunctionHelper<GridFunctionType>::R>
{
  static_assert(is_grid_function<GridFunctionType>::value);

  using Helper = internal::InverseFunctionHelper<GridFunctionType>;
  static constexpr size_t r_ = GridFunctionType::r;
  static constexpr size_t rC_ = GridFunctionType::rC;

public:
  using ThisType = InverseGridFunction;
  using BaseType = GridFunctionInterface<typename GridFunctionType::E,
                                         internal::InverseFunctionHelper<GridFunctionType>::r,
                                         internal::InverseFunctionHelper<GridFunctionType>::rC,
                                         typename internal::InverseFunctionHelper<GridFunctionType>::R>;
  using typename BaseType::E;
  using typename BaseType::LocalFunctionType;
  using typename BaseType::R;

  InverseGridFunction(GridFunction<E, r_, rC_, R> func, const int ord, const std::string& nm = "")
    : func_(func.copy_as_grid_function())
    , order_(ord)
    , name_(nm.empty() ? ("inverse of " + func_->name()) : nm)
  {}

  InverseGridFunction(GridFunctionInterface<E, r_, rC_, R>&& func, const int ord, const std::string& nm = "")
    : func_(std::move(func))
    , order_(ord)
    , name_(nm.empty() ? ("inverse of " + func_->name()) : nm)
  {}

  InverseGridFunction(const ThisType& other)
    : BaseType(other)
    , func_(other.func_->copy_as_grid_function())
    , order_(other.order_)
  {}

  InverseGridFunction(ThisType&&) = default;


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
  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    using LocalFunction = InverseElementFunction<typename GridFunctionType::LocalFunctionType>;
    return std::make_unique<LocalFunction>(func_->local_function(), order_);
  }

  std::string name() const override final
  {
    return name_;
  }

private:
  const std::unique_ptr<GridFunctionInterface<E, r_, rC_, R>> func_;
  const int order_;
  const std::string name_;
}; // class InverseGridFunction


template <class E, size_t r, size_t rC, class R>
auto inverse(ElementFunctionInterface<E, r, rC, R>& func, const int order)
{
  if constexpr (internal::InverseFunctionHelper<ElementFunctionInterface<E, r, rC, R>>::available) {
    return InverseElementFunction<ElementFunctionInterface<E, r, rC, R>>(func, order);
  } else {
    static_assert(AlwaysFalse<R>::value, "Not available for these dimensions!");
  }
}


template <size_t d, size_t r, size_t rC, class R>
auto inverse(const FunctionInterface<d, r, rC, R>& func, const int order)
{
  if constexpr (internal::InverseFunctionHelper<FunctionInterface<d, r, rC, R>>::available) {
    return InverseFunction<FunctionInterface<d, r, rC, R>>(func, order);
  } else {
    static_assert(AlwaysFalse<R>::value, "Not available for these dimensions!");
  }
}


template <class E, size_t r, size_t rC, class R>
auto inverse(const GridFunctionInterface<E, r, rC, R>& func, const int order)
{
  if constexpr (internal::InverseFunctionHelper<GridFunctionInterface<E, r, rC, R>>::available) {
    return InverseGridFunction<GridFunctionInterface<E, r, rC, R>>(func, order);
  } else {
    static_assert(AlwaysFalse<R>::value, "Not available for these dimensions!");
  }
}


} // namespace Dune::XT::Functions

#endif // DUNE_XT_FUNCTIONS_INVERSE_HH
