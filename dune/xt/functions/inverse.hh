// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2019)
//   Tobias Leibner  (2019 - 2020)

#ifndef DUNE_XT_FUNCTIONS_INVERSE_HH
#define DUNE_XT_FUNCTIONS_INVERSE_HH

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/la/matrix-inverter.hh>
#include <dune/xt/functions/interfaces/element-functions.hh>
#include <dune/xt/functions/exceptions.hh>
#include <dune/xt/functions/type_traits.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace internal {


template <class FunctionType>
class InverseFunctionHelper
{
  static_assert(is_element_function<FunctionType>::value || is_function<FunctionType>::value
                    || is_grid_function<FunctionType>::value,
                "");

public:
  static const constexpr size_t d = FunctionType::d;
  static const constexpr size_t r = FunctionType::r;
  static const constexpr size_t rC = FunctionType::rC;
  using R = typename FunctionType::R;
  using DomainType = Dune::FieldVector<double, d>;
  using RangeReturnType = typename RangeTypeSelector<R, r, rC>::return_type;

public:
  static const constexpr bool available = (FunctionType::rC == FunctionType::r);

  static RangeReturnType compute(const FunctionType& func, const DomainType& xx, const XT::Common::Parameter& param)
  {
    if constexpr (FunctionType::rC == 1 && FunctionType::r == 1) {
      auto value_to_invert = func.evaluate(xx, param);
      DUNE_THROW_IF(XT::Common::FloatCmp::eq(value_to_invert, 0.),
                    Exceptions::wrong_input_given,
                    "Scalar function value was not invertible!\n\nvalue_to_invert = " << value_to_invert);
      return 1. / value_to_invert;
    } else if constexpr (available) {
      auto matrix_to_invert = func.evaluate(xx, param);
      RangeReturnType inverse_matrix;
      try {
        inverse_matrix = XT::LA::invert_matrix(matrix_to_invert);
      } catch (const XT::LA::Exceptions::matrix_invert_failed& ee) {
        DUNE_THROW(Exceptions::wrong_input_given,
                   "Matrix-valued function value was not invertible!\n\nmatrix_to_invert = "
                       << matrix_to_invert << "\n\nThis was the original error: " << ee.what());
      }
      return inverse_matrix;
    } else {
      static_assert(AlwaysFalse<FunctionType>::value, "compute not available");
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
  using BaseType = ElementFunctionInterface<typename ElementFunctionType::E,
                                            internal::InverseFunctionHelper<ElementFunctionType>::r,
                                            internal::InverseFunctionHelper<ElementFunctionType>::rC,
                                            typename internal::InverseFunctionHelper<ElementFunctionType>::R>;

  using Helper = internal::InverseFunctionHelper<ElementFunctionType>;

public:
  using typename BaseType::DomainType;
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
  void post_bind(const ElementType& element)
  {
    func_.access().bind(element);
  }

public:
  int order(const XT::Common::Parameter& /*param*/ = {}) const final
  {
    return order_;
  }

  RangeReturnType evaluate(const DomainType& xx, const Common::Parameter& param = {}) const final
  {
    return Helper::compute(func_.access(), xx, param);
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
  using BaseType = FunctionInterface<FunctionType::d,
                                     internal::InverseFunctionHelper<FunctionType>::r,
                                     internal::InverseFunctionHelper<FunctionType>::rC,
                                     typename internal::InverseFunctionHelper<FunctionType>::R>;

  using Helper = internal::InverseFunctionHelper<FunctionType>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeReturnType;

  InverseFunction(FunctionType& func, const int ord)
    : func_(func)
    , order_(ord)
  {}

  InverseFunction(std::shared_ptr<FunctionType> func, const int ord)
    : func_(func)
    , order_(ord)
  {}

  InverseFunction(std::unique_ptr<FunctionType>&& func, const int ord)
    : func_(std::move(func))
    , order_(ord)
  {}

  int order(const XT::Common::Parameter& /*param*/ = {}) const final
  {
    return order_;
  }

  RangeReturnType evaluate(const DomainType& xx, const Common::Parameter& param = {}) const final
  {
    return Helper::compute(func_.access(), xx, param);
  }

private:
  XT::Common::StorageProvider<FunctionType> func_;
  const int order_;
}; // class InverseFunction


template <class GridFunctionType>
class InverseGridFunction
  : public GridFunctionInterface<typename GridFunctionType::E,
                                 internal::InverseFunctionHelper<GridFunctionType>::r,
                                 internal::InverseFunctionHelper<GridFunctionType>::rC,
                                 typename internal::InverseFunctionHelper<GridFunctionType>::R>
{
  using BaseType = GridFunctionInterface<typename GridFunctionType::E,
                                         internal::InverseFunctionHelper<GridFunctionType>::r,
                                         internal::InverseFunctionHelper<GridFunctionType>::rC,
                                         typename internal::InverseFunctionHelper<GridFunctionType>::R>;

  using Helper = internal::InverseFunctionHelper<GridFunctionType>;

public:
  using typename BaseType::LocalFunctionType;

  InverseGridFunction(const GridFunctionType& func, const int ord)
    : func_(func)
    , order_(ord)
  {}

  InverseGridFunction(std::shared_ptr<const GridFunctionType> func, const int ord)
    : func_(func)
    , order_(ord)
  {}

  InverseGridFunction(std::unique_ptr<const GridFunctionType>&& func, const int ord)
    : func_(std::move(func))
    , order_(ord)
  {}

  std::unique_ptr<LocalFunctionType> local_function() const final
  {
    using LocalFunction = InverseElementFunction<typename GridFunctionType::LocalFunctionType>;
    return std::unique_ptr<LocalFunction>(new LocalFunction(std::move(func_.access().local_function()), order_));
  }

  std::string name() const final
  {
    auto func_name = func_.access().name();
    return func_name.empty() ? ("dune.xt.functions.inversegridfunction") : ("inverse of " + func_name);
  }

private:
  const XT::Common::ConstStorageProvider<GridFunctionType> func_;
  const int order_;
}; // class InverseGridFunction


template <class E, size_t r, size_t rC, class R>
auto inverse(ElementFunctionInterface<E, r, rC, R>& func, const int order)
{
  if constexpr (internal::InverseFunctionHelper<ElementFunctionInterface<E, r, rC, R>>::available) {
    return InverseElementFunction<ElementFunctionInterface<E, r, rC, R>>(func, order);
  } else {
    static_assert(AlwaysFalse<R>::value, "No inverse implementation available");
  }
}


template <size_t d, size_t r, size_t rC, class R>
auto inverse(const FunctionInterface<d, r, rC, R>& func, const int order)
{
  if constexpr (internal::InverseFunctionHelper<FunctionInterface<d, r, rC, R>>::available) {
    return InverseFunction<FunctionInterface<d, r, rC, R>>(func, order);
  } else {
    static_assert(AlwaysFalse<R>::value, "No inverse implementation available");
  }
}


template <class E, size_t r, size_t rC, class R>
auto inverse(const GridFunctionInterface<E, r, rC, R>& func, const int order)
{
  if constexpr (internal::InverseFunctionHelper<GridFunctionInterface<E, r, rC, R>>::available) {
    return InverseGridFunction<GridFunctionInterface<E, r, rC, R>>(func, order);
  } else {
    static_assert(AlwaysFalse<R>::value, "No inverse implementation available");
  }
}


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_INVERSE_HH
