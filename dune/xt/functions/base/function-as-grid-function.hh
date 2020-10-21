// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2019)
//   René Fritze     (2017 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2017, 2019 - 2020)

#ifndef DUNE_XT_FUNCTIONS_BASE_FUNCTION_AS_GRID_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_BASE_FUNCTION_AS_GRID_FUNCTION_HH

#include <dune/xt/common/memory.hh>

#include <dune/xt/functions/exceptions.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/interfaces/function.hh>
#include <dune/xt/functions/type_traits.hh>

namespace Dune::XT::Functions {


template <class E, size_t r, size_t rC, class R>
class FunctionAsGridFunctionWrapper : public GridFunctionInterface<E, r, rC, R>
{
  using BaseType = GridFunctionInterface<E, r, rC, R>;
  using ThisType = FunctionAsGridFunctionWrapper;

public:
  using BaseType::d;
  using typename BaseType::ElementType;
  using typename BaseType::LocalFunctionType;
  using FunctionType = FunctionInterface<d, r, rC, R>;

  FunctionAsGridFunctionWrapper(const FunctionType& function)
    : BaseType(function.parameter_type())
    , function_(function.copy_as_function())
  {}

  FunctionAsGridFunctionWrapper(FunctionType*&& function_ptr)
    : BaseType(function_ptr->parameter_type())
    , function_(std::move(function_ptr))
  {}

  FunctionAsGridFunctionWrapper(std::unique_ptr<FunctionType>&& function_ptr)
    : BaseType(function_ptr->parameter_type())
    , function_(std::move(function_ptr))
  {}

  FunctionAsGridFunctionWrapper(const ThisType& other)
    : BaseType(other)
    , function_(other.function_->copy_as_function())
  {}

  FunctionAsGridFunctionWrapper(ThisType&&) = default;


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
    return std::make_unique<LocalFunction>(*function_);
  }

  std::string name() const override final
  {
    return function_->name();
  }

private:
  class LocalFunction : public LocalFunctionType
  {
    using BaseType = LocalFunctionType;
    using GeometryType = typename ElementType::Geometry;

  public:
    using typename BaseType::DerivativeRangeReturnType;
    using typename BaseType::DomainType;
    using typename BaseType::RangeReturnType;

    LocalFunction(const FunctionType& function)
      : BaseType()
      , function_(function.copy_as_function())
      , geometry_(nullptr)
    {}

  protected:
    void post_bind(const ElementType& el) override final
    {
      geometry_ = std::make_unique<GeometryType>(el.geometry());
    }

  public:
    int order(const Common::Parameter& param = {}) const override final
    {
      DUNE_THROW_IF(!(geometry_), Exceptions::not_bound_to_an_element_yet, function_->name());
      return function_->order(param);
    }

    using BaseType::evaluate;

    RangeReturnType evaluate(const DomainType& point_in_reference_element,
                             const Common::Parameter& param = {}) const override final
    {
      DUNE_THROW_IF(!(geometry_), Exceptions::not_bound_to_an_element_yet, function_->name());
      this->assert_inside_reference_element(point_in_reference_element);
      return function_->evaluate(geometry_->global(point_in_reference_element), param);
    }

    using BaseType::jacobian;

    DerivativeRangeReturnType jacobian(const DomainType& point_in_reference_element,
                                       const Common::Parameter& param = {}) const override final
    {
      DUNE_THROW_IF(!(geometry_), Exceptions::not_bound_to_an_element_yet, function_->name());
      this->assert_inside_reference_element(point_in_reference_element);
      return function_->jacobian(geometry_->global(point_in_reference_element), param);
    }

    using BaseType::derivative;

    DerivativeRangeReturnType derivative(const std::array<size_t, d>& alpha,
                                         const DomainType& point_in_reference_element,
                                         const Common::Parameter& param = {}) const override final
    {
      DUNE_THROW_IF(!(geometry_), Exceptions::not_bound_to_an_element_yet, function_->name());
      this->assert_inside_reference_element(point_in_reference_element);
      return function_->derivative(alpha, geometry_->global(point_in_reference_element), param);
    }

  private:
    std::unique_ptr<FunctionType> function_;
    std::unique_ptr<GeometryType> geometry_;
  }; // class LocalFunction

  std::unique_ptr<FunctionType> function_;
}; // class FunctionAsGridFunctionWrapper


} // namespace Dune::XT::Functions

#endif // DUNE_XT_FUNCTIONS_BASE_SMOOTH_LOCALIZABLE_FUNCTION_HH
