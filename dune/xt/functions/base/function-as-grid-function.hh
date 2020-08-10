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

namespace Dune {
namespace XT {
namespace Functions {


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
    : function_storage_(function)
  {}

  FunctionAsGridFunctionWrapper(FunctionType*&& function_ptr)
    : function_storage_(std::move(function_ptr))
  {}

  FunctionAsGridFunctionWrapper(std::unique_ptr<FunctionType>&& function_ptr)
    : function_storage_(std::move(function_ptr))
  {}

  FunctionAsGridFunctionWrapper(const ThisType&) = default;
  FunctionAsGridFunctionWrapper(ThisType&&) = default;

  /**
   * \name ´´This method is required by GridFunctionInterface.''
   * \{
   **/

  std::unique_ptr<LocalFunctionType> local_function() const final
  {
    return std::make_unique<LocalFunction>(function_storage_.access());
  }

  /**
   * \}
   * \name ´´These methods are optionally required by GridFunctionInterface.''
   * \{
   **/

  std::string name() const final
  {
    return function_storage_.access().name();
  }

  /**
   * \}
   **/

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
      , function_(function)
      , geometry_(nullptr)
    {}

  protected:
    void post_bind(const ElementType& el) final
    {
      geometry_ = std::make_unique<GeometryType>(el.geometry());
    }

  public:
    int order(const Common::Parameter& param = {}) const final
    {
      DUNE_THROW_IF(!(geometry_), Exceptions::not_bound_to_an_element_yet, function_.name());
      return function_.order(param);
    }

    using BaseType::evaluate;

    RangeReturnType evaluate(const DomainType& point_in_reference_element,
                             const Common::Parameter& param = {}) const final
    {
      DUNE_THROW_IF(!(geometry_), Exceptions::not_bound_to_an_element_yet, function_.name());
      this->assert_inside_reference_element(point_in_reference_element);
      return function_.evaluate(geometry_->global(point_in_reference_element), param);
    }

    using BaseType::jacobian;

    DerivativeRangeReturnType jacobian(const DomainType& point_in_reference_element,
                                       const Common::Parameter& param = {}) const final
    {
      DUNE_THROW_IF(!(geometry_), Exceptions::not_bound_to_an_element_yet, function_.name());
      this->assert_inside_reference_element(point_in_reference_element);
      return function_.jacobian(geometry_->global(point_in_reference_element), param);
    }

    using BaseType::derivative;

    DerivativeRangeReturnType derivative(const std::array<size_t, d>& alpha,
                                         const DomainType& point_in_reference_element,
                                         const Common::Parameter& param = {}) const final
    {
      DUNE_THROW_IF(!(geometry_), Exceptions::not_bound_to_an_element_yet, function_.name());
      this->assert_inside_reference_element(point_in_reference_element);
      return function_.derivative(alpha, geometry_->global(point_in_reference_element), param);
    }

  private:
    const FunctionType& function_;
    std::unique_ptr<GeometryType> geometry_;
  }; // class LocalFunction

  XT::Common::ConstStorageProvider<FunctionType> function_storage_;
}; // class FunctionAsGridFunctionWrapper


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_BASE_SMOOTH_LOCALIZABLE_FUNCTION_HH
