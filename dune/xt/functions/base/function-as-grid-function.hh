// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Tobias Leibner  (2017)

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
  using ThisType = FunctionAsGridFunctionWrapper<E, r, rC, R>;

public:
  using BaseType::d;
  using typename BaseType::ElementType;
  using typename BaseType::LocalFunctionType;
  using FunctionType = FunctionInterface<d, r, rC, R>;

  FunctionAsGridFunctionWrapper(const FunctionType& function)
    : function_storage_(function)
  {
  }

  FunctionAsGridFunctionWrapper(const FunctionType*&& function_ptr)
    : function_storage_(std::move(function_ptr))
  {
  }

  FunctionAsGridFunctionWrapper(std::unique_ptr<FunctionType>&& function_ptr)
    : function_storage_(std::move(function_ptr))
  {
  }

  /**
   * \name ´´These methods are required by GridFunctionInterface.''
   * \{
   **/

  std::unique_ptr<LocalFunctionType> local_function(const ElementType& element) const override final
  {
    return std::make_unique<LocalFunction>(function_storage_.access(), element);
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<LocalFunction>(function_storage_.access());
  }

  /**
   * \}
   * \name ´´These methods are optionally required by GridFunctionInterface.''
   * \{
   **/

  std::string type() const override final
  {
    return function_storage_.access().type();
  }

  const std::string name() const override final
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
    using typename BaseType::DomainType;
    using typename BaseType::RangeReturnType;
    using typename BaseType::DerivativeRangeReturnType;

    LocalFunction(const FunctionType& function, const ElementType& el)
      : BaseType(el)
      , function_(function)
      , geometry_(nullptr) // <- happens in post_bind
    {
      post_bind(el);
    }

    LocalFunction(const FunctionType& function)
      : BaseType()
      , function_(function)
      , geometry_(nullptr) // <- happens in post_bind
    {
    }

    using BaseType::evaluate;
    using BaseType::jacobian;
    using BaseType::derivative;

    void post_bind(const ElementType& el) override final
    {
      geometry_ = std::make_unique<GeometryType>(el.geometry());
    }

    int order(const Common::Parameter& param = {}) const override final
    {
      if (!geometry_)
        DUNE_THROW(Exceptions::not_bound_to_an_element_yet, "");
      return function_.order(param);
    }

    RangeReturnType evaluate(const DomainType& point_in_reference_element,
                             const Common::Parameter& param = {}) const override final
    {
      if (!geometry_)
        DUNE_THROW(Exceptions::not_bound_to_an_element_yet, "");
      this->assert_inside_reference_element(point_in_reference_element);
      return function_.evaluate(geometry_->global(point_in_reference_element), param);
    }

    DerivativeRangeReturnType jacobian(const DomainType& point_in_reference_element,
                                       const Common::Parameter& param = {}) const override final
    {
      if (!geometry_)
        DUNE_THROW(Exceptions::not_bound_to_an_element_yet, "");
      this->assert_inside_reference_element(point_in_reference_element);
      return function_.jacobian(geometry_->global(point_in_reference_element), param);
    }

    DerivativeRangeReturnType derivative(const std::array<size_t, d>& alpha,
                                         const DomainType& point_in_reference_element,
                                         const Common::Parameter& param = {}) const override final
    {
      if (!geometry_)
        DUNE_THROW(Exceptions::not_bound_to_an_element_yet, "");
      this->assert_inside_reference_element(point_in_reference_element);
      return function_.derivative(alpha, geometry_->global(point_in_reference_element), param);
    }

  private:
    const FunctionType& function_;
    std::unique_ptr<GeometryType> geometry_;
  }; // class LocalFunction

  const XT::Common::ConstStorageProvider<FunctionType> function_storage_;
}; // class FunctionAsGridFunctionWrapper


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_BASE_SMOOTH_LOCALIZABLE_FUNCTION_HH
