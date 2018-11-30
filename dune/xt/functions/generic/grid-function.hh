// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2018)
//   TiKeil          (2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_FUNCTIONS_GENERIC_GRID_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_GENERIC_GRID_FUNCTION_HH

#include <functional>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/type_traits.hh>

namespace Dune {
namespace XT {
namespace Functions {


/**
 * \brief A function given by a lambda expression or std::function which is evaluated locally on each element.
 *
 *        To model the function f(x) = x^p with a variable exponent p, use as
 * \code
GenericGridFunction<...> f(
                           [](const typename F::ElementType& element,
                              const typename F::DomainType& point_in_local_coordinates,
                              const XT::Common::Parameter&  param) {
                             typename F::RangeType ret(std::pow(point_in_local_coordinates[0],
param.get("power").at(0)));
                             return ret;
                           },
                           integration_order,
                           XT::Common::ParameterType("power", 1),
                           "x_power_p");
\endcode
 *        The XT::Common::ParameterType provided on construction ensures that the XT::Common::Parameter param which is
 *        passed on to the generic function is of correct type.
 * \note  The Localfunction does not implement derivative.
 */
template <class E, size_t r, size_t rC = 1, class R = double>
class GenericGridFunction : public GridFunctionInterface<E, r, rC, R>
{
  using BaseType = GridFunctionInterface<E, r, rC, R>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::LocalFunctionType;

private:
  class LocalGenericGridFunction : public ElementFunctionInterface<E, r, rC, R>
  {
    using BaseType = ElementFunctionInterface<E, r, rC, R>;

  public:
    using typename BaseType::ElementType;
    using typename BaseType::DomainType;
    using typename BaseType::RangeType;
    using typename BaseType::DerivativeRangeType;
    using typename BaseType::RangeReturnType;
    using typename BaseType::DerivativeRangeReturnType;

    using BaseType::d;

    using GenericEvaluateFunctionType = std::function<RangeReturnType(const DomainType&, const Common::Parameter&)>;
    using GenericPostBindFunctionType = std::function<void(const ElementType&)>;
    using GenericOrderFunctionType = std::function<int(const Common::Parameter&)>;
    using GenericJacobianFunctionType =
        std::function<DerivativeRangeReturnType(const ElementType&, const DomainType&, const XT::Common::Parameter&)>;
    using GenericDerivativeFunctionType = std::function<DerivativeRangeReturnType(
        const ElementType&, const std::array<size_t, d>&, const DomainType&, const XT::Common::Parameter&)>;

    LocalGenericGridFunction(const GenericOrderFunctionType& order_func,
                             const GenericPostBindFunctionType& post_bind_func,
                             const GenericEvaluateFunctionType& evalaute_func,
                             const Common::ParameterType& param_type,
                             const GenericJacobianFunctionType& jacobian_func,
                             const GenericDerivativeFunctionType& derivative_func)
      : BaseType()
      , order_(order_func)
      , post_bind_(post_bind_func)
      , evaluate_(evalaute_func)
      , param_type_(param_type)
      , jacobian_(jacobian_func)
      , derivative_(derivative_func)
    {
    }

  protected:
    void post_bind(const ElementType& element) override final
    {
      post_bind_(element);
    }

  public:
    int order(const XT::Common::Parameter& param = {}) const override final
    {
      auto parsed_param = this->parse_parameter(param);
      return order_(parsed_param);
    }

    RangeReturnType evaluate(const DomainType& point_in_local_coordinates,
                             const Common::Parameter& param = {}) const override final
    {
      auto parsed_param = this->parse_parameter(param);
      return evaluate_(point_in_local_coordinates, parsed_param);
    }

    DerivativeRangeReturnType jacobian(const DomainType& point_in_local_coordinates,
                                       const Common::Parameter& param = {}) const override final
    {
      auto parsed_param = this->parse_parameter(param);
      return jacobian_(this->element(), point_in_local_coordinates, parsed_param);
    }

    DerivativeRangeReturnType derivative(const std::array<size_t, d>& alpha,
                                         const DomainType& point_in_local_coordinates,
                                         const Common::Parameter& param = {}) const override final
    {
      auto parsed_param = this->parse_parameter(param);
      return derivative_(this->element(), alpha, point_in_local_coordinates, parsed_param);
    }


    virtual const Common::ParameterType& parameter_type() const override final
    {
      return param_type_;
    }

  private:
    const GenericOrderFunctionType& order_;
    const GenericPostBindFunctionType& post_bind_;
    const GenericEvaluateFunctionType& evaluate_;
    const Common::ParameterType& param_type_;
    const GenericJacobianFunctionType& jacobian_;
    const GenericDerivativeFunctionType& derivative_;
  }; // class LocalGenericGridFunction

public:
  using BaseType::d;
  using DomainType = typename LocalGenericGridFunction::DomainType;
  using RangeType = typename LocalGenericGridFunction::RangeType;
  using DerivativeRangeType = typename LocalGenericGridFunction::DerivativeRangeType;
  using RangeReturnType = typename LocalGenericGridFunction::RangeReturnType;
  using DerivativeRangeReturnType = typename LocalGenericGridFunction::DerivativeRangeReturnType;

  // we do not use the typedef from LocalGenericGridFunction here to document the type of the generic functions
  using GenericOrderFunctionType = std::function<int(const Common::Parameter&)>;
  using GenericPostBindFunctionType = std::function<void(const ElementType&)>;
  using GenericEvaluateFunctionType = std::function<RangeReturnType(const DomainType&, const Common::Parameter&)>;
  using GenericJacobianFunctionType =
      std::function<DerivativeRangeReturnType(const ElementType&, const DomainType&, const XT::Common::Parameter&)>;
  using GenericDerivativeFunctionType = std::function<DerivativeRangeReturnType(
      const ElementType&, const std::array<size_t, d>&, const DomainType&, const XT::Common::Parameter&)>;

  GenericGridFunction(const int ord,
                      GenericPostBindFunctionType post_bind_func = default_post_bind_function(),
                      GenericEvaluateFunctionType evaluate_func = default_evaluate_function(),
                      const Common::ParameterType& param_type = Common::ParameterType(),
                      const std::string nm = "GenericGridFunction",
                      GenericJacobianFunctionType jacobian_func = default_jacobian_function(),
                      GenericDerivativeFunctionType derivative_func = default_derivative_function())
    : order_(default_order_lambda(ord))
    , post_bind_(post_bind_func)
    , evaluate_(evaluate_func)
    , param_type_(param_type)
    , name_(nm)
    , jacobian_(jacobian_func)
    , derivative_(derivative_func)
  {
  }

  GenericGridFunction(GenericOrderFunctionType order_func,
                      GenericPostBindFunctionType post_bind_func = default_post_bind_function(),
                      GenericEvaluateFunctionType evaluate_func = default_evaluate_function(),
                      const Common::ParameterType& param_type = Common::ParameterType(),
                      const std::string nm = "GenericGridFunction",
                      GenericJacobianFunctionType jacobian_func = default_jacobian_function(),
                      GenericDerivativeFunctionType derivative_func = default_derivative_function())
    : order_(order_func)
    , post_bind_(post_bind_func)
    , evaluate_(evaluate_func)
    , param_type_(param_type)
    , name_(nm)
    , jacobian_(jacobian_func)
    , derivative_(derivative_func)
  {
  }

  const Common::ParameterType& parameter_type() const override final
  {
    return param_type_;
  }

  std::string name() const override final
  {
    return name_;
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<LocalGenericGridFunction>(
        order_, post_bind_, evaluate_, param_type_, jacobian_, derivative_);
  }

  /**
   * \}
   * \name ´´These methods may be used to provide defaults on construction.''
   * \{
   */

  static GenericOrderFunctionType default_order_lambda(const int ord)
  {
    return [=](const Common::Parameter& /*param*/ = {}) { return ord; };
  }

  static GenericPostBindFunctionType default_post_bind_function()
  {
    return [](const ElementType& /*element*/) {};
  }

  static GenericEvaluateFunctionType default_evaluate_function()
  {
    return [](const DomainType& /* point_in_local_coordinates*/, const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This GenericGridFunction does not provide evaluations, provide an evaluate_lambda on construction!");
      return RangeType();
    };
  }

  static GenericJacobianFunctionType default_jacobian_function()
  {
    return [](const DomainType& /* point_in_local_coordinates*/, const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This GenericGridFunction does not provide jacobian evaluations, provide a "
                 "jacobian_lambda on construction!");
      return DerivativeRangeReturnType();
    };
  }

  static GenericDerivativeFunctionType default_derivative_function()
  {
    return [](const std::array<size_t, d>& /*alpha*/,
              const DomainType& /* point_in_local_coordinates*/,
              const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This GenericGridFunction does not provide derivative evaluations, provide a "
                 "derivative_lambda on construction!");
      return DerivativeRangeReturnType();
    };
  }

  /**
   * \}
   */

private:
  const GenericOrderFunctionType order_;
  const GenericPostBindFunctionType post_bind_;
  const GenericEvaluateFunctionType evaluate_;
  const Common::ParameterType param_type_;
  const std::string name_;
  GenericJacobianFunctionType jacobian_;
  GenericDerivativeFunctionType derivative_;
}; // class GenericGridFunction


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_LAMBDA_LOCAL_FUNCTION_HH
