// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2018)
//   Tim Keil        (2018)
//   Tobias Leibner  (2017, 2019)

#ifndef DUNE_XT_FUNCTIONS_GENERIC_GRID_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_GENERIC_GRID_FUNCTION_HH

#include <functional>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/functions/interfaces/flux-function.hh>
#include <dune/xt/functions/type_traits.hh>

namespace Dune {
namespace XT {
namespace Functions {


/**
 * \brief A function given by a lambda expression or std::function which is evaluated locally on each element.
 *
 *        To model the function f(x, u) = x^p * u with a variable exponent p, use as
 * \code
GenericFluxFunction<...> f(
                           [](const typename F::ElementType& element,
                              const typename F::DomainType& point_in_local_coordinates,
                              const typename F::StateType& u,
                              const XT::Common::Parameter&  param) {
                             typename F::RangeType ret(std::pow(point_in_local_coordinates[0],
param.get("power").at(0)));
                             return ret * u;
                           },
                           integration_order,
                           XT::Common::ParameterType("power", 1),
                           "x_power_p");
\endcode
 *        The XT::Common::ParameterType provided on construction ensures that the XT::Common::Parameter param which is
 *        passed on to the generic function is of correct type.
 * \note  The Localfunction does not implement derivative.
 */
template <class E, size_t s, size_t r, size_t rC = 1, class R = double>
class GenericFluxFunction : public FluxFunctionInterface<E, s, r, rC, R>
{
  using BaseType = FluxFunctionInterface<E, s, r, rC, R>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::LocalFunctionType;

private:
  class LocalGenericFluxFunction : public ElementFluxFunctionInterface<E, s, r, rC, R>
  {
    using BaseType = ElementFluxFunctionInterface<E, s, r, rC, R>;

  public:
    using typename BaseType::DomainType;
    using typename BaseType::ElementType;
    using typename BaseType::JacobianRangeReturnType;
    using typename BaseType::JacobianRangeType;
    using typename BaseType::RangeReturnType;
    using typename BaseType::RangeType;
    using typename BaseType::StateType;

    using BaseType::d;

    using GenericEvaluateFunctionType =
        std::function<RangeReturnType(const DomainType&, const StateType&, const Common::Parameter&)>;
    using GenericPostBindFunctionType = std::function<void(const ElementType&)>;
    using GenericOrderFunctionType = std::function<int(const Common::Parameter&)>;
    using GenericJacobianFunctionType =
        std::function<JacobianRangeReturnType(const DomainType&, const StateType&, const XT::Common::Parameter&)>;

    LocalGenericFluxFunction(const GenericOrderFunctionType& order_func,
                             const GenericPostBindFunctionType& post_bind_func,
                             const GenericEvaluateFunctionType& evaluate_func,
                             const Common::ParameterType& param_type,
                             const GenericJacobianFunctionType& jacobian_func)
      : BaseType()
      , order_(order_func)
      , post_bind_(post_bind_func)
      , evaluate_(evaluate_func)
      , param_type_(param_type)
      , jacobian_(jacobian_func)
    {}

  protected:
    virtual void post_bind(const ElementType& element) override final
    {
      post_bind_(element);
    }

  public:
    virtual int order(const XT::Common::Parameter& param = {}) const override final
    {
      auto parsed_param = this->parse_parameter(param);
      return order_(parsed_param);
    }

    virtual RangeReturnType evaluate(const DomainType& point_in_local_coordinates,
                                     const StateType& u,
                                     const Common::Parameter& param = {}) const override final
    {
      auto parsed_param = this->parse_parameter(param);
      return evaluate_(point_in_local_coordinates, u, parsed_param);
    }

    virtual JacobianRangeReturnType jacobian(const DomainType& point_in_local_coordinates,
                                             const StateType& u,
                                             const Common::Parameter& param = {}) const override final
    {
      auto parsed_param = this->parse_parameter(param);
      return jacobian_(point_in_local_coordinates, u, parsed_param);
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
  }; // class LocalGenericFluxFunction

public:
  using BaseType::d;
  using DomainType = typename LocalGenericFluxFunction::DomainType;
  using StateType = typename LocalGenericFluxFunction::StateType;
  using RangeType = typename LocalGenericFluxFunction::RangeType;
  using JacobianRangeType = typename LocalGenericFluxFunction::JacobianRangeType;
  using RangeReturnType = typename LocalGenericFluxFunction::RangeReturnType;
  using JacobianRangeReturnType = typename LocalGenericFluxFunction::JacobianRangeReturnType;

  // we do not use the typedef from LocalGenericFluxFunction here to document the type of the generic functions
  using GenericOrderFunctionType = std::function<int(const Common::Parameter&)>;
  using GenericPostBindFunctionType = std::function<void(const ElementType&)>;
  using GenericEvaluateFunctionType =
      std::function<RangeReturnType(const DomainType&, const StateType&, const Common::Parameter&)>;
  using GenericJacobianFunctionType =
      std::function<JacobianRangeReturnType(const DomainType&, const StateType&, const XT::Common::Parameter&)>;

  GenericFluxFunction(const int ord,
                      GenericPostBindFunctionType post_bind_func = default_post_bind_function(),
                      GenericEvaluateFunctionType evaluate_func = default_evaluate_function(),
                      const Common::ParameterType& param_type = Common::ParameterType(),
                      const std::string nm = "GenericFluxFunction",
                      GenericJacobianFunctionType jacobian_func = default_jacobian_function())
    : order_(default_order_lambda(ord))
    , post_bind_(post_bind_func)
    , evaluate_(evaluate_func)
    , param_type_(param_type)
    , name_(nm)
    , jacobian_(jacobian_func)
  {}

  GenericFluxFunction(GenericOrderFunctionType order_func,
                      GenericPostBindFunctionType post_bind_func = default_post_bind_function(),
                      GenericEvaluateFunctionType evaluate_func = default_evaluate_function(),
                      const Common::ParameterType& param_type = Common::ParameterType(),
                      const std::string nm = "GenericFluxFunction",
                      GenericJacobianFunctionType jacobian_func = default_jacobian_function())
    : order_(order_func)
    , post_bind_(post_bind_func)
    , evaluate_(evaluate_func)
    , param_type_(param_type)
    , name_(nm)
    , jacobian_(jacobian_func)
  {}

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
    return std::make_unique<LocalGenericFluxFunction>(order_, post_bind_, evaluate_, param_type_, jacobian_);
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
    return [](const DomainType& /* point_in_local_coordinates*/,
              const StateType& /*u*/,
              const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This GenericFluxFunction does not provide evaluations, provide an evaluate_lambda on construction!");
      return RangeType();
    };
  }

  static GenericJacobianFunctionType default_jacobian_function()
  {
    return [](const DomainType& /* point_in_local_coordinates*/,
              const StateType& /*u*/,
              const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This GenericFluxFunction does not provide jacobian evaluations, provide a "
                 "jacobian_lambda on construction!");
      return JacobianRangeReturnType();
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
}; // class GenericFluxFunction


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_LAMBDA_LOCAL_FUNCTION_HH
