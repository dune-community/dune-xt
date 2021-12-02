// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018, 2020)
//   René Fritze     (2018 - 2020)
//   Tim Keil        (2018)
//   Tobias Leibner  (2017, 2019 - 2020)

#ifndef DUNE_XT_FUNCTIONS_GENERIC_FLUX_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_GENERIC_FLUX_FUNCTION_HH

#include <functional>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/functions/interfaces/flux-function.hh>
#include <dune/xt/functions/type_traits.hh>
#include <utility>

namespace Dune::XT::Functions {


/**
 * \brief A function given by a lambda expression or std::function which is evaluated locally on each element.
 *
 *        To model the function f(x, u) = x[0]^p * u with a variable exponent p, use as
 * \code
GenericFluxFunction<...> f([](const auto& param){ return param.get("power").at(0)}, // order
                           [](const auto& grid_element){ // do nothing }            // post_bind
                           [](const auto& x, const auto& u, const auto& param)      // evaluate
                             {
                               const auto p = param.get("power").at(0);
                               return u * std::pow(x[0], p);
                             }
                           XT::Common::ParameterType("power", 1),                   // parameter_type
                           "x_power_p");                                            // name
\endcode
 *        The XT::Common::ParameterType provided on construction ensures that the XT::Common::Parameter param which is
 *        passed on to the generic function is of correct type.
 * \note  The Localfunction does not implement derivative.
 */
template <class E, size_t s = 1, size_t r = 1, size_t rC = 1, class R = double>
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
    using typename BaseType::DynamicJacobianRangeType;
    using typename BaseType::DynamicRangeType;
    using typename BaseType::ElementType;
    using typename BaseType::JacobianRangeReturnType;
    using typename BaseType::JacobianRangeType;
    using typename BaseType::RangeReturnType;
    using typename BaseType::RangeType;
    using typename BaseType::StateType;

    using BaseType::d;

    using GenericEvaluateFunctionType =
        std::function<RangeReturnType(const DomainType&, const StateType&, const Common::Parameter&)>;
    using GenericDynamicEvaluateFunctionType =
        std::function<void(const DomainType&, const StateType&, DynamicRangeType&, const Common::Parameter&)>;
    using GenericPostBindFunctionType = std::function<void(const ElementType&)>;
    using GenericOrderFunctionType = std::function<int(const Common::Parameter&)>;
    using GenericJacobianFunctionType =
        std::function<JacobianRangeReturnType(const DomainType&, const StateType&, const XT::Common::Parameter&)>;
    using GenericDynamicJacobianFunctionType = std::function<void(
        const DomainType&, const StateType&, DynamicJacobianRangeType&, const XT::Common::Parameter&)>;

    LocalGenericFluxFunction(const GenericOrderFunctionType& order_func,
                             const GenericPostBindFunctionType& post_bind_func,
                             const GenericEvaluateFunctionType& evaluate_func,
                             const GenericDynamicEvaluateFunctionType& dynamic_evaluate_func,
                             const Common::ParameterType& param_type,
                             const GenericJacobianFunctionType& jacobian_func,
                             const GenericDynamicJacobianFunctionType& dynamic_jacobian_func)
      : BaseType()
      , order_(order_func)
      , post_bind_(post_bind_func)
      , evaluate_(evaluate_func)
      , dynamic_evaluate_(dynamic_evaluate_func)
      , param_type_(param_type)
      , jacobian_(jacobian_func)
      , dynamic_jacobian_(dynamic_jacobian_func)
    {}

  protected:
    void post_bind(const ElementType& element) final
    {
      post_bind_(element);
    }

  public:
    int order(const XT::Common::Parameter& param = {}) const final
    {
      auto parsed_param = this->parse_parameter(param);
      return order_(parsed_param);
    }

    RangeReturnType evaluate(const DomainType& point_in_local_coordinates,
                             const StateType& u,
                             const Common::Parameter& param = {}) const final
    {
      auto parsed_param = this->parse_parameter(param);
      return evaluate_(point_in_local_coordinates, u, parsed_param);
    }

    JacobianRangeReturnType jacobian(const DomainType& point_in_local_coordinates,
                                     const StateType& u,
                                     const Common::Parameter& param = {}) const final
    {
      auto parsed_param = this->parse_parameter(param);
      return jacobian_(point_in_local_coordinates, u, parsed_param);
    }

    void evaluate(const DomainType& point_in_local_coordinates,
                  const StateType& u,
                  DynamicRangeType& ret,
                  const Common::Parameter& param = {}) const final
    {
      auto parsed_param = this->parse_parameter(param);
      dynamic_evaluate_(point_in_local_coordinates, u, ret, parsed_param);
    }

    void jacobian(const DomainType& point_in_local_coordinates,
                  const StateType& u,
                  DynamicJacobianRangeType& ret,
                  const Common::Parameter& param = {}) const final
    {
      auto parsed_param = this->parse_parameter(param);
      return dynamic_jacobian_(point_in_local_coordinates, u, ret, parsed_param);
    }


    const Common::ParameterType& parameter_type() const final
    {
      return param_type_;
    }

  private:
    const GenericOrderFunctionType& order_;
    const GenericPostBindFunctionType& post_bind_;
    const GenericEvaluateFunctionType& evaluate_;
    const GenericDynamicEvaluateFunctionType& dynamic_evaluate_;
    const Common::ParameterType& param_type_;
    const GenericJacobianFunctionType& jacobian_;
    const GenericDynamicJacobianFunctionType& dynamic_jacobian_;
  }; // class LocalGenericFluxFunction

public:
  using BaseType::d;
  using DomainType = typename LocalGenericFluxFunction::DomainType;
  using StateType = typename LocalGenericFluxFunction::StateType;
  using RangeType = typename LocalGenericFluxFunction::RangeType;
  using DynamicRangeType = typename LocalGenericFluxFunction::DynamicRangeType;
  using JacobianRangeType = typename LocalGenericFluxFunction::JacobianRangeType;
  using RangeReturnType = typename LocalGenericFluxFunction::RangeReturnType;
  using JacobianRangeReturnType = typename LocalGenericFluxFunction::JacobianRangeReturnType;
  using DynamicJacobianRangeType = typename LocalGenericFluxFunction::DynamicJacobianRangeType;

  // we do not use the typedef from LocalGenericFluxFunction here to document the type of the generic functions
  using GenericOrderFunctionType = std::function<int(const Common::Parameter&)>;
  using GenericPostBindFunctionType = std::function<void(const ElementType&)>;
  using GenericEvaluateFunctionType =
      std::function<RangeReturnType(const DomainType&, const StateType&, const Common::Parameter&)>;
  using GenericDynamicEvaluateFunctionType =
      std::function<void(const DomainType&, const StateType&, DynamicRangeType&, const Common::Parameter&)>;
  using GenericJacobianFunctionType =
      std::function<JacobianRangeReturnType(const DomainType&, const StateType&, const Common::Parameter&)>;
  using GenericDynamicJacobianFunctionType =
      std::function<void(const DomainType&, const StateType&, DynamicJacobianRangeType&, const Common::Parameter&)>;

  GenericFluxFunction(const int ord,
                      GenericPostBindFunctionType post_bind_func = default_post_bind_function(),
                      GenericEvaluateFunctionType evaluate_func = default_evaluate_function(),
                      const Common::ParameterType& param_type = Common::ParameterType(),
                      const std::string& nm = "GenericFluxFunction",
                      GenericJacobianFunctionType jacobian_func = default_jacobian_function())
    : order_(default_order_lambda(ord))
    , post_bind_(post_bind_func)
    , evaluate_(evaluate_func)
    , dynamic_evaluate_(default_dynamic_evaluate_function())
    , param_type_(param_type)
    , name_(nm)
    , jacobian_(jacobian_func)
    , dynamic_jacobian_(default_dynamic_jacobian_function())
  {}

  GenericFluxFunction(GenericOrderFunctionType order_func,
                      GenericPostBindFunctionType post_bind_func = default_post_bind_function(),
                      GenericEvaluateFunctionType evaluate_func = default_evaluate_function(),
                      const Common::ParameterType& param_type = Common::ParameterType(),
                      const std::string& nm = "GenericFluxFunction",
                      GenericJacobianFunctionType jacobian_func = default_jacobian_function())
    : order_(std::move(order_func))
    , post_bind_(post_bind_func)
    , evaluate_(evaluate_func)
    , dynamic_evaluate_(default_dynamic_evaluate_function())
    , param_type_(param_type)
    , name_(nm)
    , jacobian_(jacobian_func)
    , dynamic_jacobian_(default_dynamic_jacobian_function())
  {}

  GenericFluxFunction(const int ord,
                      GenericPostBindFunctionType post_bind_func,
                      GenericDynamicEvaluateFunctionType evaluate_func,
                      const Common::ParameterType& param_type = Common::ParameterType(),
                      const std::string& nm = "GenericFluxFunction",
                      GenericDynamicJacobianFunctionType jacobian_func = default_dynamic_jacobian_function())
    : order_(default_order_lambda(ord))
    , post_bind_(post_bind_func)
    , evaluate_(evaluate_from_dynamic_evaluate(evaluate_func))
    , dynamic_evaluate_(evaluate_func)
    , param_type_(param_type)
    , name_(nm)
    , jacobian_(jacobian_from_dynamic_jacobian(jacobian_func))
    , dynamic_jacobian_(jacobian_func)
  {}

  GenericFluxFunction(GenericOrderFunctionType order_func,
                      GenericPostBindFunctionType post_bind_func,
                      GenericDynamicEvaluateFunctionType evaluate_func,
                      const Common::ParameterType& param_type = Common::ParameterType(),
                      const std::string& nm = "GenericFluxFunction",
                      GenericDynamicJacobianFunctionType jacobian_func = default_dynamic_jacobian_function())
    : order_(std::move(order_func))
    , post_bind_(post_bind_func)
    , evaluate_(evaluate_from_dynamic_evaluate(evaluate_func))
    , dynamic_evaluate_(evaluate_func)
    , param_type_(param_type)
    , name_(nm)
    , jacobian_(jacobian_from_dynamic_jacobian(jacobian_func))
    , dynamic_jacobian_(jacobian_func)
  {}

  const Common::ParameterType& parameter_type() const final
  {
    return param_type_;
  }

  std::string name() const final
  {
    return name_;
  }

  std::unique_ptr<LocalFunctionType> local_function() const final
  {
    return std::make_unique<LocalGenericFluxFunction>(
        order_, post_bind_, evaluate_, dynamic_evaluate_, param_type_, jacobian_, dynamic_jacobian_);
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

  static GenericDynamicEvaluateFunctionType default_dynamic_evaluate_function()
  {
    return [](const DomainType& /* point_in_local_coordinates*/,
              const StateType& /*u*/,
              DynamicRangeType& /*ret*/,
              const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(
          NotImplemented,
          "This GenericFluxFunction does not provide dynamic evaluations, provide an evaluate_lambda on construction!");
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

  static GenericDynamicJacobianFunctionType default_dynamic_jacobian_function()
  {
    return [](const DomainType& /* point_in_local_coordinates*/,
              const StateType& /*u*/,
              const DynamicJacobianRangeType& /*ret*/,
              const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This GenericFluxFunction does not provide dynamic jacobian evaluations, provide a "
                 "jacobian_lambda on construction!");
    };
  }

  static GenericEvaluateFunctionType
  evaluate_from_dynamic_evaluate(const GenericDynamicEvaluateFunctionType& dynamic_evaluate)
  {
    return [dynamic_evaluate](
               const DomainType& point_in_local_coordinates, const StateType& u, const Common::Parameter& param = {}) {
      DynamicRangeType dynamic_ret;
      BaseType::LocalFunctionType::RangeSelector::ensure_size(dynamic_ret);
      dynamic_evaluate(point_in_local_coordinates, u, dynamic_ret, param);
      RangeReturnType ret = dynamic_ret;
      return ret;
    };
  }

  static GenericJacobianFunctionType
  jacobian_from_dynamic_jacobian(const GenericDynamicJacobianFunctionType& dynamic_jacobian)
  {
    return [dynamic_jacobian](
               const DomainType& point_in_local_coordinates, const StateType& u, const Common::Parameter& param = {}) {
      DynamicJacobianRangeType dynamic_ret;
      BaseType::LocalFunctionType::JacobianRangeSelector::ensure_size(dynamic_ret);
      dynamic_jacobian(point_in_local_coordinates, u, dynamic_ret, param);
      JacobianRangeReturnType ret = dynamic_ret;
      return ret;
    };
  }

  /**
   * \}
   */

private:
  const GenericOrderFunctionType order_;
  const GenericPostBindFunctionType post_bind_;
  const GenericEvaluateFunctionType evaluate_;
  const GenericDynamicEvaluateFunctionType dynamic_evaluate_;
  const Common::ParameterType param_type_;
  const std::string name_;
  GenericJacobianFunctionType jacobian_;
  GenericDynamicJacobianFunctionType dynamic_jacobian_;
}; // class GenericFluxFunction


} // namespace Dune::XT::Functions

#endif // DUNE_XT_FUNCTIONS_GENERIC_FLUX_FUNCTION_HH
