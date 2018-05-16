// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_FUNCTIONS_LAMBDA_LOCAL_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_LAMBDA_LOCAL_FUNCTION_HH

#include <functional>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>
#include <dune/xt/functions/type_traits.hh>

namespace Dune {
namespace XT {
namespace Functions {


/**
 * \brief A function given by a lambda expression which is evaluated locally on each element.
 *
 *        To model the function f(x) = x^p with a variable exponent p, use as
 * \code
LocalLambdaFunction<...> f(
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
 *        The XT::Common::ParameterType provided on construction ensures that the XT::Common::Parameter  param which is
 *        passed on to the lambda is of correct type.
 * \note  The Localfunction does not implement derivative.
 */
template <class E, size_t r, size_t rC = 1, class R = double>
class LocalLambdaFunction : public LocalizableFunctionInterface<E, r, rC, R>
{
  using BaseType = LocalizableFunctionInterface<E, r, rC, R>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::LocalFunctionType;

private:
  class LocalLambdaLocalFunction : public LocalFunctionInterface<E, r, rC, R>
  {
    using BaseType = LocalFunctionInterface<E, r, rC, R>;

  public:
    using typename BaseType::ElementType;
    using typename BaseType::DomainType;
    using typename BaseType::RangeType;
    using typename BaseType::DerivativeRangeType;
    using BaseType::d;

    using EvaluateLambdaType = std::function<RangeType(const DomainType&, const Common::Parameter&)>;
    using PostBindLambdaType = std::function<void(const ElementType&)>;
    using OrderLambdaType = std::function<int(const Common::Parameter&)>;
    using JacobianLambdaType =
        std::function<DerivativeRangeType(const ElementType&, const DomainType&, const XT::Common::Parameter&)>;
    using DerivativeLambdaType = std::function<DerivativeRangeType(
        const ElementType&, const std::array<size_t, d>&, const DomainType&, const XT::Common::Parameter&)>;


    LocalLambdaLocalFunction(const ElementType& ele,
                             const OrderLambdaType& order_lambda,
                             const PostBindLambdaType& post_bind_lambda,
                             const EvaluateLambdaType& lambda,
                             const Common::ParameterType& param_type,
                             const JacobianLambdaType& jacobian_lambda,
                             const DerivativeLambdaType& derivative_lambda)
      : BaseType(ele)
      , order_lambda_(order_lambda)
      , post_bind_lambda_(post_bind_lambda)
      , evaluate_lambda_(lambda)
      , param_type_(param_type)
      , jacobian_lambda_(jacobian_lambda)
      , derivative_lambda_(derivative_lambda)
    {
      post_bind(ele);
    }

    LocalLambdaLocalFunction(const OrderLambdaType& order_lambda,
                             const PostBindLambdaType& post_bind_lambda,
                             const EvaluateLambdaType& lambda,
                             const Common::ParameterType& param_type,
                             const JacobianLambdaType& jacobian_lambda,
                             const DerivativeLambdaType& derivative_lambda)
      : BaseType()
      , order_lambda_(order_lambda)
      , post_bind_lambda_(post_bind_lambda)
      , evaluate_lambda_(lambda)
      , param_type_(param_type)
      , jacobian_lambda_(jacobian_lambda)
      , derivative_lambda_(derivative_lambda)
    {
    }

    void post_bind(const ElementType& element) override final
    {
      post_bind_lambda_(element);
    }

    int order(const XT::Common::Parameter& param = {}) const override final
    {
      auto parsed_param = this->parse_parameter(param);
      return order_lambda_(parsed_param);
    }

    RangeType evaluate(const DomainType& point_in_local_coordinates,
                       const Common::Parameter& param = {}) const override final
    {
      auto parsed_param = this->parse_parameter(param);
      return evaluate_lambda_(point_in_local_coordinates, parsed_param);
    }

    DerivativeRangeType jacobian(const DomainType& point_in_local_coordinates,
                                 const Common::Parameter& param = {}) const override final
    {
      auto parsed_param = this->parse_parameter(param);
      return jacobian_lambda_(this->element(), point_in_local_coordinates, parsed_param);
    }

    DerivativeRangeType derivative(const std::array<size_t, d>& alpha,
                                   const DomainType& point_in_local_coordinates,
                                   const Common::Parameter& param = {}) const override final
    {
      auto parsed_param = this->parse_parameter(param);
      return derivative_lambda_(this->element(), alpha, point_in_local_coordinates, parsed_param);
    }


    virtual const Common::ParameterType& parameter_type() const override final
    {
      return param_type_;
    }

  private:
    const OrderLambdaType& order_lambda_;
    const PostBindLambdaType& post_bind_lambda_;
    const EvaluateLambdaType& evaluate_lambda_;
    const Common::ParameterType param_type_;
    const JacobianLambdaType& jacobian_lambda_;
    const DerivativeLambdaType& derivative_lambda_;
  }; // class LocalLambdaLocalFunction

public:
  using BaseType::d;
  using DomainType = typename LocalLambdaLocalFunction::DomainType;
  using RangeType = typename LocalLambdaLocalFunction::RangeType;
  using DerivativeRangeType = typename LocalLambdaLocalFunction::DerivativeRangeType;
  // we do not use the typedef from LocalLambdaLocalFunction here to document the type of the lambda
  using OrderLambdaType = std::function<int(const Common::Parameter&)>;
  using PostBindLambdaType = std::function<void(const ElementType&)>;
  using EvaluateLambdaType = std::function<RangeType(const DomainType&, const Common::Parameter&)>;
  using JacobianLambdaType =
      std::function<DerivativeRangeType(const ElementType&, const DomainType&, const XT::Common::Parameter&)>;
  using DerivativeLambdaType = std::function<DerivativeRangeType(
      const ElementType&, const std::array<size_t, d>&, const DomainType&, const XT::Common::Parameter&)>;

  LocalLambdaFunction(const int ord,
                      PostBindLambdaType post_bind_lambda = default_post_bind_lambda(),
                      EvaluateLambdaType evaluate_lambda = default_evaluate_lambda(),
                      const Common::ParameterType& param_type = Common::ParameterType(),
                      const std::string nm = "locallambdafunction",
                      JacobianLambdaType jacobian_lambda = default_jacobian_lambda(),
                      DerivativeLambdaType derivative_lambda = default_derivative_lambda())
    : order_lambda_(default_order_lambda(ord))
    , post_bind_lambda_(post_bind_lambda)
    , evaluate_lambda_(evaluate_lambda)
    , param_type_(param_type)
    , name_(nm)
    , jacobian_lambda_(jacobian_lambda)
    , derivative_lambda_(derivative_lambda)
  {
  }

  LocalLambdaFunction(OrderLambdaType order_lambda,
                      PostBindLambdaType post_bind_lambda = default_post_bind_lambda(),
                      EvaluateLambdaType evaluate_lambda = default_evaluate_lambda(),
                      const Common::ParameterType& param_type = Common::ParameterType(),
                      const std::string nm = "locallambdafunction",
                      JacobianLambdaType jacobian_lambda = default_jacobian_lambda(),
                      DerivativeLambdaType derivative_lambda = default_derivative_lambda())
    : order_lambda_(order_lambda)
    , post_bind_lambda_(post_bind_lambda)
    , evaluate_lambda_(evaluate_lambda)
    , param_type_(param_type)
    , name_(nm)
    , jacobian_lambda_(jacobian_lambda)
    , derivative_lambda_(derivative_lambda)
  {
  }

  const Common::ParameterType& parameter_type() const override final
  {
    return param_type_;
  }

  std::string type() const override final
  {
    return "locallambdafunction";
  }

  std::string name() const override final
  {
    return name_;
  }

  std::unique_ptr<LocalFunctionType> local_function(const ElementType& element) const override final
  {
    return std::make_unique<LocalLambdaLocalFunction>(
        element, order_lambda_, post_bind_lambda_, evaluate_lambda_, param_type_, jacobian_lambda_, derivative_lambda_);
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<LocalLambdaLocalFunction>(
        order_lambda_, post_bind_lambda_, evaluate_lambda_, param_type_, jacobian_lambda_, derivative_lambda_);
  }

  /**
   * \}
   * \name ´´These methods may be used to provide defaults on construction.''
   * \{
   */

  static OrderLambdaType default_order_lambda(const int ord)
  {
    return [=](const Common::Parameter& /*param*/ = {}) { return ord; };
  }

  static PostBindLambdaType default_post_bind_lambda()
  {
    return [](const ElementType& /*element*/) {};
  }

  static EvaluateLambdaType default_evaluate_lambda()
  {
    return [](const DomainType& /* point_in_local_coordinates*/, const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This  LocalLambdaFunction does not provide evaluations, provide an evaluate_lambda on construction!");
      return RangeType();
    };
  }

  static JacobianLambdaType default_jacobian_lambda()
  {
    return [](const DomainType& /* point_in_local_coordinates*/, const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This  LocalLambdaFunction does not provide jacobian evaluations, provide a "
                 "jacobian_lambda on construction!");
      return DerivativeRangeType();
    };
  }

  static DerivativeLambdaType default_derivative_lambda()
  {
    return [](const std::array<size_t, d>& /*alpha*/,
              const DomainType& /* point_in_local_coordinates*/,
              const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This  LocalLambdaFunction does not provide derivative evaluations, provide a "
                 "derivative_lambda on construction!");
      return DerivativeRangeType();
    };
  }

  /**
   * \}
   */

private:
  const OrderLambdaType order_lambda_;
  const PostBindLambdaType post_bind_lambda_;
  const EvaluateLambdaType evaluate_lambda_;
  const Common::ParameterType param_type_;
  const std::string name_;
  JacobianLambdaType jacobian_lambda_;
  DerivativeLambdaType derivative_lambda_;
}; // class LocalLambdaFunction


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_LAMBDA_LOCAL_FUNCTION_HH
