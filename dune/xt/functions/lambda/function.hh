// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2013 - 2016, 2018)
//   Sven Kaulmann   (2013)
//   Tobias Leibner  (2014, 2017)

#ifndef DUNE_XT_FUNCTIONS_LAMBDA_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_LAMBDA_FUNCTION_HH

#include <functional>

#include <dune/xt/common/memory.hh>

#include <dune/xt/functions/interfaces/function.hh>

namespace Dune {
namespace XT {
namespace Functions {


/**
 * Smooth function you can pass lambda expressions to that gets evaluated.
 *
 * \example LambdaType lambda(1, [](const auto& x, const auto& param = {}) { return x;});
 */
template <size_t domain_dim, size_t range_dim = 1, size_t range_dim_cols = 1, class RangeField = double>
class LambdaFunction : public FunctionInterface<domain_dim, range_dim, range_dim_cols, RangeField>
{
  using BaseType = FunctionInterface<domain_dim, range_dim, range_dim_cols, RangeField>;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::RangeReturnType;
  using typename BaseType::DerivativeRangeReturnType;

  using OrderLambdaType = std::function<int(const Common::Parameter&)>;
  using EvaluateLambdaType = std::function<RangeReturnType(const DomainType&, const Common::Parameter&)>;
  using JacobianLambdaType = std::function<DerivativeRangeReturnType(const DomainType&, const Common::Parameter&)>;
  using DerivativeLambdaType = std::function<DerivativeRangeReturnType(
      const std::array<size_t, d>&, const DomainType&, const Common::Parameter&)>;

  LambdaFunction(OrderLambdaType order_lambda,
                 EvaluateLambdaType evaluate_lambda = default_evaluate_lambda(),
                 const std::string nm = "smooth_lambda_function",
                 const Common::ParameterType& param_type = {},
                 JacobianLambdaType jacobian_lambda = default_jacobian_lambda(),
                 DerivativeLambdaType derivative_lambda = default_derivative_lambda())
    : BaseType(param_type)
    , order_lambda_(order_lambda)
    , evaluate_lambda_(evaluate_lambda)
    , jacobian_lambda_(jacobian_lambda)
    , derivative_lambda_(derivative_lambda)
    , name_(nm)
  {
  }

  LambdaFunction(int ord,
                 EvaluateLambdaType evaluate_lambda = default_evaluate_lambda(),
                 const std::string nm = "smooth_lambda_function",
                 const Common::ParameterType& param_type = {},
                 JacobianLambdaType jacobian_lambda = default_jacobian_lambda(),
                 DerivativeLambdaType derivative_lambda = default_derivative_lambda())
    : BaseType(param_type)
    , order_lambda_([=](const auto& /*param*/) { return ord; })
    , evaluate_lambda_(evaluate_lambda)
    , jacobian_lambda_(jacobian_lambda)
    , derivative_lambda_(derivative_lambda)
    , name_(nm)
  {
  }

  /**
   * \name ´´These methods are required by FunctionInterface.''
   * \{
   */

  int order(const Common::Parameter& param = {}) const override final
  {
    return order_lambda_(this->parse_parameter(param));
  }

  RangeReturnType evaluate(const DomainType& point_in_global_coordinates,
                           const Common::Parameter& param = {}) const override final
  {
    return evaluate_lambda_(point_in_global_coordinates, this->parse_parameter(param));
  }

  DerivativeRangeReturnType jacobian(const DomainType& point_in_global_coordinates,
                                     const Common::Parameter& param = {}) const override final
  {
    return jacobian_lambda_(point_in_global_coordinates, this->parse_parameter(param));
  }

  DerivativeRangeReturnType derivative(const std::array<size_t, d>& alpha,
                                       const DomainType& point_in_global_coordinates,
                                       const Common::Parameter& param = {}) const override final
  {
    return derivative_lambda_(alpha, point_in_global_coordinates, this->parse_parameter(param));
  }

  std::string type() const override final
  {
    return "smooth_lambda_function";
  }

  std::string name() const override final
  {
    return name_;
  }

  /**
   * \}
   * \name ´´These methods may be used to provide defaults on construction.''
   * \{
   */

  static EvaluateLambdaType default_evaluate_lambda()
  {
    return [](const DomainType& /*point_in_global_coordinates*/, const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This LambdaFunction does not provide evaluations, provide an evaluate_lambda on construction!");
      return RangeReturnType();
    };
  }

  static JacobianLambdaType default_jacobian_lambda()
  {
    return [](const DomainType& /*point_in_global_coordinates*/, const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This LambdaFunction does not provide jacobian evaluations, provide a "
                 "jacobian_lambda on construction!");
      return DerivativeRangeReturnType();
    };
  }

  static DerivativeLambdaType default_derivative_lambda()
  {
    return [](const std::array<size_t, d>& /*alpha*/,
              const DomainType& /*point_in_global_coordinates*/,
              const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This LambdaFunction does not provide derivative evaluations, provide a "
                 "derivative_lambda on construction!");
      return DerivativeRangeReturnType();
    };
  }

  /**
   * \}
   */

  const OrderLambdaType order_lambda_;
  const EvaluateLambdaType evaluate_lambda_;
  const JacobianLambdaType jacobian_lambda_;
  const DerivativeLambdaType derivative_lambda_;
  const std::string name_;
}; // class LambdaFunction


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_LAMBDA_SMOOTH_FUNCTION_HH
