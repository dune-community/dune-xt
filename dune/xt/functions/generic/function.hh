// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2020)
//   René Fritze     (2013 - 2016, 2018 - 2020)
//   Sven Kaulmann   (2013)
//   Tim Keil        (2018)
//   Tobias Leibner  (2014, 2017, 2019 - 2020)

#ifndef DUNE_XT_FUNCTIONS_GENERIC_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_GENERIC_FUNCTION_HH

#include <functional>

#include <dune/xt/common/memory.hh>

#include <dune/xt/functions/interfaces/function.hh>
#include <utility>

namespace Dune::XT::Functions {


/**
 * Smooth function you can pass lambda expressions or std::functions to that gets evaluated.
 *
 * \example LambdaType lambda(1, [](const auto& x, const auto& param = {}) { return x;});
 */
template <size_t domain_dim, size_t range_dim = 1, size_t range_dim_cols = 1, class RangeField = double>
class GenericFunction : public FunctionInterface<domain_dim, range_dim, range_dim_cols, RangeField>
{
  using ThisType = GenericFunction;
  using BaseType = FunctionInterface<domain_dim, range_dim, range_dim_cols, RangeField>;

public:
  using BaseType::d;
  using typename BaseType::DerivativeRangeReturnType;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicRangeType;
  using typename BaseType::RangeReturnType;

  using GenericOrderFunctionType = std::function<int(const Common::Parameter&)>;
  using GenericEvaluateFunctionType = std::function<RangeReturnType(const DomainType&, const Common::Parameter&)>;
  using GenericDynamicEvaluateFunctionType =
      std::function<void(const DomainType&, DynamicRangeType&, const Common::Parameter&)>;
  using GenericJacobianFunctionType =
      std::function<DerivativeRangeReturnType(const DomainType&, const Common::Parameter&)>;
  using GenericDerivativeFunctionType = std::function<DerivativeRangeReturnType(
      const std::array<size_t, d>&, const DomainType&, const Common::Parameter&)>;

  GenericFunction(GenericOrderFunctionType order_func,
                  GenericEvaluateFunctionType evaluate_func = default_evaluate_function(),
                  const std::string& nm = "smooth_lambda_function",
                  const Common::ParameterType& param_type = {},
                  GenericJacobianFunctionType jacobian_func = default_jacobian_function(),
                  GenericDerivativeFunctionType derivative_func = default_derivative_function())
    : BaseType(param_type)
    , order_(std::move(order_func))
    , evaluate_(evaluate_func)
    , dynamic_evaluate_(static_to_dynamic_evaluate(evaluate_))
    , jacobian_(jacobian_func)
    , derivative_(derivative_func)
    , name_(nm)
  {}

  GenericFunction(int ord,
                  GenericEvaluateFunctionType evaluate_func = default_evaluate_function(),
                  const std::string& nm = "smooth_lambda_function",
                  const Common::ParameterType& param_type = {},
                  GenericJacobianFunctionType jacobian_func = default_jacobian_function(),
                  GenericDerivativeFunctionType derivative_func = default_derivative_function())
    : BaseType(param_type)
    , order_([=](const auto& /*param*/) { return ord; })
    , evaluate_(evaluate_func)
    , dynamic_evaluate_(static_to_dynamic_evaluate(evaluate_))
    , jacobian_(jacobian_func)
    , derivative_(derivative_func)
    , name_(nm)
  {}

  GenericFunction(GenericOrderFunctionType order_func,
                  GenericDynamicEvaluateFunctionType dynamic_evaluate_func,
                  const std::string& nm = "smooth_lambda_function",
                  const Common::ParameterType& param_type = {},
                  GenericJacobianFunctionType jacobian_func = default_jacobian_function(),
                  GenericDerivativeFunctionType derivative_func = default_derivative_function())
    : BaseType(param_type)
    , order_(std::move(order_func))
    , evaluate_(dynamic_to_static_evaluate(dynamic_evaluate_func))
    , dynamic_evaluate_(dynamic_evaluate_func)
    , jacobian_(jacobian_func)
    , derivative_(derivative_func)
    , name_(nm)
  {}

  GenericFunction(int ord,
                  GenericDynamicEvaluateFunctionType dynamic_evaluate_func,
                  const std::string& nm = "smooth_lambda_function",
                  const Common::ParameterType& param_type = {},
                  GenericJacobianFunctionType jacobian_func = default_jacobian_function(),
                  GenericDerivativeFunctionType derivative_func = default_derivative_function())
    : BaseType(param_type)
    , order_([=](const auto& /*param*/) { return ord; })
    , evaluate_(dynamic_to_static_evaluate(dynamic_evaluate_func))
    , dynamic_evaluate_(dynamic_evaluate_func)
    , jacobian_(jacobian_func)
    , derivative_(derivative_func)
    , name_(nm)
  {}

  GenericFunction(const ThisType&) = default;

  GenericFunction(ThisType&&) = default;

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

  int order(const Common::Parameter& param = {}) const override final
  {
    return order_(this->parse_parameter(param));
  }

  RangeReturnType evaluate(const DomainType& point_in_global_coordinates,
                           const Common::Parameter& param = {}) const override final
  {
    return evaluate_(point_in_global_coordinates, this->parse_parameter(param));
  }

  void evaluate(const DomainType& point_in_global_coordinates,
                DynamicRangeType& ret,
                const Common::Parameter& param = {}) const override final
  {
    dynamic_evaluate_(point_in_global_coordinates, ret, this->parse_parameter(param));
  }

  DerivativeRangeReturnType jacobian(const DomainType& point_in_global_coordinates,
                                     const Common::Parameter& param = {}) const override final
  {
    return jacobian_(point_in_global_coordinates, this->parse_parameter(param));
  }

  DerivativeRangeReturnType derivative(const std::array<size_t, d>& alpha,
                                       const DomainType& point_in_global_coordinates,
                                       const Common::Parameter& param = {}) const override final
  {
    return derivative_(alpha, point_in_global_coordinates, this->parse_parameter(param));
  }

  std::string name() const override final
  {
    return name_;
  }

  /**
   * \name ´´These methods may be used to provide defaults on construction.''
   * \{
   */

  static GenericEvaluateFunctionType default_evaluate_function()
  {
    return [](const DomainType& /*point_in_global_coordinates*/, const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This GenericFunction does not provide evaluations, provide an evaluate_lambda on construction!");
      return RangeReturnType{};
    };
  }

  static GenericDynamicEvaluateFunctionType default_dynamic_evaluate_function()
  {
    return [](const DomainType& /*point_in_global_coordinates*/,
              DynamicRangeType& /*ret*/,
              const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This GenericFunction does not provide evaluations, provide an evaluate_lambda on construction!");
    };
  }

  static GenericEvaluateFunctionType dynamic_to_static_evaluate(GenericDynamicEvaluateFunctionType dynamic_eval)
  {
    return [=](const DomainType& point_in_global_coordinates, const Common::Parameter& param = {}) {
      DynamicRangeType ret(range_dim);
      dynamic_eval(point_in_global_coordinates, ret, param);
      return XT::Common::convert_to<RangeReturnType>(ret);
    };
  }

  static GenericDynamicEvaluateFunctionType static_to_dynamic_evaluate(GenericEvaluateFunctionType static_eval)
  {
    return
        [=](const DomainType& point_in_global_coordinates, DynamicRangeType& ret, const Common::Parameter& param = {}) {
          ret = XT::Common::convert_to<DynamicRangeType>(static_eval(point_in_global_coordinates, param));
        };
  }

  static GenericJacobianFunctionType default_jacobian_function()
  {
    return [](const DomainType& /*point_in_global_coordinates*/, const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This GenericFunction does not provide jacobian evaluations, provide a "
                 "jacobian_lambda on construction!");
      return DerivativeRangeReturnType();
    };
  }

  static GenericDerivativeFunctionType default_derivative_function()
  {
    return [](const std::array<size_t, d>& /*alpha*/,
              const DomainType& /*point_in_global_coordinates*/,
              const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This GenericFunction does not provide derivative evaluations, provide a "
                 "derivative_lambda on construction!");
      return DerivativeRangeReturnType();
    };
  }

  /**
   * \}
   */

  const GenericOrderFunctionType order_;
  const GenericEvaluateFunctionType evaluate_;
  const GenericDynamicEvaluateFunctionType dynamic_evaluate_;
  const GenericJacobianFunctionType jacobian_;
  const GenericDerivativeFunctionType derivative_;
  const std::string name_;
}; // class GenericFunction


} // namespace Dune::XT::Functions

#endif // DUNE_XT_FUNCTIONS_GENERIC_FUNCTION_HH
