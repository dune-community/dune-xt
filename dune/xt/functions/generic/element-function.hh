// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_XT_FUNCTIONS_GENERIC_ELEMENT_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_GENERIC_ELEMENT_FUNCTION_HH

#include <functional>

#include <dune/xt/common/memory.hh>

#include <dune/xt/functions/interfaces/element-functions.hh>

namespace Dune {
namespace XT {
namespace Functions {


/**
 * Set of element functions you can pass lambda expressions or std::functions to that get evaluated.
 */
template <class E, size_t r, size_t rC, class R = double>
class GenericElementFunctionSet : public ElementFunctionSetInterface<E, r, rC, R>
{
  using BaseType = ElementFunctionSetInterface<E, r, rC, R>;

public:
  using BaseType::d;
  using typename BaseType::DerivativeRangeType;
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::RangeType;

  using GenericSizeFunctionType = std::function<size_t(const Common::Parameter&)>;
  using GenericOrderFunctionType = std::function<int(const Common::Parameter&)>;
  using GenericEvaluateFunctionType =
      std::function<void(const DomainType&, std::vector<RangeType>&, const Common::Parameter&)>;
  using GenericJacobianFunctionType =
      std::function<void(const DomainType&, std::vector<DerivativeRangeType>&, const Common::Parameter&)>;
  using GenericDerivativeFunctionType = std::function<void(
      const std::array<size_t, d>&, const DomainType&, std::vector<DerivativeRangeType>&, const Common::Parameter&)>;
  using GenericPostBindFunctionType = std::function<void(const ElementType&)>;

  GenericElementFunctionSet(GenericSizeFunctionType size_func,
                            GenericSizeFunctionType max_size_func,
                            GenericOrderFunctionType order_func,
                            GenericEvaluateFunctionType evaluate_func = default_evaluate_function(),
                            const Common::ParameterType& param_type = {},
                            GenericJacobianFunctionType jacobian_func = default_jacobian_function(),
                            GenericDerivativeFunctionType derivative_func = default_derivative_function(),
                            GenericPostBindFunctionType post_bind_func = [](const auto&) {})
    : BaseType(param_type)
    , size_(size_func)
    , max_size_(max_size_func)
    , order_(order_func)
    , evaluate_(evaluate_func)
    , jacobian_(jacobian_func)
    , derivative_(derivative_func)
    , post_bind_(post_bind_func)
  {}

  GenericElementFunctionSet(const size_t fixed_size,
                            const int ord,
                            GenericEvaluateFunctionType evaluate_lambda = default_evaluate_function(),
                            const Common::ParameterType& param_type = {},
                            GenericJacobianFunctionType jacobian_lambda = default_jacobian_function(),
                            GenericDerivativeFunctionType derivative_lambda = default_derivative_function(),
                            GenericPostBindFunctionType post_bind_func = [](const auto&) {})
    : BaseType(param_type)
    , size_([=](const auto& /*param*/) { return fixed_size; })
    , max_size_([=](const auto& /*param*/) { return fixed_size; })
    , order_([=](const auto& /*param*/) { return ord; })
    , evaluate_(evaluate_lambda)
    , jacobian_(jacobian_lambda)
    , derivative_(derivative_lambda)
    , post_bind_(post_bind_func)
  {}

  /**
   * \name ´´These methods are required by XT::Grid::ElementBoundObject.''
   * \{
   */

protected:
  void post_bind(const ElementType& element) override final
  {
    post_bind_(element);
  }

public:
  /**
   * \}
   * \name ´´These methods are required by ElementFunctionSetInterface.''
   * \{
   */

  size_t size(const Common::Parameter& param = {}) const override final
  {
    return size_(this->parse_parameter(param));
  }

  size_t max_size(const Common::Parameter& param = {}) const override final
  {
    return max_size_(this->parse_parameter(param));
  }

  int order(const Common::Parameter& param = {}) const override final
  {
    return order_(this->parse_parameter(param));
  }

  void evaluate(const DomainType& point_in_reference_element,
                std::vector<RangeType>& result,
                const Common::Parameter& param = {}) const override final
  {
    const auto sz = this->size(param);
    if (result.size() < sz)
      result.resize(sz);
    evaluate_(point_in_reference_element, result, this->parse_parameter(param));
  }

  void jacobians(const DomainType& point_in_reference_element,
                 std::vector<DerivativeRangeType>& result,
                 const Common::Parameter& param = {}) const override final
  {
    const auto sz = this->size(param);
    if (result.size() < sz)
      result.resize(sz);
    jacobian_(point_in_reference_element, result, this->parse_parameter(param));
  }

  void derivatives(const std::array<size_t, d>& alpha,
                   const DomainType& point_in_reference_element,
                   std::vector<DerivativeRangeType>& result,
                   const Common::Parameter& param = {}) const override final
  {
    const auto sz = this->size(param);
    if (result.size() < sz)
      result.resize(sz);
    derivative_(alpha, point_in_reference_element, result, this->parse_parameter(param));
  }

  /**
   * \}
   * \name ´´These methods may be used to provide defaults on construction.''
   * \{
   */

  static GenericEvaluateFunctionType default_evaluate_function()
  {
    return [](const DomainType& /*point_in_reference_element*/,
              std::vector<RangeType>& /*result*/,
              const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(
          NotImplemented,
          "This GenericElementFunctionSet does not provide evaluations, provide an evaluate_lambda on construction!");
    };
  }

  static GenericJacobianFunctionType default_jacobian_function()
  {
    return [](const DomainType& /*point_in_reference_element*/,
              std::vector<DerivativeRangeType>& /*result*/,
              const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This GenericElementFunctionSet does not provide jacobian evaluations, provide a "
                 "jacobian_lambda on construction!");
    };
  }

  static GenericDerivativeFunctionType default_derivative_function()
  {
    return [](const std::array<size_t, d>& /*alpha*/,
              const DomainType& /*point_in_reference_element*/,
              std::vector<DerivativeRangeType>& /*result*/,
              const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This GenericElementFunctionSet does not provide derivative evaluations, provide a "
                 "derivative_lambda on construction!");
    };
  }

  /**
   * \}
   */

  const GenericSizeFunctionType size_;
  const GenericSizeFunctionType max_size_;
  const GenericOrderFunctionType order_;
  const GenericEvaluateFunctionType evaluate_;
  const GenericJacobianFunctionType jacobian_;
  const GenericDerivativeFunctionType derivative_;
  const GenericPostBindFunctionType post_bind_;
}; // class GenericElementFunctionSet


/**
 * Element function you can pass lambda expressions or std::functions to that get evaluated.
 */
template <class E, size_t r, size_t rC, class R = double>
class GenericElementFunction : public ElementFunctionInterface<E, r, rC, R>
{
  using BaseType = ElementFunctionInterface<E, r, rC, R>;

public:
  using BaseType::d;
  using typename BaseType::DerivativeRangeReturnType;
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::RangeReturnType;

  using GenericOrderFunctionType = std::function<int(const Common::Parameter&)>;
  using GenericEvaluateFunctionType = std::function<RangeReturnType(const DomainType&, const Common::Parameter&)>;
  using GenericJacobianFunctionType =
      std::function<DerivativeRangeReturnType(const DomainType&, const Common::Parameter&)>;
  using GenericDerivativeFunctionType = std::function<DerivativeRangeReturnType(
      const std::array<size_t, d>&, const DomainType&, const Common::Parameter&)>;
  using GenericPostBindFunctionType = std::function<void(const ElementType&)>;

  GenericElementFunction(GenericOrderFunctionType order_func,
                         GenericEvaluateFunctionType evaluate_func = default_evaluate_function(),
                         const Common::ParameterType& param_type = {},
                         GenericJacobianFunctionType jacobian_func = default_jacobian_function(),
                         GenericDerivativeFunctionType derivative_func = default_derivative_function(),
                         GenericPostBindFunctionType post_bind_func = [](const auto&) {})
    : BaseType(param_type)
    , order_(order_func)
    , evaluate_(evaluate_func)
    , jacobian_(jacobian_func)
    , derivative_(derivative_func)
    , post_bind_(post_bind_func)
  {}

  GenericElementFunction(const int ord,
                         GenericEvaluateFunctionType evaluate_lambda = default_evaluate_function(),
                         const Common::ParameterType& param_type = {},
                         GenericJacobianFunctionType jacobian_lambda = default_jacobian_function(),
                         GenericDerivativeFunctionType derivative_lambda = default_derivative_function(),
                         GenericPostBindFunctionType post_bind_func = [](const auto&) {})
    : BaseType(param_type)
    , order_([=](const auto& /*param*/) { return ord; })
    , evaluate_(evaluate_lambda)
    , jacobian_(jacobian_lambda)
    , derivative_(derivative_lambda)
    , post_bind_(post_bind_func)
  {}

  /**
   * \name ´´These methods are required by XT::Grid::ElementBoundObject.''
   * \{
   */

protected:
  void post_bind(const ElementType& element) override final
  {
    post_bind_(element);
  }

public:
  /**
   * \}
   * \name ´´These methods are required by ElementFunctionInterface.''
   * \{
   */

  int order(const Common::Parameter& param = {}) const override final
  {
    return order_(this->parse_parameter(param));
  }

  RangeReturnType evaluate(const DomainType& point_in_reference_element,
                           const Common::Parameter& param = {}) const override final
  {
    return evaluate_(point_in_reference_element, this->parse_parameter(param));
  }

  DerivativeRangeReturnType jacobian(const DomainType& point_in_reference_element,
                                     const Common::Parameter& param = {}) const override final
  {
    return jacobian_(point_in_reference_element, this->parse_parameter(param));
  }

  DerivativeRangeReturnType derivative(const std::array<size_t, d>& alpha,
                                       const DomainType& point_in_reference_element,
                                       const Common::Parameter& param = {}) const override final
  {
    return derivative_(alpha, point_in_reference_element, this->parse_parameter(param));
  }

  /**
   * \}
   * \name ´´These methods may be used to provide defaults on construction.''
   * \{
   */

  static GenericEvaluateFunctionType default_evaluate_function()
  {
    return [](const DomainType& /*point_in_reference_element*/, const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(
          NotImplemented,
          "This GenericElementFunction does not provide evaluations, provide an evaluate_lambda on construction!");
      return RangeReturnType();
    };
  }

  static GenericJacobianFunctionType default_jacobian_function()
  {
    return [](const DomainType& /*point_in_reference_element*/, const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This GenericElementFunction does not provide jacobian evaluations, provide a "
                 "jacobian_lambda on construction!");
      return DerivativeRangeReturnType();
    };
  }

  static GenericDerivativeFunctionType default_derivative_function()
  {
    return [](const std::array<size_t, d>& /*alpha*/,
              const DomainType& /*point_in_reference_element*/,
              const Common::Parameter& /*param*/ = {}) {
      DUNE_THROW(NotImplemented,
                 "This GenericElementFunction does not provide derivative evaluations, provide a "
                 "derivative_lambda on construction!");
      return DerivativeRangeReturnType();
    };
  }

  /**
   * \}
   */

  const GenericOrderFunctionType order_;
  const GenericEvaluateFunctionType evaluate_;
  const GenericJacobianFunctionType jacobian_;
  const GenericDerivativeFunctionType derivative_;
  const GenericPostBindFunctionType post_bind_;
}; // class GenericElementFunction


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_LAMBDA_SMOOTH_FUNCTION_HH
