// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2020)
//   René Fritze     (2018 - 2020)
//   Tim Keil        (2018)
//   Tobias Leibner  (2017, 2019 - 2020)

#ifndef DUNE_XT_FUNCTIONS_GENERIC_GRID_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_GENERIC_GRID_FUNCTION_HH

#include <functional>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/type_traits.hh>
#include <utility>

namespace Dune::XT::Functions {


/**
 * \brief A function given by a lambda expression or std::function which is evaluated locally on each element.
 *
 *        To model the function f(x) = x^p (where x is the local coordinate on the reference element) with a variable
 *        exponent p, use as
 *        \code
 *        GenericGridFunction<...> f([](const typename F::DomainType& point_in_local_coordinates,
 *                                      const XT::Common::Parameter&  param) {
 *                                     typename F::RangeType ret(std::pow(point_in_local_coordinates[0],
 *                                             param.get("power").at(0)));
 *                                     return ret;
 *                                   },
 *                                   integration_order,
 *                                   XT::Common::ParameterType("power", 1),
 *                                   "x_power_p");
 *        \endcode
 *        The XT::Common::ParameterType provided on construction ensures that the XT::Common::Parameter param which is
 *        passed on to the generic function is of correct type.
 *        If you want the function to have Jacobians, you have to provide an additional generic function to the
 *        constructor.
 *        Note that this generic jacobian function should describe the local jacobian on the reference element (for the
 *        x^p function in one dimension, this would be p x^{p-1}). The jacobian(...) method of the
 *        LocalGenericGridFunction will do the transformation for you and return the jacobian on the actual grid
 *        element.
 */
template <class E, size_t r = 1, size_t rC = 1, class R = double>
class GenericGridFunction : public GridFunctionInterface<E, r, rC, R>
{
  using ThisType = GenericGridFunction;
  using BaseType = GridFunctionInterface<E, r, rC, R>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::LocalFunctionType;

private:
  class LocalGenericGridFunction : public ElementFunctionInterface<E, r, rC, R>
  {
    using BaseType = ElementFunctionInterface<E, r, rC, R>;

  public:
    using typename BaseType::DerivativeRangeReturnType;
    using typename BaseType::DerivativeRangeType;
    using typename BaseType::DomainType;
    using typename BaseType::ElementType;
    using typename BaseType::RangeReturnType;
    using typename BaseType::RangeType;

    using BaseType::d;

    using GenericEvaluateFunctionType = std::function<RangeReturnType(const DomainType&, const Common::Parameter&)>;
    using GenericPostBindFunctionType = std::function<void(const ElementType&)>;
    using GenericOrderFunctionType = std::function<int(const Common::Parameter&)>;
    using GenericJacobianFunctionType =
        std::function<DerivativeRangeReturnType(const DomainType&, const XT::Common::Parameter&)>;
    using GenericDerivativeFunctionType = std::function<DerivativeRangeReturnType(
        const std::array<size_t, d>&, const DomainType&, const XT::Common::Parameter&)>;

    LocalGenericGridFunction(const GenericOrderFunctionType& order_func,
                             const GenericPostBindFunctionType& post_bind_func,
                             const GenericEvaluateFunctionType& evaluate_func,
                             const Common::ParameterType& param_type,
                             const GenericJacobianFunctionType& jacobian_func,
                             const GenericDerivativeFunctionType& derivative_func)
      : BaseType(param_type)
      , order_(order_func)
      , post_bind_(post_bind_func)
      , evaluate_(evaluate_func)
      , jacobian_(jacobian_func)
      , derivative_(derivative_func)
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
                             const Common::Parameter& param = {}) const final
    {
      auto parsed_param = this->parse_parameter(param);
      return evaluate_(point_in_local_coordinates, parsed_param);
    }

    DerivativeRangeReturnType jacobian(const DomainType& point_in_local_coordinates,
                                       const Common::Parameter& param = {}) const final
    {
      auto parsed_param = this->parse_parameter(param);
      auto local_jacobian = jacobian_(point_in_local_coordinates, parsed_param);
      const auto J_inv_T = this->element().geometry().jacobianInverseTransposed(point_in_local_coordinates);
      DerivativeRangeReturnType global_jacobian;
      if constexpr (rC == 1) {
        for (size_t rr = 0; rr < r; ++rr)
          J_inv_T.mv(local_jacobian[rr], global_jacobian[rr]);
      } else {
        for (size_t rr = 0; rr < r; ++rr)
          for (size_t ii = 0; ii < rC; ++ii)
            J_inv_T.mv(local_jacobian[rr][ii], global_jacobian[rr][ii]);
      }
      return global_jacobian;
    }

    DerivativeRangeReturnType derivative(const std::array<size_t, d>& alpha,
                                         const DomainType& point_in_local_coordinates,
                                         const Common::Parameter& param = {}) const final
    {
      DUNE_THROW(Dune::NotImplemented,
                 "This function should also transform the derivatives (like the jacobian method), go ahead and "
                 "implement if you want to use this method!");
      auto parsed_param = this->parse_parameter(param);
      return derivative_(alpha, point_in_local_coordinates, parsed_param);
    }

  private:
    const GenericOrderFunctionType order_;
    const GenericPostBindFunctionType post_bind_;
    const GenericEvaluateFunctionType evaluate_;
    const GenericJacobianFunctionType jacobian_;
    const GenericDerivativeFunctionType derivative_;
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
      std::function<DerivativeRangeReturnType(const DomainType&, const XT::Common::Parameter&)>;
  using GenericDerivativeFunctionType = std::function<DerivativeRangeReturnType(
      const std::array<size_t, d>&, const DomainType&, const XT::Common::Parameter&)>;

  GenericGridFunction(const int ord,
                      GenericPostBindFunctionType post_bind_func = default_post_bind_function(),
                      GenericEvaluateFunctionType evaluate_func = default_evaluate_function(),
                      const Common::ParameterType& param_type = Common::ParameterType(),
                      const std::string& nm = "GenericGridFunction",
                      GenericJacobianFunctionType jacobian_func = default_jacobian_function(),
                      GenericDerivativeFunctionType derivative_func = default_derivative_function())
    : BaseType(param_type)
    , order_(default_order_lambda(ord))
    , post_bind_(post_bind_func)
    , evaluate_(evaluate_func)
    , name_(nm)
    , jacobian_(jacobian_func)
    , derivative_(derivative_func)
  {}

  GenericGridFunction(GenericOrderFunctionType order_func,
                      GenericPostBindFunctionType post_bind_func = default_post_bind_function(),
                      GenericEvaluateFunctionType evaluate_func = default_evaluate_function(),
                      const Common::ParameterType& param_type = Common::ParameterType(),
                      const std::string& nm = "GenericGridFunction",
                      GenericJacobianFunctionType jacobian_func = default_jacobian_function(),
                      GenericDerivativeFunctionType derivative_func = default_derivative_function())
    : BaseType(param_type)
    , order_(std::move(order_func))
    , post_bind_(post_bind_func)
    , evaluate_(evaluate_func)
    , name_(nm)
    , jacobian_(jacobian_func)
    , derivative_(derivative_func)
  {}

  GenericGridFunction(const ThisType&) = default;

  GenericGridFunction(ThisType&&) = default;

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

  std::string name() const final
  {
    return name_;
  }

  std::unique_ptr<LocalFunctionType> local_function() const final
  {
    return std::make_unique<LocalGenericGridFunction>(
        order_, post_bind_, evaluate_, this->parameter_type(), jacobian_, derivative_);
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
  const std::string name_;
  const GenericJacobianFunctionType jacobian_;
  const GenericDerivativeFunctionType derivative_;
}; // class GenericGridFunction


} // namespace Dune::XT::Functions

#endif // DUNE_XT_FUNCTIONS_LAMBDA_LOCAL_FUNCTION_HH
