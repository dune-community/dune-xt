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

#ifndef DUNE_XT_FUNCTIONS_EXPRESSION_PARAMETRIC_HH
#define DUNE_XT_FUNCTIONS_EXPRESSION_PARAMETRIC_HH

#include <limits>

#include <dune/xt/common/parameter.hh>
#include "dune/xt/functions/interfaces/function.hh"

#include "base.hh"

namespace Dune {
namespace XT {
namespace Functions {


template <size_t d, size_t r = 1, size_t rC = 1, class R = double>
class ParametricExpressionFunction : public FunctionInterface<d, r, rC, R>
{
public:
  ParametricExpressionFunction()
  {
    static_assert(AlwaysFalse<R>::value, "Not available for these dimension!");
  }
};


template <size_t d, size_t r, class R>
class ParametricExpressionFunction<d, r, 1, R> : public FunctionInterface<d, r, 1, R>
{
  using BaseType = FunctionInterface<d, r, 1, R>;
  using typename BaseType::D;
  using ActualFunctionType = DynamicMathExpressionBase<D, R, r>;

public:
  using typename BaseType::DomainType;
  using BaseType::domain_dim;
  using typename BaseType::RangeReturnType;
  using typename BaseType::DerivativeRangeReturnType;
  using BaseType::range_dim;

  static std::string static_id()
  {
    return BaseType::static_id() + ".parametricexpression";
  }

  ParametricExpressionFunction(const std::string& variable,
                               const Common::ParameterType& param_type,
                               const Common::FieldVector<std::string, r>& expressions,
                               const size_t ord = 0,
                               const std::string nm = static_id())
    : order_(ord)
    , name_(nm)
    , param_type_(param_type)
    , num_parameter_variables_(0)
  {
    if (variable.empty())
      DUNE_THROW(Common::Exceptions::wrong_input_given, "Given variable must not be empty!");

    std::vector<std::string> variables;
    std::vector<std::string> expression;
    for (size_t rr = 0; rr < r; ++rr)
      expression.emplace_back(expressions[rr]);
    for (const auto& key : param_type_.keys()) {
      const size_t value_size = param_type_.get(key);
      if (value_size == 1) {
        variables.push_back(key);
        ++num_parameter_variables_;
      } else {
        for (size_t ii = 0; ii < value_size; ++ii) {
          variables.push_back(key + "[" + Common::to_string(ii) + "]");
          ++num_parameter_variables_;
        }
      }
    }
    for (size_t ii = 0; ii < domain_dim; ++ii)
      variables.push_back(variable + "[" + Common::to_string(ii) + "]");
    function_ = std::make_shared<ActualFunctionType>(variables, expression);
  }

  std::string type() const override final
  {
    return BaseType::static_id() + ".parametricexpression";
  }

  std::string name() const override final
  {
    return name_;
  }

  virtual int order(const Common::Parameter& /*param*/ = {}) const override final
  {
    return order_;
  }

  bool is_parametric() const override final
  {
    return !param_type_.empty();
  }

  const Common::ParameterType& parameter_type() const override final
  {
    return param_type_;
  }

  RangeReturnType evaluate(const DomainType& point_in_global_coordinates,
                           const Common::Parameter& param = {}) const override final
  {
    RangeReturnType ret(0.);
    Common::Parameter parsed_param;
    if (!param_type_.empty()) {
      parsed_param = this->parse_parameter(param);
      if (parsed_param.type() != param_type_)
        DUNE_THROW(Common::Exceptions::parameter_error,
                   "parameter_type(): " << param_type_ << "\n   "
                                        << "param.type(): "
                                        << param.type());
    }
    DynamicVector<D> args(num_parameter_variables_ + domain_dim);
    size_t II = 0;
    for (const auto& key : param_type_.keys()) {
      for (const auto& value : parsed_param.get(key)) {
        args[II] = value;
        ++II;
      }
    }
    for (size_t ii = 0; ii < domain_dim; ++ii) {
      args[II] = point_in_global_coordinates[ii];
      ++II;
    }
    function_->evaluate(args, ret);

#ifndef NDEBUG
#ifndef DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS
    bool failure = false;
    std::string error_type;
    for (size_t rr = 0; rr < range_dim; ++rr) {
      if (Common::isnan(ret[rr])) {
        failure = true;
        error_type = "NaN";
      } else if (Common::isinf(ret[rr])) {
        failure = true;
        error_type = "inf";
      } else if (std::abs(ret[rr]) > (0.9 * std::numeric_limits<R>::max())) {
        failure = true;
        error_type = "an unlikely value";
      }
      if (failure)
        DUNE_THROW(Common::Exceptions::internal_error,
                   "evaluating this function yielded:     "
                       << error_type
                       << "\n   "
                       << "The variables of this function are:   "
                       << function_->variables()
                       << "\n   "
                       << "The expressions of this function are: "
                       << function_->expressions()
                       << "\n   "
                       << "You evaluated it with            point_in_global_coordinates : "
                       << point_in_global_coordinates
                       << "\n   "
                       << "                                 param : "
                       << param
                       << "\n   "
                       << "The result was:                       "
                       << ret[rr]
                       << "\n\n"
                       << "You can disable this check by defining DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS\n");
    }
#endif // DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS
#endif // NDEBUG
    return ret;
  } // ... evaluate(...)

  DerivativeRangeReturnType jacobian(const DomainType& /*point_in_global_coordinates*/,
                                     const Common::Parameter& /*param*/ = {}) const override final
  {
    DUNE_THROW(NotImplemented, "Not yet, at least...");
  }

private:
  size_t order_;
  std::string name_;
  Common::ParameterType param_type_;
  size_t num_parameter_variables_;
  std::shared_ptr<const ActualFunctionType> function_;
}; // class ParametricExpressionFunction


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_EXPRESSION_PARAMETRIC_HH
