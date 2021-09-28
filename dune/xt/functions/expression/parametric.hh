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
//   Tobias Leibner  (2017 - 2020)

#ifndef DUNE_XT_FUNCTIONS_EXPRESSION_PARAMETRIC_HH
#define DUNE_XT_FUNCTIONS_EXPRESSION_PARAMETRIC_HH

#include <limits>

#include <dune/xt/common/parameter.hh>
#include "dune/xt/functions/interfaces/function.hh"

#include "base.hh"

namespace Dune::XT::Functions {


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
  using ThisType = ParametricExpressionFunction;
  using BaseType = FunctionInterface<d, r, 1, R>;

public:
  using BaseType::domain_dim;
  using BaseType::range_dim;
  using typename BaseType::D;
  using typename BaseType::DerivativeRangeReturnType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeReturnType;

private:
  using ActualFunctionType = DynamicMathExpressionBase<D, R, r>;

public:
  static std::string static_id()
  {
    return BaseType::static_id() + ".parametricexpression";
  }

  ParametricExpressionFunction(const std::string& variable,
                               const Common::ParameterType& param_type,
                               const Common::FieldVector<std::string, r>& expressions,
                               const size_t ord = 0,
                               const std::string& nm = static_id())
    : BaseType(param_type)
    , order_(ord)
    , name_(nm)
    , num_parameter_variables_(0)
  {
    DUNE_THROW_IF(variable.empty(), Common::Exceptions::wrong_input_given, "Given variable must not be empty!");

    std::vector<std::string> variables;
    std::vector<std::string> expression;
    for (size_t rr = 0; rr < r; ++rr)
      expression.emplace_back(expressions[rr]);
    for (const auto& key : this->parameter_type().keys()) {
      const size_t value_size = this->parameter_type().get(key);
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
    function_ = std::make_unique<ActualFunctionType>(variables, expression);
  }

  ParametricExpressionFunction(const ThisType& other)
    : BaseType(other)
    , order_(other.order_)
    , name_(other.name_)
    , num_parameter_variables_(other.num_parameter_variables_)
    , function_(new ActualFunctionType(*other.function_))
  {}

  ParametricExpressionFunction(ThisType&&) = default;

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

  std::string name() const override final
  {
    return name_;
  }

  int order(const Common::Parameter& /*param*/ = {}) const override final
  {
    return static_cast<int>(order_);
  }

  RangeReturnType evaluate(const DomainType& point_in_global_coordinates,
                           const Common::Parameter& param = {}) const override final
  {
    RangeReturnType ret(0.);
    Common::Parameter parsed_param;
    if (!this->parameter_type().empty()) {
      parsed_param = this->parse_parameter(param);
      if (parsed_param.type() != this->parameter_type())
        DUNE_THROW(Common::Exceptions::parameter_error,
                   "parameter_type(): " << this->parameter_type() << "\n   "
                                        << "param.type(): " << param.type());
    }
    DynamicVector<D> args(num_parameter_variables_ + domain_dim);
    size_t II = 0;
    for (const auto& key : this->parameter_type().keys()) {
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
#  ifndef DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS
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
                       << error_type << "\n   "
                       << "The variables of this function are:   " << function_->variables() << "\n   "
                       << "The expressions of this function are: " << function_->expressions() << "\n   "
                       << "You evaluated it with            point_in_global_coordinates : "
                       << point_in_global_coordinates << "\n   "
                       << "                                 param : " << param << "\n   "
                       << "The result was:                       " << ret[rr] << "\n\n"
                       << "You can disable this check by defining DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS\n");
    }
#  endif // DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS
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
  size_t num_parameter_variables_;
  std::unique_ptr<ActualFunctionType> function_;
}; // class ParametricExpressionFunction


} // namespace Dune::XT::Functions

#endif // DUNE_XT_FUNCTIONS_EXPRESSION_PARAMETRIC_HH
