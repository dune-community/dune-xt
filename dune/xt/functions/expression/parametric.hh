// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_FUNCTIONS_EXPRESSION_PARAMETRIC_HH
#define DUNE_XT_FUNCTIONS_EXPRESSION_PARAMETRIC_HH

#include <string>
#include <limits>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/parameter.hh>
#include <dune/xt/functions/interfaces/global-function.hh>

#include "base.hh"

namespace Dune {
namespace XT {
namespace Functions {


template <class E, class D, size_t d, class R, size_t r, size_t rC = 1>
class ParametricExpressionFunction
{
  static_assert(AlwaysFalse<E>::value, "Not available for these dimension!");
};


template <class E, class D, size_t d, class R, size_t r>
class ParametricExpressionFunction<E, D, d, R, r, 1> : public GlobalFunctionInterface<E, D, d, R, r, 1>
{
  typedef GlobalFunctionInterface<E, D, d, R, r, 1> BaseType;
  typedef DynamicMathExpressionBase<D, R, r> ActualFunctionType;

public:
  using typename BaseType::DomainType;
  using BaseType::dimDomain;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;
  using BaseType::dimRange;

  static std::string static_id()
  {
    return BaseType::static_id() + ".parametricexpression";
  }

  ParametricExpressionFunction(const std::string& variable,
                               const Common::ParameterType& param_type,
                               const std::vector<std::string>& expressions,
                               const size_t ord = 0,
                               const std::string nm = static_id())
    : order_(ord)
    , name_(nm)
    , param_type_(param_type)
    , num_parameter_variables_(0)
  {
    if (variable.empty())
      DUNE_THROW(Common::Exceptions::wrong_input_given, "Given variable must not be empty!");
    if (expressions.size() != dimRange)
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "dimRange: " << size_t(dimRange) << "\n   "
                              << "expressions.size(): "
                              << expressions.size());
    std::vector<std::string> variables;
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
    for (size_t ii = 0; ii < dimDomain; ++ii)
      variables.push_back(variable + "[" + Common::to_string(ii) + "]");
    function_ = std::make_shared<ActualFunctionType>(variables, expressions);
  }

  std::string type() const override final
  {
    return BaseType::static_id() + ".parametricexpression";
  }

  std::string name() const override final
  {
    return name_;
  }

  virtual size_t order(const Common::Parameter& /*mu*/ = {}) const override final
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

  void evaluate(const DomainType& xx, RangeType& ret, const Common::Parameter& mu = {}) const override final
  {
    Common::Parameter parsed_mu;
    if (!param_type_.empty()) {
      parsed_mu = this->parse_parameter(mu);
      if (parsed_mu.type() != param_type_)
        DUNE_THROW(Common::Exceptions::parameter_error,
                   "parameter_type(): " << param_type_ << "\n   "
                                        << "mu.type(): "
                                        << mu.type());
    }
    DynamicVector<D> args(num_parameter_variables_ + dimDomain);
    size_t II = 0;
    for (const auto& key : param_type_.keys()) {
      for (const auto& value : parsed_mu.get(key)) {
        args[II] = value;
        ++II;
      }
    }
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      args[II] = xx[ii];
      ++II;
    }
    function_->evaluate(args, ret);

#ifndef NDEBUG
#ifndef DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS
    bool failure = false;
    std::string error_type;
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (Dune::XT::Common::isnan(ret[rr])) {
        failure = true;
        error_type = "NaN";
      } else if (Dune::XT::Common::isinf(ret[rr])) {
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
                       << "You evaluated it with            xx : "
                       << xx
                       << "\n   "
                       << "                                 mu : "
                       << mu
                       << "\n   "
                       << "The result was:                       "
                       << ret[rr]
                       << "\n\n"
                       << "You can disable this check by defining DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS\n");
    }
#endif // DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS
#endif // NDEBUG
  } // ... evaluate(...)

  void jacobian(const DomainType& /*xx*/,
                JacobianRangeType& /*ret*/,
                const Common::Parameter& /*mu*/ = {}) const override final
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
