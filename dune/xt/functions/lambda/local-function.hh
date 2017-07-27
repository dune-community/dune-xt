// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

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
 * \brief A function given by a lambda expression which is evaluated locally on each entity.
 *
 *        To model the function f(x) = x^p with a variable exponent p, use as
 * \code
LocalLambdaFunction<...> f(
                           [](const typename F::EntityType& entity,
                              const typename F::DomainType& xx,
                              const XT::Common::Parameter& mu) {
                             typename F::RangeType ret(std::pow(xx[0], mu.get("power").at(0)));
                             return ret;
                           },
                           integration_order,
                           XT::Common::ParameterType("power", 1),
                           "x_power_p");
\endcode
 *        The XT::Common::ParameterType provided on construction ensures that the XT::Common::Parameter mu which is
 *        passed on to the lambda is of correct type.
 * \note  The Localfunction does not implement jacobian.
 */
template <class E, class D, size_t d, class R, size_t r, size_t rC = 1>
class LocalLambdaFunction : public LocalizableFunctionInterface<E, D, d, R, r, rC>
{
  typedef LocalizableFunctionInterface<E, D, d, R, r, rC> BaseType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::LocalfunctionType;

private:
  class LocalLambdaLocalFunction : public LocalfunctionInterface<E, D, d, R, r, rC>
  {
    typedef LocalfunctionInterface<E, D, d, R, r, rC> BaseType;

  public:
    using typename BaseType::EntityType;
    using typename BaseType::DomainType;
    using typename BaseType::RangeType;
    using typename BaseType::JacobianRangeType;

    typedef std::function<RangeType(const EntityType&, const DomainType&, const Common::Parameter&)> LambdaType;
    typedef std::function<size_t(const Common::Parameter&)> OrderLambdaType;

    LocalLambdaLocalFunction(const EntityType& ent,
                             const LambdaType& lambda,
                             const OrderLambdaType& order_lambda,
                             const Common::ParameterType& param_type)
      : BaseType(ent)
      , lambda_(lambda)
      , order_lambda_(order_lambda)
      , param_type_(param_type)
    {
    }

    virtual size_t order(const XT::Common::Parameter& mu = {}) const override final
    {
      auto parsed_mu = this->parse_and_check(mu);
      return order_lambda_(parsed_mu);
    }

    void evaluate(const DomainType& x, RangeType& ret, const Common::Parameter& mu = {}) const override final
    {
      auto parsed_mu = this->parse_and_check(mu);
      ret = lambda_(this->entity(), x, parsed_mu);
    } // ... evaluate(...)

    void jacobian(const DomainType& /*x*/,
                  JacobianRangeType& /*ret*/,
                  const Common::Parameter& /*mu*/ = {}) const override final
    {
      DUNE_THROW(NotImplemented, "");
    }

    const Common::ParameterType& parameter_type() const override final
    {
      return param_type_;
    }

  private:
    const LambdaType& lambda_;
    const OrderLambdaType& order_lambda_;
    const Common::ParameterType param_type_;
  }; // class LocalLambdaLocalFunction

public:
  typedef typename LocalLambdaLocalFunction::DomainType DomainType;
  typedef typename LocalLambdaLocalFunction::RangeType RangeType;
  // we do not use the typedef from LocalLambdaLocalFunction here to document the type of the lambda
  typedef std::function<RangeType(const EntityType&, const DomainType&, const Common::Parameter&)> LambdaType;
  typedef std::function<size_t(const Common::Parameter&)> OrderLambdaType;

  LocalLambdaFunction(LambdaType lambda,
                      const size_t ord,
                      const Common::ParameterType& param_type = Common::ParameterType(),
                      const std::string nm = "locallambdafunction")
    : lambda_(lambda)
    , order_lambda_([=](const Common::Parameter&) { return ord; })
    , param_type_(param_type)
    , name_(nm)
  {
  }

  LocalLambdaFunction(LambdaType lambda,
                      OrderLambdaType order_lambda,
                      const Common::ParameterType& param_type = Common::ParameterType(),
                      const std::string nm = "locallambdafunction")
    : lambda_(lambda)
    , order_lambda_(order_lambda)
    , param_type_(param_type)
    , name_(nm)
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

  std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override final
  {
    return std::make_unique<LocalLambdaLocalFunction>(entity, lambda_, order_lambda_, param_type_);
  }

private:
  const LambdaType lambda_;
  const OrderLambdaType order_lambda_;
  const Common::ParameterType param_type_;
  const std::string name_;
}; // class LocalLambdaFunction


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_LAMBDA_LOCAL_FUNCTION_HH
