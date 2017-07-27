// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2017)

#ifndef DUNE_XT_FUNCTIONS_LAMBDA_GLOBAL_FLUX_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_LAMBDA_GLOBAL_FLUX_FUNCTION_HH

#include <functional>

#include <dune/xt/common/memory.hh>

#include <dune/xt/functions/interfaces.hh>

namespace Dune {
namespace XT {
namespace Functions {


template <class U, size_t state_derivative_order = 0, class R = typename U::RangeFieldType, size_t r = 1, size_t rC = 1>
class GlobalLambdaFluxFunction
{
  static_assert(AlwaysFalse<U>::value, "Not available for these dimensions or this state_derivative_order!");
};

/**
 * \brief A flux function given by a lambda expression which can be evaluated without an entity.
 *
 *        To model the flux F in a Burgers equation, F(u) := u^p, with a variable exponent p given a localizable
 *        function u, use as
 * \code
typedef some_localizable_function U;
typedef GlobalLambdaFluxFunction<U> FluxType;

FluxType F([](const typename F::DomainType& x,
              const typename F::StateRangeType& u,
              const XT::Common::Parameter& mu) { return std::pow(u[0], mu.get("power").at(0)); },
           XT::Common::ParameterType("power", 1),
           "burgers_flux",
           [](const XT::Common::Parameter& mu) { return mu.get("power").at(0);});
\endcode
 *        The XT::Common::ParameterType provided on construction ensures that the XT::Common::Parameter mu which is
 *        passed on to the lambda is of correct type.
 *        The (optional) fourth argument is a lambda function that takes a parameter and returns the order of the
 *        GlobalLambdaFluxFunction with this parameter. If you do no provide this order lambda, you can't call
 *        the order() method, everything else will work as expected.
 * \note  This function does not implement jacobians.
 */
template <class U, class R, size_t r, size_t rC>
class GlobalLambdaFluxFunction<U, 0, R, r, rC>
    : public GlobalFluxFunctionInterface<typename U::E, typename U::D, U::d, U, 0, R, r, rC>
{
  typedef GlobalFluxFunctionInterface<typename U::E, typename U::D, U::d, U, 0, R, r, rC> BaseType;

public:
  using typename BaseType::DomainType;
  using typename BaseType::StateRangeType;
  using typename BaseType::RangeType;

private:
  typedef std::function<RangeType(const DomainType&, const StateRangeType&, const Common::Parameter&)> LambdaType;
  typedef std::function<size_t(const Common::Parameter&)> OrderLambdaType;

public:
  GlobalLambdaFluxFunction(LambdaType lambda,
                           const Common::ParameterType& param_type,
                           const std::string nm = "globallambdafunction",
                           OrderLambdaType order_lambda =
                               [](const Common::Parameter&) {
                                 DUNE_THROW(
                                     NotImplemented,
                                     "To call the order method, you have to provide an order lambda on construction!");
                                 return 0;
                               })
    : lambda_(lambda)
    , param_type_(param_type)
    , name_(nm)
    , order_lambda_(order_lambda)
  {
  }

  virtual size_t order(const Common::Parameter& mu) const override final
  {
    auto parsed_mu = this->parse_and_check(mu);
    return order_lambda_(parsed_mu);
  }

  using BaseType::evaluate;

  void evaluate(const DomainType& xx,
                const StateRangeType& uu,
                RangeType& ret,
                const Common::Parameter& mu = {}) const override final
  {
    auto parsed_mu = this->parse_and_check(mu);
    ret = lambda_(xx, uu, parsed_mu);
  }

  std::string type() const override final
  {
    return "globallambdafluxfunction";
  }

  std::string name() const override final
  {
    return name_;
  }

  const Common::ParameterType& parameter_type() const override final
  {
    return param_type_;
  }

private:
  const LambdaType lambda_;
  const Common::ParameterType param_type_;
  const std::string name_;
  const OrderLambdaType order_lambda_;
};


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_LAMBDA_GLOBAL_FLUX_FUNCTION_HH
