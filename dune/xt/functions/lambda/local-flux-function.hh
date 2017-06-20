// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_FUNCTIONS_LAMBDA_LOCAL_FLUX_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_LAMBDA_LOCAL_FLUX_FUNCTION_HH

#include <functional>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>
#include <dune/xt/functions/interfaces/localizable-flux-function.hh>
#include <dune/xt/functions/type_traits.hh>

namespace Dune {
namespace XT {
namespace Functions {


template <class U, size_t state_derivative_order = 0, class R = typename U::RangeFieldType, size_t r = 1, size_t rC = 1>
class LocalLambdaFluxFunction
{
  static_assert(AlwaysFalse<U>::value, "Not available for these dimensions or this state_derivative_order!");
};


/**
 * \brief A flux function given by a lambda expression which is evaluated locally on each entity.
 *
 *        To model the flux F in a Burgers equation, F(u) := u^p, with a variable exponent p given a localizable
 *        function u, use as
 * \code
typedef some_localizable_function U;
typedef LocalLambdaFluxFunction<U> FluxType;

FluxType F([](const typename F::EntityType& entity,
              const typename F::DomainType& x,
              const typename F::StateRangeType& u,
              const XT::Common::Parameter& mu) { return std::pow(u[0], mu.get("power").at(0)); },
           XT::Common::ParameterType("power", 1),
           "burgers_flux",
           [](const XT::Common::Parameter& mu) { return mu.get("power").at(0);});
\endcode
 *        The XT::Common::ParameterType provided on construction ensures that the XT::Common::Parameter mu which is
 *        passed on to the lambda is of correct type.
 *        The (optional) fourth argument is a lambda function that takes a parameter and returns the order of the
 *        LocalLambdaFluxLocalFunction with this parameter. If you do no provide this order lambda, you can't call
 *        the order() method of the LocalLambdaFluxLocalFunction, everything else will work as expected.
 * \note  The Localfunction does not implement jacobians.
 */
template <class U, class R, size_t r, size_t rC>
class LocalLambdaFluxFunction<U, 0, R, r, rC>
    : public LocalizableFluxFunctionInterface<typename U::E, typename U::D, U::d, U, 0, R, r, rC>
{
  typedef LocalizableFluxFunctionInterface<typename U::E, typename U::D, U::d, U, 0, R, r, rC> BaseType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::E;
  using typename BaseType::D;
  using BaseType::d;
  using typename BaseType::LocalfunctionType;

private:
  class LocalLambdaFluxLocalFunction
      : public LocalFluxFunctionInterface<E, D, d, typename U::LocalfunctionType, 0, R, r, rC>
  {
    typedef LocalFluxFunctionInterface<E, D, d, typename U::LocalfunctionType, 0, R, r, rC> BaseType;

  public:
    using typename BaseType::EntityType;
    using typename BaseType::DomainType;
    using typename BaseType::StateRangeType;
    using typename BaseType::RangeType;
    using typename BaseType::JacobianWrtURangeType;
    using typename BaseType::JacobianWrtXRangeType;

    typedef std::function<RangeType(
        const EntityType&, const DomainType&, const StateRangeType&, const Common::Parameter&)>
        LambdaType;

    typedef std::function<size_t(const Common::Parameter&)> OrderLambdaType;

    LocalLambdaFluxLocalFunction(const EntityType& ent,
                                 const LambdaType& lambda,
                                 const OrderLambdaType& order_lambda,
                                 const Common::ParameterType& param_type)
      : BaseType(ent)
      , lambda_(lambda)
      , order_lambda_(order_lambda)
      , param_type_(param_type)
    {
    }

    virtual size_t order(const XT::Common::Parameter& mu = Common::Parameter()) const override
    {
      auto parsed_mu = parse_and_check(mu);
      return order_lambda_(parsed_mu);
    }

    void evaluate(const DomainType& x,
                  const StateRangeType& u,
                  RangeType& ret,
                  const Common::Parameter& mu = Common::Parameter()) const override final
    {
      auto parsed_mu = parse_and_check(mu);
      ret = lambda_(this->entity(), x, u, parsed_mu);
    } // ... evaluate(...)

    void jacobian_wrt_x(const DomainType& /*x*/,
                        const StateRangeType& /*u*/,
                        JacobianWrtXRangeType& /*ret*/,
                        const Common::Parameter& /*mu*/ = Common::Parameter()) const override final
    {
      DUNE_THROW(NotImplemented, "");
    }

    void jacobian_wrt_u(const DomainType& /*x*/,
                        const StateRangeType& /*u*/,
                        JacobianWrtURangeType& /*ret*/,
                        const Common::Parameter& /*mu*/ = Common::Parameter()) const override final
    {
      DUNE_THROW(NotImplemented, "");
    }

    const Common::ParameterType& parameter_type() const override final
    {
      return param_type_;
    }

  private:
    Common::Parameter parse_and_check(const Common::Parameter& mu) const
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
      return parsed_mu;
    }

    const LambdaType lambda_;
    const OrderLambdaType order_lambda_;
    const Common::ParameterType param_type_;
  }; // class LocalLambdaFluxLocalFunction

public:
  typedef typename LocalLambdaFluxLocalFunction::DomainType DomainType;
  typedef typename LocalLambdaFluxLocalFunction::StateRangeType StateRangeType;
  typedef typename LocalLambdaFluxLocalFunction::RangeType RangeType;
  typedef std::function<RangeType(
      const EntityType&, const DomainType&, const StateRangeType&, const Common::Parameter&)>
      LambdaType; // we do not use the typedef from LocalLambdaFluxLocalFunction here to document the type of the lambda
  typedef std::function<size_t(const Common::Parameter&)> OrderLambdaType;

  LocalLambdaFluxFunction(LambdaType lambda,
                          const Common::ParameterType& param_type = Common::ParameterType(),
                          const std::string nm = "locallambdafluxfunction",
                          OrderLambdaType order_lambda =
                              [](const Common::Parameter&) {
                                DUNE_THROW(NotImplemented, "");
                                return 0;
                              })
    : lambda_(lambda)
    , param_type_(param_type)
    , name_(nm)
    , order_lambda_(order_lambda)
  {
  }

  const Common::ParameterType& parameter_type() const override final
  {
    return param_type_;
  }

  std::string type() const override final
  {
    return "locallambdafluxfunction";
  }

  std::string name() const override final
  {
    return name_;
  }

  std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override final
  {
    return std::make_unique<LocalLambdaFluxLocalFunction>(entity, lambda_, order_lambda_, param_type_);
  }

private:
  const LambdaType lambda_;
  const Common::ParameterType param_type_;
  const std::string name_;
  const OrderLambdaType order_lambda_;
}; // class LocalLambdaFluxFunction


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_LAMBDA_LOCAL_FLUX_FUNCTION_HH
