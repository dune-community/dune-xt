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

#ifndef DUNE_XT_FUNCTIONS_LAMBDA_LOCAL_FLUX_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_LAMBDA_LOCAL_FLUX_FUNCTION_HH
#if 0
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
 *        The optional fourth argument is a lambda function that takes a parameter and returns the order of the
 *        LocalLambdaFluxLocalFunction with this parameter. If you do no provide this order lambda, you can't call
 *        the order() method of the LocalLambdaFluxLocalFunction, everything else will work as expected.
 *        The optional fifth and sixth argument are lambda functions that represent the partial derivative with
 *        respect to x and u, respectively. You can either provide a lambda that returns Partial{X,U}RangeType or
 *        a FieldVector<ColPartial{X,U}LambdaType, dimRangeCols> where each entry is a lambda returning the partial
 *        derivative for the respective column. If you do no provide these derivative lambdas, you can't call
 *        the partial_x(...) and partial_u(...) methods of the LocalLambdaFluxLocalFunction.
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
  using BaseType::dimRangeCols;

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
    using typename BaseType::PartialURangeType;
    using typename BaseType::PartialXRangeType;
    using typename BaseType::ColPartialXRangeType;
    using typename BaseType::ColPartialURangeType;

    typedef std::function<RangeType(
        const EntityType&, const DomainType&, const StateRangeType&, const Common::Parameter&)>
        LambdaType;
    typedef std::function<size_t(const Common::Parameter&)> OrderLambdaType;
    typedef std::function<PartialXRangeType(
        const EntityType&, const DomainType&, const StateRangeType&, const Common::Parameter&)>
        PartialXLambdaType;
    typedef std::function<PartialURangeType(
        const EntityType&, const DomainType&, const StateRangeType&, const Common::Parameter&)>
        PartialULambdaType;
    typedef std::function<ColPartialXRangeType(
        const EntityType&, const DomainType&, const StateRangeType&, const Common::Parameter&)>
        ColPartialXLambdaType;
    typedef std::function<ColPartialURangeType(
        const EntityType&, const DomainType&, const StateRangeType&, const Common::Parameter&)>
        ColPartialULambdaType;

    LocalLambdaFluxLocalFunction(const EntityType& ent,
                                 const LambdaType& lambda,
                                 const OrderLambdaType& order_lambda,
                                 const Common::ParameterType& param_type,
                                 const PartialXLambdaType& partial_x_lambda,
                                 const PartialULambdaType& partial_u_lambda,
                                 const FieldVector<ColPartialXLambdaType, dimRangeCols>& partial_x_col_lambdas,
                                 const FieldVector<ColPartialULambdaType, dimRangeCols>& partial_u_col_lambdas)

      : BaseType(ent)
      , lambda_(lambda)
      , order_lambda_(order_lambda)
      , param_type_(param_type)
      , partial_x_lambda_(partial_x_lambda)
      , partial_u_lambda_(partial_u_lambda)
      , partial_x_col_lambdas_(partial_x_col_lambdas)
      , partial_u_col_lambdas_(partial_u_col_lambdas)
    {
    }

    virtual size_t order(const XT::Common::Parameter& mu = {}) const override
    {
      auto parsed_mu = this->parse_parameter(mu);
      return order_lambda_(parsed_mu);
    }

    void evaluate(const DomainType& x,
                  const StateRangeType& u,
                  RangeType& ret,
                  const Common::Parameter& mu = {}) const override final
    {
      auto parsed_mu = this->parse_parameter(mu);
      ret = lambda_(this->entity(), x, u, parsed_mu);
    } // ... evaluate(...)

    void partial_x(const DomainType& /*x*/,
                   const StateRangeType& /*u*/,
                   PartialXRangeType& /*ret*/,
                   const Common::Parameter& /*mu*/ = {}) const override final
    {
      DUNE_THROW(NotImplemented, "");
    }

    void partial_u(const DomainType& /*x*/,
                   const StateRangeType& /*u*/,
                   PartialURangeType& /*ret*/,
                   const Common::Parameter& /*mu*/ = {}) const override final
    {
      DUNE_THROW(NotImplemented, "");
    }

    const Common::ParameterType& parameter_type() const override final
    {
      return param_type_;
    }

  private:
    const LambdaType lambda_;
    const OrderLambdaType order_lambda_;
    const Common::ParameterType param_type_;
    const PartialXLambdaType& partial_x_lambda_;
    const PartialULambdaType& partial_u_lambda_;
    const FieldVector<ColPartialXLambdaType, dimRangeCols>& partial_x_col_lambdas_;
    const FieldVector<ColPartialULambdaType, dimRangeCols>& partial_u_col_lambdas_;
  }; // class LocalLambdaFluxLocalFunction

public:
  typedef typename LocalLambdaFluxLocalFunction::DomainType DomainType;
  typedef typename LocalLambdaFluxLocalFunction::StateRangeType StateRangeType;
  typedef typename LocalLambdaFluxLocalFunction::RangeType RangeType;
  typedef typename LocalLambdaFluxLocalFunction::PartialXRangeType PartialXRangeType;
  typedef typename LocalLambdaFluxLocalFunction::PartialURangeType PartialURangeType;
  typedef typename LocalLambdaFluxLocalFunction::ColPartialXRangeType ColPartialXRangeType;
  typedef typename LocalLambdaFluxLocalFunction::ColPartialURangeType ColPartialURangeType;
  // we do not use the typedef from LocalLambdaFluxLocalFunction here to document the type of the lambda
  typedef std::function<RangeType(
      const EntityType&, const DomainType&, const StateRangeType&, const Common::Parameter&)>
      LambdaType;
  typedef std::function<size_t(const Common::Parameter&)> OrderLambdaType;
  typedef std::function<PartialXRangeType(
      const EntityType&, const DomainType&, const StateRangeType&, const Common::Parameter&)>
      PartialXLambdaType;
  typedef std::function<PartialURangeType(
      const EntityType&, const DomainType&, const StateRangeType&, const Common::Parameter&)>
      PartialULambdaType;
  typedef std::function<ColPartialXRangeType(
      const EntityType&, const DomainType&, const StateRangeType&, const Common::Parameter&)>
      ColPartialXLambdaType;
  typedef std::function<ColPartialURangeType(
      const EntityType&, const DomainType&, const StateRangeType&, const Common::Parameter&)>
      ColPartialULambdaType;

  LocalLambdaFluxFunction(
      LambdaType lambda,
      const Common::ParameterType& param_type = Common::ParameterType(),
      const std::string nm = "locallambdafluxfunction",
      OrderLambdaType order_lambda =
          [](const Common::Parameter&) {
            DUNE_THROW(NotImplemented,
                       "To call the order method, you have to provide an order lambda on construction!");
            return 0;
          },
      PartialXLambdaType partial_x_lambda =
          [](const EntityType&, const DomainType&, const StateRangeType&, const Common::Parameter&) {
            DUNE_THROW(NotImplemented,
                       "You need to provide a lambda for the partial x derivative if you want to use it!");
            return PartialXRangeType();
          },
      PartialULambdaType partial_u_lambda =
          [](const EntityType&, const DomainType&, const StateRangeType&, const Common::Parameter&) {
            DUNE_THROW(NotImplemented,
                       "You need to provide a lambda for the partial u derivative if you want to use it!");
            return PartialURangeType();
          })
    : lambda_(lambda)
    , param_type_(param_type)
    , name_(nm)
    , order_lambda_(order_lambda)
    , partial_x_lambda_(partial_x_lambda)
    , partial_u_lambda_(partial_u_lambda)
  {
    create_col_lambdas(partial_x_lambda_, partial_u_lambda_);
  }

  LocalLambdaFluxFunction(LambdaType lambda,
                          const Common::ParameterType& param_type,
                          const std::string nm,
                          OrderLambdaType order_lambda,
                          FieldVector<ColPartialXLambdaType, dimRangeCols> partial_x_col_lambdas,
                          FieldVector<ColPartialULambdaType, dimRangeCols> partial_u_col_lambdas)
    : lambda_(lambda)
    , param_type_(param_type)
    , name_(nm)
    , order_lambda_(order_lambda)
    , partial_x_col_lambdas_(partial_x_col_lambdas)
    , partial_u_col_lambdas_(partial_u_col_lambdas)
  {
    create_lambdas(partial_x_col_lambdas_, partial_u_col_lambdas_);
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
    return std::make_unique<LocalLambdaFluxLocalFunction>(entity,
                                                          lambda_,
                                                          order_lambda_,
                                                          param_type_,
                                                          partial_x_lambda_,
                                                          partial_u_lambda_,
                                                          partial_x_col_lambdas_,
                                                          partial_u_col_lambdas_);
  }

private:
  void create_lambdas(const FieldVector<ColPartialXLambdaType, dimRangeCols>& partial_x_col_lambdas,
                      const FieldVector<ColPartialXLambdaType, dimRangeCols>& partial_u_col_lambdas)
  {
    partial_x_lambda_ =
        [&](const EntityType& ent, const DomainType& x, const StateRangeType& u, const Common::Parameter& mu) {
          PartialXRangeType ret;
          for (size_t ii = 0; ii < dimRangeCols; ++ii)
            helper<dimRangeCols, PartialXRangeType, ColPartialXRangeType>::set_col_jacobian(
                ii, ret, partial_x_col_lambdas[ii](ent, x, u, mu));
        };
    partial_u_lambda_ =
        [&](const EntityType& ent, const DomainType& x, const StateRangeType& u, const Common::Parameter& mu) {
          PartialURangeType ret;
          for (size_t ii = 0; ii < dimRangeCols; ++ii)
            helper<dimRangeCols, PartialURangeType, ColPartialURangeType>::set_col_jacobian(
                ii, ret, partial_u_col_lambdas[ii](ent, x, u, mu));
        };
  }

  void create_col_lambdas(const PartialXLambdaType& partial_x_lambda, const PartialULambdaType& partial_u_lambda)
  {
    for (size_t ii = 0; ii < dimRangeCols; ++ii) {
      partial_x_col_lambdas_[ii] =
          [&](const EntityType& ent, const DomainType& x, const StateRangeType& u, const Common::Parameter& mu) {
            return helper<dimRangeCols, PartialXRangeType, ColPartialXRangeType>::get_col_jacobian(
                ii, partial_x_lambda(ent, x, u, mu));
          };
      partial_u_col_lambdas_[ii] =
          [&](const EntityType& ent, const DomainType& x, const StateRangeType& u, const Common::Parameter& mu) {
            return helper<dimRangeCols, PartialURangeType, ColPartialURangeType>::get_col_jacobian(
                ii, partial_u_lambda(ent, x, u, mu));
          };
    }
  }

  template <size_t rangeDimCols, class JacobianRangeType, class ColJacobianRangeType>
  struct helper
  {
    static void set_col_jacobian(size_t col, JacobianRangeType& jacobian, const ColJacobianRangeType& jacobian_col)
    {
      jacobian[col] = jacobian_col;
    }

    static const ColJacobianRangeType& get_col_jacobian(size_t col, const JacobianRangeType& jacobian)
    {
      return jacobian[col];
    }
  };

  template <class JacobianRangeType, class ColJacobianRangeType>
  struct helper<1, JacobianRangeType, ColJacobianRangeType>
  {
    static void set_col_jacobian(size_t col, JacobianRangeType& jacobian, const ColJacobianRangeType& jacobian_col)
    {
      assert(col == 0);
      jacobian = jacobian_col;
    }

    static const ColJacobianRangeType& get_col_jacobian(size_t col, const JacobianRangeType& jacobian)
    {
      assert(col == 0);
      return jacobian;
    }
  };

  const LambdaType lambda_;
  const Common::ParameterType param_type_;
  const std::string name_;
  const OrderLambdaType order_lambda_;
  PartialXLambdaType partial_x_lambda_;
  PartialULambdaType partial_u_lambda_;
  FieldVector<ColPartialXLambdaType, dimRangeCols> partial_x_col_lambdas_;
  FieldVector<ColPartialULambdaType, dimRangeCols> partial_u_col_lambdas_;
}; // class LocalLambdaFluxFunction


} // namespace Functions
} // namespace XT
} // namespace Dune
#endif
#endif // DUNE_XT_FUNCTIONS_LAMBDA_LOCAL_FLUX_FUNCTION_HH
