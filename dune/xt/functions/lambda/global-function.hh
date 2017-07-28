// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2013 - 2016)
//   Sven Kaulmann   (2013)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_FUNCTIONS_LAMBDA_GLOBAL_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_LAMBDA_GLOBAL_FUNCTION_HH

#include <functional>

#include <dune/xt/common/memory.hh>

#include <dune/xt/functions/interfaces.hh>

namespace Dune {
namespace XT {
namespace Functions {


/**
 * Global-valued function you can pass a lambda expression to that gets evaluated
 * \example LambdaType lambda([](DomainType x) { return x;}, 1 );
 */
template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class GlobalLambdaFunction
    : public GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> BaseType;

public:
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

private:
  typedef std::function<RangeType(DomainType, XT::Common::Parameter)> LambdaType;
  typedef std::function<size_t(const Common::Parameter&)> OrderLambdaType;
  typedef std::function<JacobianRangeType(DomainType, XT::Common::Parameter)> JacobianLambdaType;

public:
  GlobalLambdaFunction(LambdaType lambda,
                       const size_t order_in,
                       const Common::ParameterType& param_type = {},
                       const std::string nm = "globallambdafunction",
                       JacobianLambdaType jacobian_lambda =
                           [](const DomainType&, const Common::Parameter&) {
                             DUNE_THROW(NotImplemented,
                                        "You need to provide a lambda for the jacobian if you want to use it!");
                             return JacobianRangeType();
                           })
    : lambda_(lambda)
    , order_lambda_([=](const Common::Parameter&) { return order_in; })
    , param_type_(param_type)
    , name_(nm)
    , jacobian_lambda_(jacobian_lambda)
  {
  }

  GlobalLambdaFunction(LambdaType lambda,
                       OrderLambdaType order_lambda,
                       const Common::ParameterType& param_type = {},
                       const std::string nm = "globallambdafunction",
                       JacobianLambdaType jacobian_lambda =
                           [](const DomainType&, const Common::Parameter&) {
                             DUNE_THROW(NotImplemented,
                                        "You need to provide a lambda for the jacobian if you want to use it!");
                             return JacobianRangeType();
                           })
    : lambda_(lambda)
    , order_lambda_(order_lambda)
    , param_type_(param_type)
    , name_(nm)
    , jacobian_lambda_(jacobian_lambda)
  {
  }

  virtual size_t order(const Common::Parameter& mu = {}) const override final
  {
    return order_lambda_(mu);
  }

  virtual void evaluate(const DomainType& xx, RangeType& ret, const Common::Parameter& mu = {}) const override final
  {
    auto parsed_mu = this->parse_and_check(mu);
    ret = lambda_(xx, parsed_mu);
  }

  virtual RangeType evaluate(const DomainType& xx, const Common::Parameter& mu = {}) const override final
  {
    auto parsed_mu = this->parse_and_check(mu);
    return lambda_(xx, parsed_mu);
  }

  virtual void
  jacobian(const DomainType& xx, JacobianRangeType& ret, const Common::Parameter& mu = {}) const override final
  {
    auto parsed_mu = this->parse_and_check(mu);
    ret = jacobian_lambda_(xx, parsed_mu);
  }

  virtual std::string type() const override final
  {
    return "globallambdafunction";
  }

  virtual std::string name() const override final
  {
    return name_;
  }

private:
  LambdaType lambda_;
  OrderLambdaType order_lambda_;
  Common::ParameterType param_type_;
  std::string name_;
  JacobianLambdaType jacobian_lambda_;
};


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_LAMBDA_GLOBAL_FUNCTION_HH
