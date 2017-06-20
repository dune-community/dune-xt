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

#ifndef DUNE_XT_FUNCTIONS_LAMBDA_GLOBAL_FLUX_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_LAMBDA_GLOBAL_FLUX_FUNCTION_HH

#include <functional>

#include <dune/xt/common/memory.hh>

#include <dune/xt/functions/interfaces.hh>

namespace Dune {
namespace XT {
namespace Functions {


/**
 * Global-valued flux function you can pass a lambda expression to that gets evaluated
 * \example LambdaType lambda([](DomainType x, StateRangeType u) { return u;}, 1 );
 */
template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class U_,
          size_t state_derivative_order,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class GlobalLambdaFluxFunction : public GlobalFluxFunctionInterface<EntityImp,
                                                                    DomainFieldImp,
                                                                    domainDim,
                                                                    U_,
                                                                    state_derivative_order,
                                                                    RangeFieldImp,
                                                                    rangeDim,
                                                                    rangeDimCols>
{
  typedef GlobalFluxFunctionInterface<EntityImp,
                                      DomainFieldImp,
                                      domainDim,
                                      U_,
                                      state_derivative_order,
                                      RangeFieldImp,
                                      rangeDim,
                                      rangeDimCols>
      BaseType;

public:
  using typename BaseType::DomainType;
  using typename BaseType::StateRangeType;
  using typename BaseType::RangeType;

private:
  typedef std::function<RangeType(DomainType, StateRangeType)> LambdaType;

public:
  GlobalLambdaFluxFunction(LambdaType lambda,
                           const size_t order_in,
                           const Common::ParameterType& param_type,
                           const std::string nm = "globallambdafunction")
    : lambda_(lambda)
    , order_(order_in)
    , param_type_(param_type)
    , name_(nm)
  {
  }

  size_t order() const override final
  {
    return order_;
  }

  void evaluate(const DomainType& xx,
                const StateRangeType& uu,
                RangeType& ret,
                const Common::Parameter& /*mu*/ = Common::Parameter()) const override final
  {
    ret = lambda_(xx, uu);
  }

  RangeType evaluate(const DomainType& xx,
                     const StateRangeType& uu,
                     const Common::Parameter& /*mu*/ = Common::Parameter()) const override final
  {
    return lambda_(xx, uu);
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
  const size_t order_;
  const Common::ParameterType param_type_;
  const std::string name_;
};


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_LAMBDA_GLOBAL_FLUX_FUNCTION_HH
