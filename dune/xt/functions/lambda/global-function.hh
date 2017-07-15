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

private:
  typedef std::function<RangeType(DomainType, XT::Common::Parameter)> LambdaType;

public:
  GlobalLambdaFunction(LambdaType lambda, const size_t order_in, const std::string nm = "globallambdafunction")
    : lambda_(lambda)
    , order_(order_in)
    , name_(nm)
  {
  }

  size_t order() const override final
  {
    return order_;
  }

  void
  evaluate(const DomainType& xx, RangeType& ret, const Common::Parameter& mu = Common::Parameter()) const override final
  {
    ret = lambda_(xx, mu);
  }

  RangeType evaluate(const DomainType& xx, const Common::Parameter& mu = Common::Parameter()) const override final
  {
    return lambda_(xx, mu);
  }

  std::string type() const override final
  {
    return "globallambdafunction";
  }

  std::string name() const override final
  {
    return name_;
  }

private:
  LambdaType lambda_;
  size_t order_;
  std::string name_;
};


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_LAMBDA_GLOBAL_FUNCTION_HH
