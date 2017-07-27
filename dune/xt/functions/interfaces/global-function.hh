// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_FUNCTIONS_INTERFACES_GLOBAL_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_INTERFACES_GLOBAL_FUNCTION_HH

#include "localizable-function.hh"

namespace Dune {
namespace XT {
namespace Functions {


template <class OtherEntityImp, class GlobalFunctionImp>
class TransferredGlobalFunction;


/**
 * base class for global scalar-, vector- or matrix-valued valued functions that provides automatic local functions via
 * LocalizableFunctionInterface
 */
template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class GlobalFunctionInterface
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      BaseType;
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;

public:
  typedef typename BaseType::LocalfunctionType LocalfunctionType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  virtual ~GlobalFunctionInterface()
  {
  }

  virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const = 0;

  virtual void evaluate(const DomainType& /*xx*/, RangeType& /*ret*/, const Common::Parameter& /*mu*/ = {}) const = 0;

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& /*ret*/, const Common::Parameter& /*mu*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "You have to implement it if you intend to use it!");
  }

  virtual RangeType evaluate(const DomainType& xx, const Common::Parameter& mu = {}) const
  {
    RangeType ret;
    evaluate(xx, ret, mu);
    return ret;
  }

  virtual JacobianRangeType jacobian(const DomainType& xx, const Common::Parameter& mu = {}) const
  {
    JacobianRangeType ret;
    jacobian(xx, ret, mu);
    return ret;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityImp& entity) const override final
  {
    return Common::make_unique<Localfunction>(entity, *this);
  }

  virtual std::string type() const override
  {
    return "globalfunction";
  }

  virtual std::string name() const override
  {
    return "globalfunction";
  }

private:
  class Localfunction : public LocalfunctionType
  {
  public:
    Localfunction(const EntityImp& entity_in, const ThisType& global_function)
      : LocalfunctionType(entity_in)
      , geometry_(entity_in.geometry())
      , global_function_(global_function)
    {
    }

    virtual ~Localfunction()
    {
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret, const Common::Parameter& mu = {}) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.evaluate(xx_global, ret, mu);
    }

    virtual void
    jacobian(const DomainType& xx, JacobianRangeType& ret, const Common::Parameter& mu = {}) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.jacobian(xx_global, ret, mu);
    }

    virtual size_t order(const Common::Parameter& /*mu*/ = {}) const override final
    {
      return global_function_.order();
    }

  private:
    const typename EntityImp::Geometry geometry_;
    const ThisType& global_function_;
  }; // class Localfunction

public:
  template <class OtherEntityImp>
  struct Transfer
  {
    typedef TransferredGlobalFunction<OtherEntityImp, ThisType> Type;
  };

  template <class OtherEntityImp>
  typename Transfer<OtherEntityImp>::Type transfer() const
  {
    return typename Transfer<OtherEntityImp>::Type(*this);
  }
}; // class GlobalFunctionInterface


template <class OtherEntityImp, class GlobalFunctionImp>
class TransferredGlobalFunction : public GlobalFunctionInterface<OtherEntityImp,
                                                                 typename GlobalFunctionImp::DomainFieldType,
                                                                 GlobalFunctionImp::dimDomain,
                                                                 typename GlobalFunctionImp::RangeFieldType,
                                                                 GlobalFunctionImp::dimRange,
                                                                 GlobalFunctionImp::dimRangeCols>
{
  typedef GlobalFunctionInterface<OtherEntityImp,
                                  typename GlobalFunctionImp::DomainFieldType,
                                  GlobalFunctionImp::dimDomain,
                                  typename GlobalFunctionImp::RangeFieldType,
                                  GlobalFunctionImp::dimRange,
                                  GlobalFunctionImp::dimRangeCols>
      BaseType;

public:
  TransferredGlobalFunction(const GlobalFunctionImp& function)
    : function_(function)
  {
  }

  virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const
  {
    return function_.order();
  }

  virtual void evaluate(const typename BaseType::DomainType& x,
                        typename BaseType::RangeType& ret,
                        const Common::Parameter& mu = {}) const
  {
    function_.evaluate(x, ret, mu);
  }

private:
  const GlobalFunctionImp& function_;
}; // class TransferredGlobalFunction


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_INTERFACES_GLOBAL_FUNCTION_HH
