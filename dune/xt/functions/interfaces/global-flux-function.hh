// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_FUNCTIONS_INTERFACES_GLOBAL_FLUX_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_INTERFACES_GLOBAL_FLUX_FUNCTION_HH

#include "localizable-function.hh"

namespace Dune {
namespace XT {
namespace Functions {


template <class OtherEntityImp, class GlobalFluxFunctionImp>
class TransferredGlobalFluxFunction;


/**
 * base class for global functions that provides automatic local functions via
 * LocalizableFunctionInterface
 */
template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class U_,
          size_t state_derivative_order,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class GlobalFluxFunctionInterface : public LocalizableFluxFunctionInterface<EntityImp,
                                                                            DomainFieldImp,
                                                                            domainDim,
                                                                            U_,
                                                                            state_derivative_order,
                                                                            RangeFieldImp,
                                                                            rangeDim,
                                                                            rangeDimCols>
{
  typedef LocalizableFluxFunctionInterface<EntityImp,
                                           DomainFieldImp,
                                           domainDim,
                                           U_,
                                           state_derivative_order,
                                           RangeFieldImp,
                                           rangeDim,
                                           rangeDimCols>
      BaseType;
  typedef GlobalFluxFunctionInterface<EntityImp,
                                      DomainFieldImp,
                                      domainDim,
                                      U_,
                                      state_derivative_order,
                                      RangeFieldImp,
                                      rangeDim,
                                      rangeDimCols>
      ThisType;

public:
  typedef typename BaseType::LocalfunctionType LocalfunctionType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::StateType StateType;
  typedef typename BaseType::StateRangeType StateRangeType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename LocalfunctionType::JacobianWrtXRangeType JacobianWrtXRangeType;
  typedef typename LocalfunctionType::JacobianWrtURangeType JacobianWrtURangeType;

  virtual ~GlobalFluxFunctionInterface()
  {
  }

  virtual size_t order(const Common::Parameter& /*mu*/) const = 0;

  virtual void evaluate(const DomainType& /*x*/,
                        const StateRangeType& /*u*/,
                        RangeType& /*ret*/,
                        const Common::Parameter& /*mu*/ = Common::Parameter()) const = 0;

  virtual void jacobian_wrt_x(const DomainType& /*x*/,
                              const StateRangeType& /*u*/,
                              JacobianWrtXRangeType& /*ret*/,
                              const Common::Parameter& /*mu*/ = Common::Parameter()) const
  {
    DUNE_THROW(NotImplemented, "This does not make sense yet for matrix-valued functions!");
  }

  virtual void jacobian_wrt_u(const DomainType& /*x*/,
                              const StateRangeType& /*u*/,
                              JacobianWrtURangeType& /*ret*/,
                              const Common::Parameter& /*mu*/ = Common::Parameter()) const
  {
    DUNE_THROW(NotImplemented, "This does not make sense yet for matrix-valued functions!");
  }

  virtual RangeType
  evaluate(const DomainType& xx, const StateRangeType& uu, const Common::Parameter& mu = Common::Parameter()) const
  {
    RangeType ret;
    evaluate(xx, uu, ret, mu);
    return ret;
  }

  virtual JacobianWrtXRangeType jacobian_wrt_x(const DomainType& xx,
                                               const StateRangeType& uu,
                                               const Common::Parameter& mu = Common::Parameter()) const
  {
    JacobianWrtXRangeType ret;
    jacobian_wrt_x(xx, uu, ret, mu);
    return ret;
  }

  virtual JacobianWrtURangeType jacobian_wrt_u(const DomainType& xx,
                                               const StateRangeType& uu,
                                               const Common::Parameter& mu = Common::Parameter()) const
  {
    JacobianWrtURangeType ret;
    jacobian_wrt_u(xx, uu, ret, mu);
    return ret;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityImp& entity) const override final
  {
    return Common::make_unique<Localfunction>(entity, *this);
  }

  virtual std::string type() const override
  {
    return "globalfluxfunction";
  }

  virtual std::string name() const override
  {
    return "globalfluxfunction";
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

    virtual size_t order(const Common::Parameter& mu) const override final
    {
      return global_function_.order(mu);
    }

    virtual void evaluate(const DomainType& xx,
                          const StateRangeType& uu,
                          RangeType& ret,
                          const Common::Parameter& mu = Common::Parameter()) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.evaluate(xx_global, uu, ret, mu);
    }

    virtual void jacobian_wrt_x(const DomainType& xx,
                                const StateRangeType& uu,
                                JacobianWrtXRangeType& ret,
                                const Common::Parameter& mu = Common::Parameter()) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.jacobian_wrt_x(xx_global, uu, ret, mu);
    }

    virtual void jacobian_wrt_u(const DomainType& xx,
                                const StateRangeType& uu,
                                JacobianWrtURangeType& ret,
                                const Common::Parameter& mu = Common::Parameter()) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.jacobian_wrt_u(xx_global, uu, ret, mu);
    }

  private:
    const typename EntityImp::Geometry geometry_;
    const ThisType& global_function_;
  }; // class Localfunction

public:
  template <class OtherEntityImp>
  struct Transfer
  {
    typedef TransferredGlobalFluxFunction<OtherEntityImp, ThisType> Type;
  };

  template <class OtherEntityImp>
  typename Transfer<OtherEntityImp>::Type transfer() const
  {
    return typename Transfer<OtherEntityImp>::Type(*this);
  }
}; // class GlobalFluxFunctionInterface


template <class OtherEntityImp, class GlobalFluxFunctionImp>
class TransferredGlobalFluxFunction
    : public GlobalFluxFunctionInterface<OtherEntityImp,
                                         typename GlobalFluxFunctionImp::DomainFieldType,
                                         GlobalFluxFunctionImp::dimDomain,
                                         typename GlobalFluxFunctionImp::StateType,
                                         GlobalFluxFunctionImp::state_derivative_order,
                                         typename GlobalFluxFunctionImp::RangeFieldType,
                                         GlobalFluxFunctionImp::dimRange,
                                         GlobalFluxFunctionImp::dimRangeCols>
{
  typedef GlobalFluxFunctionInterface<OtherEntityImp,
                                      typename GlobalFluxFunctionImp::DomainFieldType,
                                      GlobalFluxFunctionImp::dimDomain,
                                      typename GlobalFluxFunctionImp::StateType,
                                      GlobalFluxFunctionImp::state_derivative_order,
                                      typename GlobalFluxFunctionImp::RangeFieldType,
                                      GlobalFluxFunctionImp::dimRange,
                                      GlobalFluxFunctionImp::dimRangeCols>
      BaseType;

public:
  TransferredGlobalFluxFunction(const GlobalFluxFunctionImp& function)
    : function_(function)
  {
  }

  virtual size_t order() const
  {
    return function_.order();
  }

  virtual void evaluate(const typename BaseType::DomainType& x,
                        const typename BaseType::StateRangeType& u,
                        typename BaseType::RangeType& ret,
                        const Common::Parameter& mu = Common::Parameter()) const
  {
    function_.evaluate(x, u, ret, mu);
  }

private:
  const GlobalFluxFunctionImp& function_;
}; // class TransferredFluxGlobalFunction


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_INTERFACES_GLOBAL_FLUX_FUNCTION_HH
