// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2017)

#ifndef DUNE_XT_FUNCTIONS_INTERFACES_GLOBAL_FLUX_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_INTERFACES_GLOBAL_FLUX_FUNCTION_HH

#include "localizable-function.hh"

namespace Dune {
namespace XT {
namespace Functions {


template <class OtherEntityImp, class GlobalFluxFunctionImp>
class TransferredGlobalFluxFunction;


/**
 * Base class for flux functions that can be evaluated in global coordinates without an entity. Provides automatic local
 * functions via LocalizableFluxFunctionInterface. Used to model functions that are continuous in x and u.
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
  typedef typename LocalfunctionType::PartialXRangeType PartialXRangeType;
  typedef typename LocalfunctionType::PartialURangeType PartialURangeType;
  typedef typename BaseType::ColRangeType ColRangeType;
  typedef typename LocalfunctionType::ColPartialXRangeType ColPartialXRangeType;
  typedef typename LocalfunctionType::ColPartialURangeType ColPartialURangeType;

  virtual ~GlobalFluxFunctionInterface()
  {
  }

  virtual size_t order(const Common::Parameter& /*mu*/ = {}) const = 0;

  virtual void evaluate(const DomainType& /*x*/,
                        const StateRangeType& /*u*/,
                        RangeType& /*ret*/,
                        const Common::Parameter& /*mu*/ = {}) const = 0;

  virtual void evaluate_col(const size_t /*col*/,
                            const DomainType& /*x*/,
                            const StateRangeType& /*u*/,
                            ColRangeType& /*ret*/,
                            const Common::Parameter& /*mu*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  virtual void partial_x(const DomainType& /*x*/,
                         const StateRangeType& /*u*/,
                         PartialXRangeType& /*ret*/,
                         const Common::Parameter& /*mu*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  virtual void partial_x_col(const size_t /*col*/,
                             const DomainType& /*x*/,
                             const StateRangeType& /*u*/,
                             ColPartialXRangeType& /*ret*/,
                             const Common::Parameter& /*mu*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  virtual void partial_u(const DomainType& /*x*/,
                         const StateRangeType& /*u*/,
                         PartialURangeType& /*ret*/,
                         const Common::Parameter& /*mu*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  virtual void partial_u_col(const size_t /*col*/,
                             const DomainType& /*x*/,
                             const StateRangeType& /*u*/,
                             ColPartialURangeType& /*ret*/,
                             const Common::Parameter& /*mu*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  virtual RangeType evaluate(const DomainType& xx, const StateRangeType& uu, const Common::Parameter& mu = {}) const
  {
    RangeType ret;
    evaluate(xx, uu, ret, mu);
    return ret;
  }

  virtual ColRangeType
  evaluate_col(const size_t col, const DomainType& xx, const StateRangeType& uu, const Common::Parameter& mu = {}) const
  {
    ColRangeType ret;
    evaluate_col(col, xx, uu, ret, mu);
    return ret;
  }

  virtual PartialXRangeType
  partial_x(const DomainType& xx, const StateRangeType& uu, const Common::Parameter& mu = {}) const
  {
    PartialXRangeType ret;
    partial_x(xx, uu, ret, mu);
    return ret;
  }

  virtual ColPartialXRangeType partial_x_col(const size_t col,
                                             const DomainType& xx,
                                             const StateRangeType& uu,
                                             const Common::Parameter& mu = {}) const
  {
    ColPartialXRangeType ret;
    partial_x_col(col, xx, uu, ret, mu);
    return ret;
  }

  virtual PartialURangeType
  partial_u(const DomainType& xx, const StateRangeType& uu, const Common::Parameter& mu = {}) const
  {
    PartialURangeType ret;
    partial_u(xx, uu, ret, mu);
    return ret;
  }

  virtual ColPartialURangeType partial_u_col(const size_t col,
                                             const DomainType& xx,
                                             const StateRangeType& uu,
                                             const Common::Parameter& mu = {}) const
  {
    ColPartialURangeType ret;
    partial_u_col(col, xx, uu, ret, mu);
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
                          const Common::Parameter& mu = {}) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.evaluate(xx_global, uu, ret, mu);
    }

    virtual void evaluate_col(const size_t col,
                              const DomainType& xx,
                              const StateRangeType& uu,
                              ColRangeType& ret,
                              const Common::Parameter& mu = {}) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.evaluate_col(col, xx_global, uu, ret, mu);
    }


    virtual void partial_x(const DomainType& xx,
                           const StateRangeType& uu,
                           PartialXRangeType& ret,
                           const Common::Parameter& mu = {}) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.partial_x(xx_global, uu, ret, mu);
    }

    virtual void partial_u(const DomainType& xx,
                           const StateRangeType& uu,
                           PartialURangeType& ret,
                           const Common::Parameter& mu = {}) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.partial_u(xx_global, uu, ret, mu);
    }

    virtual void partial_u_col(const size_t col,
                               const DomainType& xx,
                               const StateRangeType& uu,
                               ColPartialURangeType& ret,
                               const Common::Parameter& mu = {}) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.partial_u_col(col, xx_global, uu, ret, mu);
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

  virtual size_t order(const XT::Common::Parameter& mu = {}) const
  {
    return function_.order(mu);
  }

  virtual void evaluate(const typename BaseType::DomainType& x,
                        const typename BaseType::StateRangeType& u,
                        typename BaseType::RangeType& ret,
                        const Common::Parameter& mu = {}) const
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
