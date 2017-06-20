// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_FUNCTIONS_INTERFACES_LOCAL_FLUX_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_INTERFACES_LOCAL_FLUX_FUNCTION_HH

#include <dune/common/typetraits.hh>

#include <dune/xt/common/parameter.hh>
#include <dune/xt/functions/type_traits.hh>

#include "local-functions.hh"

namespace Dune {
namespace XT {
namespace Functions {


template <class E, class D, size_t d, class U, size_t state_derivative_order, class R, size_t r, size_t rC = 1>
class LocalFluxFunctionInterface : public Common::ParametricInterface
{
  static_assert(AlwaysFalse<E>::value, "Not available yet for this state_derivative_order!");
};


template <class E_, class D_, size_t d_, class U_, class R_, size_t r_, size_t rC_>
class LocalFluxFunctionInterface<E_, D_, d_, U_, 0, R_, r_, rC_> : public Common::ParametricInterface
{
  static_assert(is_localfunction<U_>::value, "");

public:
  // domain of physical space
  typedef E_ E;
  typedef E_ EntityType;
  typedef D_ D;
  typedef D_ DomainFieldType;
  static const constexpr size_t d = d_;
  static const constexpr size_t dimDomain = d_;
  typedef Dune::FieldVector<D, d> DomainType;
  // domain of state space
  typedef U_ U;
  typedef U_ StateType;
  typedef typename StateType::RangeType StateRangeType;
  // range
  typedef R_ R;
  typedef R_ RangeFieldType;
  static const constexpr size_t r = r_;
  static const constexpr size_t dimRange = r_;
  static const constexpr size_t rC = rC_;
  static const constexpr size_t dimRangeCols = rC_;
  typedef typename LocalfunctionSetInterface<E, D, d, R, r, rC>::RangeType RangeType;
  typedef typename LocalfunctionSetInterface<E, D, d, R, r, rC>::JacobianRangeType PartialXRangeType;
  typedef typename JacobianRangeTypeSelector<StateType::dimRange, R, r, rC>::type PartialURangeType;

  LocalFluxFunctionInterface(const EntityType& en)
    : entity_(en)
  {
  }

  virtual ~LocalFluxFunctionInterface() = default;

  virtual const EntityType& entity() const
  {
    return entity_;
  }

  virtual size_t order(const Common::Parameter& /*mu*/) const = 0;

  virtual void evaluate(const DomainType& /*x*/,
                        const StateRangeType& /*u*/,
                        RangeType& /*ret*/,
                        const Common::Parameter& /*mu*/ = Common::Parameter()) const = 0;

  virtual void partial_x(const DomainType& /*x*/,
                         const StateRangeType& /*u*/,
                         PartialXRangeType& /*ret*/,
                         const Common::Parameter& /*mu*/ = Common::Parameter()) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  virtual void partial_u(const DomainType& /*x*/,
                         const StateRangeType& /*u*/,
                         PartialURangeType& /*ret*/,
                         const Common::Parameter& /*mu*/ = Common::Parameter()) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  /**
   * \name ´´These methods are provided by the interface.''
   * \{
   **/
  RangeType
  evaluate(const DomainType& x, const StateRangeType& u, const Common::Parameter& mu = Common::Parameter()) const
  {
    RangeType ret(0);
    evaluate(x, u, ret, mu);
    return ret;
  }

  PartialXRangeType
  partial_x(const DomainType& x, const StateRangeType& u, const Common::Parameter& mu = Common::Parameter()) const
  {
    PartialXRangeType ret(0);
    partial_x(x, u, ret, mu);
    return ret;
  }

  PartialURangeType
  partial_u(const DomainType& x, const StateRangeType& u, const Common::Parameter& mu = Common::Parameter()) const
  {
    PartialURangeType ret(0);
    partial_u(x, u, ret, mu);
    return ret;
  }

  /**
   * \}
   **/
private:
  const EntityType& entity_;
}; // class LocalFluxFunctionInterface< ..., 0, ...>


template <class E_, class D_, size_t d_, class U_, class R_, size_t r_, size_t rC_>
class LocalFluxFunctionInterface<E_, D_, d_, U_, 1, R_, r_, rC_> : public Common::ParametricInterface
{
  static_assert(is_localfunction<U_>::value, "");

public:
  // domain of physical space
  typedef E_ E;
  typedef E_ EntityType;
  typedef D_ D;
  typedef D_ DomainFieldType;
  static const constexpr size_t d = d_;
  static const constexpr size_t dimDomain = d_;
  typedef Dune::FieldVector<D, d> DomainType;
  // domain of state space
  typedef U_ U;
  typedef U_ StateType;
  typedef typename StateType::RangeType StateRangeType;
  typedef typename StateType::JacobianRangeType StateJacobianRangeType;
  // range
  typedef R_ R;
  typedef R_ RangeFieldType;
  static const constexpr size_t r = r_;
  static const constexpr size_t dimRange = r_;
  static const constexpr size_t rC = rC_;
  static const constexpr size_t dimRangeCols = rC_;
  typedef typename LocalfunctionSetInterface<E, D, d, R, r, rC>::RangeType RangeType;
  typedef typename LocalfunctionSetInterface<E, D, d, R, r, rC>::JacobianRangeType JacobianRangeType;

  LocalFluxFunctionInterface(const EntityType& en)
    : entity_(en)
  {
  }

  virtual ~LocalFluxFunctionInterface() = default;

  virtual const EntityType& entity() const
  {
    return entity_;
  }

  virtual void evaluate(const DomainType& /*x*/,
                        const StateRangeType& /*u*/,
                        const StateJacobianRangeType& /*grad_u*/,
                        RangeType& /*ret*/,
                        const Common::Parameter& /*mu*/ = Common::Parameter()) const = 0;

  virtual void jacobian(const DomainType& /*x*/,
                        const StateRangeType& /*u*/,
                        const StateJacobianRangeType& /*grad_u*/,
                        JacobianRangeType& /*ret*/,
                        const Common::Parameter& /*mu*/ = Common::Parameter()) const = 0;

  /**
   * \name ´´These methods are provided by the interface.''
   * \{
   **/
  RangeType evaluate(const DomainType& x,
                     const StateType& u,
                     const StateJacobianRangeType& grad_u,
                     const Common::Parameter& mu = Common::Parameter()) const
  {
    RangeType ret(0);
    evaluate(x, u, grad_u, ret, mu);
    return ret;
  }

  JacobianRangeType jacobian(const DomainType& x,
                             const StateType& u,
                             const StateJacobianRangeType& grad_u,
                             const Common::Parameter& mu = Common::Parameter()) const
  {
    JacobianRangeType ret(0);
    jacobian(x, u, grad_u, ret, mu);
    return ret;
  }
  /**
   * \}
   **/
private:
  const EntityType& entity_;
}; // class LocalFluxFunctionInterface< ..., 1, ... >


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_INTERFACES_LOCAL_FLUX_FUNCTION_HH
