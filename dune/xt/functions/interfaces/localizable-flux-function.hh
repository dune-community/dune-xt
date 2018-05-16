// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017 - 2018)

#ifndef DUNE_XT_FUNCTIONS_INTERFACES_LOCALIZABLE_FLUX_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_INTERFACES_LOCALIZABLE_FLUX_FUNCTION_HH

#include <dune/common/typetraits.hh>

#include <dune/xt/common/parameter.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/type_traits.hh>

#include "local-flux-function.hh"

namespace Dune {
namespace XT {
namespace Functions {


template <class E_, class D_, size_t d_, class U_, size_t s_, class R_, size_t r_, size_t rC_ = 1>
class LocalizableFluxFunctionInterface : public Common::ParametricInterface
{
  static_assert(is_localizable_function<U_>::value, "");

public:
  typedef E_ E;
  typedef E_ EntityType;
  typedef D_ D;
  typedef D_ DomainFieldType;
  static const constexpr size_t d = d_;
  static const constexpr size_t dimDomain = d_;
  typedef U_ U;
  typedef U_ StateType;
  static const constexpr size_t s = s_;
  static const constexpr size_t stateDerivativeOrder = s_;
  typedef R_ R;
  typedef R_ RangeFieldType;
  static const constexpr size_t r = r_;
  static const constexpr size_t dimRange = r_;
  static const constexpr size_t rC = rC_;
  static const constexpr size_t dimRangeCols = rC_;
  typedef LocalFluxFunctionInterface<E, D, d, typename U::LocalfunctionType, s, R, r, rC> LocalfunctionType;

  typedef typename LocalfunctionType::DomainType DomainType;
  typedef typename LocalfunctionType::RangeType RangeType;
  typedef typename LocalfunctionType::StateRangeType StateRangeType;
  typedef typename LocalfunctionType::PartialXRangeType PartialXRangeType;
  typedef typename LocalfunctionType::PartialURangeType PartialURangeType;
  typedef typename LocalfunctionType::ColRangeType ColRangeType;
  typedef typename LocalfunctionType::ColPartialXRangeType ColPartialXRangeType;
  typedef typename LocalfunctionType::ColPartialURangeType ColPartialURangeType;

  static const bool available = false;
  virtual ~LocalizableFluxFunctionInterface() = default;

  static std::string static_id()
  {
    return "fluxfunction";
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& /*entity*/) const = 0;

  virtual bool is_affine() const
  {
    return false;
  }

  /**
   * \name ´´These methods should be implemented in order to identify the function.''
   * \{
   */
  virtual std::string type() const
  {
    return "fluxfunction";
  }

  virtual std::string name() const
  {
    return "fluxfunction";
  }
  /**
    * \}
    */
}; // class LocalizableFluxFunctionInterface


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_INTERFACES_LOCALIZABLE_FLUX_FUNCTION_HH
