// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2019)
//   Kirsten Weber   (2013)
//   René Fritze     (2013 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2014 - 2015, 2017, 2019 - 2020)

#ifndef DUNE_XT_FUNCTIONS_CONSTANT_HH
#define DUNE_XT_FUNCTIONS_CONSTANT_HH

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/string.hh>

#include <dune/xt/functions/interfaces/function.hh>
#include <dune/xt/functions/interfaces/flux-function.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

namespace Dune {
namespace XT {
namespace Functions {


template <size_t d, size_t r = 1, size_t rC = 1, class R = double>
class ConstantFunction : public FunctionInterface<d, r, rC, R>
{
  using BaseType = FunctionInterface<d, r, rC, R>;
  using ThisType = ConstantFunction;

public:
  using DerivativeRangeReturnType = typename BaseType::DerivativeRangeReturnType;
  using DomainType = typename BaseType::DomainType;
  using RangeReturnType = typename BaseType::RangeReturnType;
  using SingleDerivativeRangeReturnType = typename BaseType::SingleDerivativeRangeReturnType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".constant";
  }

  static Common::Configuration defaults()
  {
    Common::Configuration config;
    config["type"] = static_id();
    config["value"] = "[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]";
    config["name"] = static_id();
    return config;
  } // ... defaults(...)

  explicit ConstantFunction(const RangeReturnType& value, const std::string nm = "")
    : value_(value)
    , name_(nm.empty() ? ((r == 1 && rC == 1) ? std::string("ConstantFunction(" + Common::to_string(value) + ")")
                                              : "ConstantFunction")
                       : nm)
  {}

  ConstantFunction(const ThisType&) = default;

  ConstantFunction(ThisType&) = default;

private:
  ThisType* copy_as_function_impl() const override
  {
    return new ThisType(*this);
  }

public:
  std::unique_ptr<ThisType> copy_as_function() const
  {
    return std::unique_ptr<ThisType>(this->copy_as_function_impl());
  }

  int order(const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    return 0;
  }

  RangeReturnType evaluate(const DomainType& /*point_in_global_coordinates*/,
                           const Common::Parameter& /*param*/ = {}) const override final
  {
    return value_;
  }

  DerivativeRangeReturnType jacobian(const DomainType& /*point_in_global_coordinates*/,
                                     const Common::Parameter& /*param*/ = {}) const override final
  {
    return DerivativeRangeReturnType(); // defaults to 0
  }

  std::string name() const override final
  {
    return name_;
  }

  const RangeReturnType value_;
  const std::string name_;
}; // class ConstantFunction


template <class E, size_t r = 1, size_t rC = 1, class R = double>
class ConstantGridFunction : public FunctionAsGridFunctionWrapper<E, r, rC, R>
{
  using ThisType = ConstantGridFunction;
  using BaseType = FunctionAsGridFunctionWrapper<E, r, rC, R>;

public:
  using typename BaseType::LocalFunctionType;

  ConstantGridFunction(const typename LocalFunctionType::RangeReturnType value, const std::string nm = "")
    : BaseType(ConstantFunction<E::dimension, r, rC, R>(
        value,
        nm.empty() ? ((r == 1 && rC == 1) ? std::string("ConstantGridFunction(" + Common::to_string(value) + ")")
                                          : "ConstantGridFunction")
                   : nm))
  {}
  std::unique_ptr<ThisType> copy_as_grid_function() const
  {
    return std::unique_ptr<ThisType>(this->copy_as_grid_function_impl());
  }

private:
  ThisType* copy_as_grid_function_impl() const override
  {
    return new ThisType(*this);
  }

}; // class ConstantGridFunction


template <class E, size_t stateDim, size_t r = 1, size_t rC = 1, class R = double>
class ConstantFluxFunction : public FluxFunctionInterface<E, stateDim, r, rC, R>
{
  using BaseType = FluxFunctionInterface<E, stateDim, r, rC, R>;

public:
  using typename BaseType::LocalFunctionType;

  ConstantFluxFunction(const typename LocalFunctionType::RangeReturnType constant,
                       const std::string name_in = static_id())
    : constant_function_(constant, name_in)
    , constant_flux_function_(constant_function_)
  {}

  static std::string static_id()
  {
    return "dune.xt.functions.constantfluxfunction";
  }

  bool x_dependent() const override final
  {
    return false;
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return constant_flux_function_.local_function();
  }

  std::string name() const override final
  {
    return constant_function_.name();
  }

private:
  ConstantFunction<stateDim, r, rC, R> constant_function_;
  StateFunctionAsFluxFunctionWrapper<E, stateDim, r, rC, R> constant_flux_function_;
}; // class ConstantGridFunction


} // namespace Functions
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_FUNCTIONS_CONSTANT_HH
