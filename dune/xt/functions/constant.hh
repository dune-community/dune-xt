// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2018)
//   Kirsten Weber   (2013)
//   René Fritze     (2013 - 2018)
//   Tim Keil        (2018)
//   Tobias Leibner  (2014 - 2015, 2017)

#ifndef DUNE_XT_FUNCTIONS_CONSTANT_HH
#define DUNE_XT_FUNCTIONS_CONSTANT_HH

#include <dune/xt/common/configuration.hh>

#include <dune/xt/functions/interfaces/function.hh>

namespace Dune {
namespace XT {
namespace Functions {


template <size_t d, size_t r = 1, size_t rC = 1, class RangeField = double>
class ConstantFunction : public FunctionInterface<d, r, rC, RangeField>
{
  using BaseType = FunctionInterface<d, r, rC, RangeField>;
  using ThisType = ConstantFunction<d, r, rC, RangeField>;

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

  ConstantFunction(const RangeReturnType constant, const std::string name_in = static_id())
    : constant_(constant)
    , name_(name_in)
  {}

#if !DUNE_XT_WITH_PYTHON_BINDINGS
  ConstantFunction(const ThisType& other) = default;
#endif

  int order(const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    return 0;
  }

  RangeReturnType evaluate(const DomainType& /*point_in_global_coordinates*/,
                           const Common::Parameter& /*param*/ = {}) const override final
  {
    return constant_;
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

  const RangeReturnType constant_;
  const std::string name_;
};


template <class Element, size_t rangeDim = 1, size_t rangeDimCols = 1, class RangeField = double>
class ConstantGridFunction : public GridFunctionInterface<Element, rangeDim, rangeDimCols, RangeField>
{
  using BaseType = GridFunctionInterface<Element, rangeDim, rangeDimCols, RangeField>;

public:
  using typename BaseType::LocalFunctionType;

  ConstantGridFunction(const typename LocalFunctionType::RangeReturnType constant,
                       const std::string name_in = static_id())
    : constant_function_(constant, name_in)
    , constant_grid_function_(constant_function_)
  {}

  static std::string static_id()
  {
    return "dune.xt.functions.constantgridfunction";
  }

  virtual std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return constant_grid_function_.local_function();
  }

  virtual std::string name() const override final
  {
    return constant_function_.name();
  }

private:
  ConstantFunction<BaseType::domain_dim, rangeDim, rangeDimCols, RangeField> constant_function_;
  FunctionAsGridFunctionWrapper<Element, rangeDim, rangeDimCols, RangeField> constant_grid_function_;
};


} // namespace Functions
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_FUNCTIONS_CONSTANT_HH
