// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2019)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_FUNCTIONS_BASE_DERIVATIVES_OF_GRID_FUNCTIONS_HH
#define DUNE_XT_FUNCTIONS_BASE_DERIVATIVES_OF_GRID_FUNCTIONS_HH

#include <dune/xt/common/memory.hh>

#include <dune/xt/functions/grid-function.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/type_traits.hh>

#include "derivatives-of-element-functions.hh"


namespace Dune {
namespace XT {
namespace Functions {


template <class GridFunctionType, DerivativeType derivative>
class DerivativeGridFunction
  : public GridFunctionInterface<typename internal::DerivativeElementFunctionHelper<GridFunctionType, derivative>::E,
                                 internal::DerivativeElementFunctionHelper<GridFunctionType, derivative>::r,
                                 internal::DerivativeElementFunctionHelper<GridFunctionType, derivative>::rC,
                                 typename internal::DerivativeElementFunctionHelper<GridFunctionType, derivative>::R>
{
  static_assert(is_grid_function<GridFunctionType>::value, "");

  using BaseType =
      GridFunctionInterface<typename internal::DerivativeElementFunctionHelper<GridFunctionType, derivative>::E,
                            internal::DerivativeElementFunctionHelper<GridFunctionType, derivative>::r,
                            internal::DerivativeElementFunctionHelper<GridFunctionType, derivative>::rC,
                            typename internal::DerivativeElementFunctionHelper<GridFunctionType, derivative>::R>;

  using Select = internal::DerivativeElementFunctionHelper<GridFunctionType, derivative>;

  static const constexpr size_t r_ = GridFunctionType::r;
  static const constexpr size_t rC_ = GridFunctionType::rC;

public:
  using BaseType::r;
  using BaseType::rC;
  using typename BaseType::E;
  using typename BaseType::LocalFunctionType;
  using typename BaseType::R;

  DerivativeGridFunction(GridFunction<E, r_, rC_, R> grid_function, const std::string& nm = "")
    : BaseType()
    , grid_function_(grid_function)
    , name_(nm.empty() ? "dune.xt.functions.derivativegridfunction" : nm)
  {}

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<DerivativeElementFunction<typename GridFunctionType::LocalFunctionType, derivative>>(
        grid_function_.local_function());
  }

  std::string name() const override final
  {
    return name_;
  }

private:
  GridFunction<E, r_, rC_, R> grid_function_;
  const std::string name_;
}; // class DerivativeGridFunction


template <class GridFunctionType>
class DivergenceGridFunction : public DerivativeGridFunction<GridFunctionType, DerivativeType::divergence>
{
  using BaseType = DerivativeGridFunction<GridFunctionType, DerivativeType::divergence>;

public:
  template <class... Args>
  explicit DivergenceGridFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
}; // class DivergenceGridFunction


template <class GridFunctionType>
class GradientGridFunction : public DerivativeGridFunction<GridFunctionType, DerivativeType::gradient>
{
  using BaseType = DerivativeGridFunction<GridFunctionType, DerivativeType::gradient>;

public:
  template <class... Args>
  explicit GradientGridFunction(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
}; // class GradientGridFunction


} // namespace Functions
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_FUNCTIONS_DERIVED_HH
