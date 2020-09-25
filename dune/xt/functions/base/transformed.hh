// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2018 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_FUNCTIONS_BASE_TRANSFORMED_HH
#define DUNE_XT_FUNCTIONS_BASE_TRANSFORMED_HH

#include <dune/xt/functions/exceptions.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/type_traits.hh>


namespace Dune {
namespace XT {
namespace Functions {


/**
 * \brief Can be used to apply a transformation to the evaluation of a given function.
 * \todo  In general, it would be nicer to have a concatenation of two functions with matching dimensions.
 *
 *        Consider a given function u: R^d -> R^{d+2}, i.e. modeling the conservative variables [density (scalar),
 *        density*velocity (d-dimensional), energy (scalar)] in the context of the euler equations and one wishes to
 *        obtain a function modeling the primitive variables [density (scalar), velocity (d-dimensional), pressure
 *        (scalar)], the following code would achieve just that (with u_conservative already given):
\code
const double gamma = 1.4;
const auto to_primitive = [&](const auto& conservative_variables) {
  const auto& density = conservative_variables[0];
  DomainType velocity;
  for (size_t ii = 0; ii < d; ++ii)
    velocity = conservative_variables[ii + 1] / density;
  const auto& energy = conservative_variables[m - 1];
  FieldVector<R, d + 2> primitive_variables;
  primitive_variables[0] = density;
  for (size_t ii = 0; ii < d; ++ii)
    primitive_variables[ii + 1] = velocity[ii];
  primitive_variables[m - 1] = (gamma - 1.) * (energy - 0.5 * density * velocity.two_norm2());
  return primitive_variables;
};
auto u_primitive = XT::Functions::make_transformed_function<d + 2, 1, R>(u_conservative, to_primitive);
\endcode
 */
template <class GF, size_t r = GF::r, size_t rC = GF::rC, class R = typename GF::R>
class TransformedGridFunction : public XT::Functions::GridFunctionInterface<typename GF::E, r, rC, R>
{
  static_assert(is_grid_function<GF>::value, "");

  using ThisType = TransformedGridFunction;
  using BaseType = XT::Functions::GridFunctionInterface<typename GF::E, r, rC, R>;

  class TransformedLocalFunction : public XT::Functions::ElementFunctionInterface<typename GF::E, r, rC, R>
  {
    using BaseType = XT::Functions::ElementFunctionInterface<typename GF::E, r, rC, R>;
    using UntransformedLocalFunctionType = typename GF::LocalFunctionType;

  public:
    using UntransformedRangeType = typename UntransformedLocalFunctionType::RangeType;
    using typename BaseType::DerivativeRangeReturnType;
    using typename BaseType::DerivativeRangeType;
    using typename BaseType::DomainType;
    using typename BaseType::ElementType;
    using typename BaseType::RangeReturnType;
    using typename BaseType::RangeType;
    using Transformation = std::function<RangeType(const UntransformedRangeType&)>;

    TransformedLocalFunction(const GF& function, const Transformation& transformation)
      : BaseType(function.parameter_type())
      , function_(function.copy_as_grid_function())
      , local_function_(function_->local_function())
      , transformation_(transformation)
    {}

  protected:
    void post_bind(const ElementType& element) override final
    {
      local_function_->bind(element);
    }

  public:
    int order(const XT::Common::Parameter& param = {}) const override final
    {
      return local_function_->order(param);
    }

    RangeReturnType evaluate(const DomainType& xx, const XT::Common::Parameter& param = {}) const override final
    {
      return transformation_(local_function_->evaluate(xx, param));
    }

    DerivativeRangeReturnType jacobian(const DomainType& /*xx*/,
                                       const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      DUNE_THROW(NotImplemented, "TransformedLocalFunction does not provide jacobian evaluations (yet)!");
    }

  private:
    std::unique_ptr<GridFunctionInterface<typename GF::E, GF::r, GF::rC, typename GF::R>> function_;
    std::unique_ptr<UntransformedLocalFunctionType> local_function_;
    const Transformation& transformation_;
  }; // class TransformedLocalFunction

public:
  using typename BaseType::ElementType;
  using typename BaseType::LocalFunctionType;
  using UntransformedRangeType = typename TransformedLocalFunction::UntransformedRangeType;
  using TransformedRangeType = typename TransformedLocalFunction::RangeType;

  TransformedGridFunction(const GF& func,
                          std::function<TransformedRangeType(const UntransformedRangeType&)> transformation,
                          const std::string& nm = "")
    : BaseType(func.parameter_type())
    , function_(func.copy_as_grid_function())
    , transformation_(transformation)
    , name_(nm.empty() ? "transformed " + function_->name() : nm)
  {}

  TransformedGridFunction(const ThisType& other)
    : BaseType(other)
    , function_(other.function_->copy_as_grid_function())
    , transformation_(other.transformation_)
    , name_(other.name_)
  {}

  TransformedGridFunction(ThisType&&) = default;


private:
  ThisType* copy_as_grid_function_impl() const override
  {
    return new ThisType(*this);
  }

public:
  std::unique_ptr<ThisType> copy_as_grid_function() const
  {
    return std::unique_ptr<ThisType>(this->copy_as_grid_function_impl());
  }
  std::string name() const override final
  {
    return name_;
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<TransformedLocalFunction>(*function_, transformation_);
  }

private:
  std::unique_ptr<GridFunctionInterface<typename GF::E, GF::r, GF::rC, typename GF::R>> function_;
  const typename TransformedLocalFunction::Transformation transformation_;
  const std::string name_;
}; // class TransformedGridFunction


template <size_t new_r, size_t new_rC, class new_R, class E, size_t r, size_t rC, class R>
TransformedGridFunction<GridFunctionInterface<E, r, rC, R>, new_r, new_rC, new_R> make_transformed_function(
    const GridFunctionInterface<E, r, rC, R>& function,
    std::function<
        typename TransformedGridFunction<GridFunctionInterface<E, r, rC, R>, new_r, new_rC, new_R>::
            TransformedRangeType(
                const typename TransformedGridFunction<GridFunctionInterface<E, r, rC, R>, new_r, new_rC, new_R>::
                    UntransformedRangeType&)> transformation,
    const std::string& name = "xt.functions.transformed")
{
  return TransformedGridFunction<GridFunctionInterface<E, r, rC, R>, new_r, new_rC, new_R>(
      function, transformation, name);
}


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_TRANSFORMED_HH
