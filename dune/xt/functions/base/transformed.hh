// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_FUNCTIONS_BASE_TRANSFORMED_HH
#define DUNE_XT_FUNCTIONS_BASE_TRANSFORMED_HH

#include <dune/xt/functions/type_traits.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>


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
template <class LF, size_t r = LF::r, size_t rC = LF::rC, class R = typename LF::R>
class TransformedGridFunction : public XT::Functions::GridFunctionInterface<typename LF::E, r, rC, R>
{
  static_assert(is_grid_function<LF>::value, "");
  using BaseType = XT::Functions::GridFunctionInterface<typename LF::E, r, rC, R>;

  class TransformedLocalFunction : public XT::Functions::ElementFunctionInterface<typename LF::E, r, rC, R>
  {
    using BaseType = XT::Functions::ElementFunctionInterface<typename LF::E, r, rC, R>;
    using UntransformedLocalFunctionType = typename LF::LocalFunctionType;

  public:
    using UntransformedRangeType = typename UntransformedLocalFunctionType::RangeType;
    using typename BaseType::DomainType;
    using typename BaseType::RangeType;
    using typename BaseType::DerivativeRangeType;
    using typename BaseType::ElementType;
    using Transformation = std::function<RangeType(const UntransformedRangeType&)>;

    TransformedLocalFunction(const LF& function, const Transformation& transformation)
      : BaseType()
      , local_function_(function.local_function())
      , transformation_(transformation)
    {
    }

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

    RangeType evaluate(const DomainType& xx, const XT::Common::Parameter& param = {}) const override final
    {
      return transformation_(local_function_->evaluate(xx, param));
    }

    DerivativeRangeType jacobian(const DomainType& /*xx*/,
                                 const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      DUNE_THROW(NotImplemented, "Yet!");
    }

  private:
    std::unique_ptr<UntransformedLocalFunctionType> local_function_;
    const Transformation& transformation_;
  }; // class TransformedLocalFunction

public:
  using typename BaseType::ElementType;
  using typename BaseType::LocalFunctionType;
  using UntransformedRangeType = typename TransformedLocalFunction::UntransformedRangeType;
  using TransformedRangeType = typename TransformedLocalFunction::RangeType;

  TransformedGridFunction(const LF& f,
                          std::function<TransformedRangeType(const UntransformedRangeType&)> transformation,
                          const std::string& nm = "")
    : function_(f)
    , transformation_(transformation)
    , name_(nm)
  {
  }

  std::string name() const override final
  {
    return name_.empty() ? "transformed " + function_.name() : name_;
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<TransformedLocalFunction>(function_, transformation_);
  }

private:
  const LF& function_;
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
