// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef DUNE_XT_FUNCTIONS_ELEMENTWISE_MINIMUM_HH
#define DUNE_XT_FUNCTIONS_ELEMENTWISE_MINIMUM_HH

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/la/eigen-solver.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/grid-function.hh>
#include <dune/xt/functions/type_traits.hh>

namespace Dune::XT::Functions {
namespace internal {


template <class FunctionType>
class ElementwiseMinimumFunctionHelper
{
  static_assert(is_element_function<FunctionType>::value);

public:
  static const constexpr size_t d = FunctionType::d;
  static const constexpr size_t r = FunctionType::r;
  static const constexpr size_t rC = FunctionType::rC;
  using E = typename FunctionType::E;
  using R = typename FunctionType::R;
  using DomainType = Dune::FieldVector<double, d>;

  template <size_t r_ = r, size_t rC_ = rC, bool anything = true>
  struct dim_switch
  {};

  template <bool anything>
  struct dim_switch<1, 1, anything>
  {
    static R compute(const FunctionType& func, const int order, const XT::Common::Parameter& param)
    {
      // approximate minimum over the element (evaluate at some points)
      double min = std::numeric_limits<double>::max();
      for (auto&& quadrature_point : QuadratureRules<double, d>::rule(func.element().type(), order))
        min = std::min(min, func.evaluate(quadrature_point.position(), param)[0]);
      return min;
    }
  };

  template <size_t r_, bool anything>
  struct dim_switch<r_, r_, anything>
  {
    static R compute(const FunctionType& func, const int order, const XT::Common::Parameter& param)
    {
      // approximate minimum eigenvalue over the element (evaluate at some points)
      double min_EV = std::numeric_limits<double>::max();
      for (auto&& quadrature_point : QuadratureRules<double, d>::rule(func.element().type(), order)) {
        auto value = func.evaluate(quadrature_point.position(), param);
        auto eigen_solver =
            XT::LA::make_eigen_solver(value,
                                      {{"type", XT::LA::EigenSolverOptions<decltype(value)>::types().at(0)},
                                       {"assert_positive_eigenvalues", "1e-15"}});
        min_EV = std::min(min_EV, eigen_solver.min_eigenvalues(1).at(0));
      }
      return min_EV;
    }
  };

public:
  static R compute(const FunctionType& func, const int order, const XT::Common::Parameter& param)
  {
    return dim_switch<>::compute(func, order, param);
  }
}; // class ElementwiseMinimumFunctionHelper


} // namespace internal


/// \todo Consider searching the elements corners!
template <class SomeFunction>
class ElementwiseMinimumFunction : public GridFunctionInterface<typename SomeFunction::E>
{
  static_assert(is_function<SomeFunction>::value || is_grid_function<SomeFunction>::value);

  using ThisType = ElementwiseMinimumFunction;
  using BaseType = GridFunctionInterface<typename SomeFunction::E>;

  static const constexpr size_t r_ = SomeFunction::r;
  static const constexpr size_t rC_ = SomeFunction::rC;

public:
  using BaseType::r;
  using BaseType::rC;
  using typename BaseType::E;
  using typename BaseType::LocalFunctionType;
  using typename BaseType::R;

private:
  static_assert(r == rC);

  class LocalFunction : public XT::Functions::ElementFunctionInterface<E>
  {
    using BaseType = XT::Functions::ElementFunctionInterface<E>;

  public:
    using typename BaseType::DerivativeRangeReturnType;
    using typename BaseType::DomainType;
    using typename BaseType::ElementType;
    using typename BaseType::RangeReturnType;

    LocalFunction(const GridFunctionInterface<E, r_, rC_, R>& some_func, const int search_quadrature_order)
      : BaseType()
      , some_func_(some_func.copy_as_grid_function())
      , some_lf_(some_func_->local_function())
      , search_quadrature_order_(search_quadrature_order)
      , min_(0)
    {}

  protected:
    void post_bind(const ElementType& element) override final
    {
      some_lf_->bind(element);
      min_ = internal::ElementwiseMinimumFunctionHelper<typename SomeFunction::LocalFunctionType>::compute(
          *some_lf_, search_quadrature_order_, {});
    }

  public:
    int order(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return 0;
    }

    RangeReturnType evaluate(const DomainType& /*xx*/, const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return min_;
    }

  private:
    std::unique_ptr<GridFunctionInterface<E, r_, rC_, R>> some_func_;
    std::unique_ptr<typename GridFunctionInterface<E, r_, rC_, R>::LocalFunctionType> some_lf_;
    const int search_quadrature_order_;
    double min_;
  }; // class LocalFunction

public:
  ElementwiseMinimumFunction(GridFunction<E, r_, rC_> some_func,
                             const int search_quadrature_order,
                             const std::string nm = "ElementwiseMinimumFunction")
    : BaseType()
    , some_func_(some_func.copy_as_grid_function())
    , search_quadrature_order_(search_quadrature_order)
    , name_(nm)
  {
    DUNE_THROW_IF(!some_func_->parameter_type().empty(),
                  Exceptions::parameter_error,
                  "Not available for parametric functions yet!");
  }

  ElementwiseMinimumFunction(const ThisType& other)
    : BaseType(other)
    , some_func_(other.some_func_->copy_as_grid_function())
    , search_quadrature_order_(other.search_quadrature_order_)
    , name_(other.name_)
  {}

  ElementwiseMinimumFunction(ThisType&&) = default;


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
  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<LocalFunction>(*some_func_, search_quadrature_order_);
  }

  std::string name() const override final
  {
    return name_;
  }

private:
  std::unique_ptr<GridFunctionInterface<E, r_, rC_, R>> some_func_;
  const int search_quadrature_order_;
  const std::string name_;
}; // class ElementwiseMinimumFunction


} // namespace Dune::XT::Functions

#endif // DUNE_XT_FUNCTIONS_ELEMENTWISE_MINIMUM_HH
