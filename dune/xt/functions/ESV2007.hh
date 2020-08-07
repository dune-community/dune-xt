// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2018)
//   René Fritze     (2014 - 2016, 2018 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2014 - 2020)

#ifndef DUNEXT_FUNCTIONS_ESV2007_HH
#define DUNEXT_FUNCTIONS_ESV2007_HH

#include <type_traits>

#include <dune/geometry/quadraturerules.hh>

#if HAVE_EIGEN
#  include <dune/xt/common/disable_warnings.hh>
#  include <Eigen/Eigenvalues>
#  include <dune/xt/common/reenable_warnings.hh>
#endif

#include <dune/xt/common/configuration.hh>
#include <dune/xt/la/eigen-solver.hh>

#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/interfaces/function.hh>


namespace Dune::XT::Functions::ESV2007 {


template <size_t d, size_t r, size_t rC = 1, class R = double>
class Testcase1Force : public FunctionInterface<d, r, rC, R>
{
  Testcase1Force()
  {
    static_assert(AlwaysFalse<R>::value, "Not available for these dimensions!");
  }
};


template <class R>
class Testcase1Force<2, 1, 1, R> : public FunctionInterface<2, 1, 1, R>
{
  using ThisType = Testcase1Force;
  using BaseType = FunctionInterface<2, 1, 1, R>;

public:
  using typename BaseType::DerivativeRangeReturnType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeReturnType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".ESV2007.testcase1.force";
  }

  static Common::Configuration defaults()
  {
    Common::Configuration config;
    config["integration_order"] = "3";
    config["name"] = static_id();
    return config;
  } // ... defaults(...)

  Testcase1Force(const size_t ord = defaults().template get<int>("integration_order"),
                 const std::string nm = static_id())
    : order_(static_cast<int>(ord))
    , name_(nm)
  {}

  Testcase1Force(const ThisType& /*other*/) = default;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  std::string name() const final
  {
    return name_;
  }

  int order(const XT::Common::Parameter& /*param*/ = {}) const final
  {
    return order_;
  }

  /**
   * \brief "0.5 * pi * pi * cos(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])"
   */
  RangeReturnType evaluate(const DomainType& xx, const Common::Parameter& /*param*/ = {}) const final
  {
    return (M_PI_2l * M_PIl * cos(M_PI_2l * xx[0]) * cos(M_PI_2l * xx[1]));
  }

  /**
   * \brief ["-0.25 * pi * pi * pi * sin(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])"
   *         "-0.25 * pi * pi * pi * cos(0.5 * pi * x[0]) * sin(0.5 * pi * x[1])"]
   */
  DerivativeRangeReturnType jacobian(const DomainType& xx, const Common::Parameter& /*param*/ = {}) const final
  {
    DerivativeRangeReturnType ret(0.);
    const DomainFieldType pre = -0.25 * M_PIl * M_PIl * M_PIl;
    const DomainFieldType x_arg = M_PI_2l * xx[0];
    const DomainFieldType y_arg = M_PI_2l * xx[1];
    ret[0][0] = pre * sin(x_arg) * cos(y_arg);
    ret[0][1] = pre * cos(x_arg) * sin(y_arg);
    return ret;
  } // ... jacobian(...)

private:
  const int order_;
  const std::string name_;
}; // class Testcase1Force


template <size_t d, size_t r, size_t rC = 1, class R = double>
class Testcase1ExactSolution : public FunctionInterface<d, r, rC, R>
{
  Testcase1ExactSolution()
  {
    static_assert(AlwaysFalse<R>::value, "Not available for these dimensions!");
  }
};


template <class R>
class Testcase1ExactSolution<2, 1, 1, R> : public FunctionInterface<2, 1, 1, R>
{
  using ThisType = Testcase1ExactSolution;
  using BaseType = FunctionInterface<2, 1, 1, R>;

public:
  using typename BaseType::DerivativeRangeReturnType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeReturnType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".ESV2007.testcase1.exactsolution";
  }

  static Common::Configuration defaults()
  {
    Common::Configuration config;
    config["integration_order"] = "3";
    config["name"] = static_id();
    return config;
  } // ... defaults(...)

  Testcase1ExactSolution(const size_t ord = defaults().template get<int>("integration_order"),
                         const std::string nm = static_id())
    : order_(static_cast<int>(ord))
    , name_(nm)
  {}

  Testcase1ExactSolution(const ThisType& /*other*/) = default;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  std::string name() const final
  {
    return name_;
  }

  int order(const XT::Common::Parameter& /*param*/ = {}) const final
  {
    return order_;
  }

  /**
   * \brief "cos(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])"
   */
  RangeReturnType evaluate(const DomainType& xx, const Common::Parameter& /*param*/ = {}) const final
  {
    return (cos(M_PI_2l * xx[0]) * cos(M_PI_2l * xx[1]));
  }

  /**
   * \brief ["-0.5 * pi * sin(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])"
   *         "-0.5 * pi * cos(0.5 * pi * x[0]) * sin(0.5 * pi * x[1])"]
   */
  DerivativeRangeReturnType jacobian(const DomainType& xx, const Common::Parameter& /*param*/ = {}) const final
  {
    DerivativeRangeReturnType ret(0.);
    const DomainFieldType pre = -0.5 * M_PIl;
    const DomainFieldType x_arg = M_PI_2l * xx[0];
    const DomainFieldType y_arg = M_PI_2l * xx[1];
    ret[0][0] = pre * sin(x_arg) * cos(y_arg);
    ret[0][1] = pre * cos(x_arg) * sin(y_arg);
    return ret;
  } // ... jacobian(...)

private:
  const int order_;
  const std::string name_;
}; // class Testcase1ExactSolution


template <class E, class R = double>
class CutoffFunction : public GridFunctionInterface<E, 1, 1, R>
{
  using BaseType = GridFunctionInterface<E, 1, 1, R>;
  using ThisType = CutoffFunction;

public:
  using BaseType::d;
  using typename BaseType::D;
  using DiffusionType = GridFunctionInterface<E, d, d, R>;

private:
  class LocalCutoffFunction : public ElementFunctionInterface<E, 1, 1, R>
  {
    using BaseType = ElementFunctionInterface<E, 1, 1, R>;


  public:
    using typename BaseType::DerivativeRangeReturnType;
    using typename BaseType::DomainFieldType;
    using typename BaseType::DomainType;
    using typename BaseType::ElementType;
    using typename BaseType::RangeFieldType;
    using typename BaseType::RangeReturnType;

    LocalCutoffFunction(const DiffusionType& diffusion, const RangeFieldType poincare_constant)
      : BaseType()
      , value_(0)
      , local_diffusion_(diffusion.local_function())
      , poincare_constant_(poincare_constant)
    {}

  protected:
    void post_bind(const ElementType& ele) final
    {
      post_bind_helper<ElementType, d>::post_bind(ele, value_, local_diffusion_, poincare_constant_);
    }

  public:
    int order(const XT::Common::Parameter& /*param*/ = {}) const final
    {
      return 0;
    }

    RangeReturnType evaluate(const DomainType& DXTC_DEBUG_ONLY(xx), const Common::Parameter& /*param*/ = {}) const final
    {
      this->assert_inside_reference_element(xx);
      return value_;
    }

    DerivativeRangeReturnType jacobian(const DomainType& DXTC_DEBUG_ONLY(xx),
                                       const Common::Parameter& /*param*/ = {}) const final
    {
      this->assert_inside_reference_element(xx);
      return DerivativeRangeReturnType();
    }

  private:
    template <class ElementType, size_t d>
    struct post_bind_helper
    {
      static void post_bind(const ElementType& ele,
                            RangeFieldType& value_,
                            std::unique_ptr<typename DiffusionType::LocalFunctionType>& local_diffusion_,
                            const RangeFieldType& poincare_constant_)
      {
        local_diffusion_->bind(ele);
        value_ = std::numeric_limits<R>::max();
        for (const auto& quadrature_point :
             Dune::QuadratureRules<double, d>::rule(ele.type(), local_diffusion_->order() + 2)) {
          const auto point_in_reference_element = quadrature_point.position();
          Common::FieldMatrix<R, d, d> diffusion_tensor_value(0.);
          diffusion_tensor_value = local_diffusion_->evaluate(point_in_reference_element);
          value_ = std::min(value_,
                            LA::make_eigen_solver(diffusion_tensor_value,
                                                  {{"type", LA::eigen_solver_types(diffusion_tensor_value).at(0)},
                                                   {"assert_positive_eigenvalues", "1e-15"}})
                                .min_eigenvalues(1)
                                .at(0));
        }
        const auto hh = Grid::entity_diameter(ele);
        value_ = (poincare_constant_ * hh * hh) / value_;
      }
    };

    template <class ElementType>
    struct post_bind_helper<ElementType, 1>
    {
      static void post_bind(const ElementType& ele,
                            RangeFieldType& value_,
                            std::unique_ptr<typename DiffusionType::LocalFunctionType>& local_diffusion_,
                            const RangeFieldType& poincare_constant_)
      {
        local_diffusion_->bind(ele);
        value_ = std::numeric_limits<R>::max();
        for (const auto& quadrature_point :
             Dune::QuadratureRules<double, d>::rule(ele.type(), local_diffusion_->order() + 2)) {
          const auto point_in_reference_element = quadrature_point.position();
          value_ = std::min(value_, local_diffusion_->evaluate(point_in_reference_element)[0]);
        }
        const auto hh = Grid::entity_diameter(ele);
        value_ = (poincare_constant_ * hh * hh) / value_;
      }
    };

    RangeFieldType value_;
    std::unique_ptr<typename DiffusionType::LocalFunctionType> local_diffusion_;
    const RangeFieldType poincare_constant_;
  }; // class LocalCutoffFunction

public:
  using ElementType = typename BaseType::ElementType;
  using LocalFunctionType = typename BaseType::LocalFunctionType;
  using RangeFieldType = typename BaseType::RangeFieldType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".ESV2007.cutoff";
  }

  CutoffFunction(const DiffusionType& diffusion,
                 const RangeFieldType poincare_constant = 1.0 / (M_PIl * M_PIl),
                 const std::string nm = static_id())
    : diffusion_(diffusion)
    , poincare_constant_(poincare_constant)
    , name_(nm)
  {}

  CutoffFunction(const ThisType& other) = default;

  ThisType& operator=(const ThisType& other) = delete;

  std::string name() const final
  {
    return name_;
  }

  std::unique_ptr<LocalFunctionType> local_function() const final
  {
    return std::make_unique<LocalCutoffFunction>(diffusion_, poincare_constant_);
  }

private:
  const DiffusionType& diffusion_;
  const RangeFieldType poincare_constant_;
  std::string name_;
}; // class Cutoff


} // namespace Dune::XT::Functions::ESV2007

#endif // DUNEXT_FUNCTIONS_ESV2007_HH
