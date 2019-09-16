// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2019 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2015 - 2016, 2018)
//   Tobias Leibner  (2015 - 2016)

#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/grid/grids.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/ESV2007.hh>
#include <dune/xt/functions/expression.hh>

using namespace Dune::XT;

{% for GRIDNAME, GRID in config['types'] %}


struct ESV2007CutoffFunction_from_{{GRIDNAME}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static const constexpr size_t d = GridType::dimension;
  static const size_t r = d;
  static const size_t rC = d;

  using DiffusionType = Functions::ExpressionFunction<d, r, rC>;
  using FunctionType = Functions::ESV2007::CutoffFunction<ElementType>;

  using ScalarRangeExpressionType = typename Dune::XT::Common::FieldVector<std::string, 1>;

  using RangeReturnType = typename FunctionType::LocalFunctionType::RangeReturnType;
  using DomainType = typename FunctionType::LocalFunctionType::DomainType;
  using DerivativeRangeReturnType = typename FunctionType::LocalFunctionType::DerivativeRangeReturnType;

  ESV2007CutoffFunction_from_{{GRIDNAME}}()
    : grid_(Dune::XT::Grid::make_cube_grid<GridType>(1, 3., 4))
  {
    grid_.visualize("grid");
  }

  const Grid::GridProvider<GridType> grid_;
};

TEST_F(ESV2007CutoffFunction_from_{{GRIDNAME}}, is_constructible)
{
  ScalarRangeExpressionType expr(std::string("x[0]"));
  DiffusionType diffusion_function("x", expr, 1);
  FunctionType function(diffusion_function.template as_grid_function<ElementType>());
}

TEST_F(ESV2007CutoffFunction_from_{{GRIDNAME}}, is_visualizable)
{
  const auto leaf_view = grid_.leaf_view();
  ScalarRangeExpressionType expr(std::string("x[0]"));
  DiffusionType diffusion_function("x", expr, 1);
  FunctionType function(diffusion_function.template as_grid_function<ElementType>());
  function.visualize(leaf_view, "test__ESV2007CutoffFunction_from_{{GRIDNAME}}__is_visualizable");
}

TEST_F(ESV2007CutoffFunction_from_{{GRIDNAME}}, is_bindable)
{
  ScalarRangeExpressionType expr(std::string("x[0]"));
  DiffusionType diffusion_function("x", expr, 1);
  FunctionType function(diffusion_function.template as_grid_function<ElementType>());
  auto local_f = function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
  }
}


TEST_F(ESV2007CutoffFunction_from_{{GRIDNAME}}, local_order)
{
  const int expected_order = 0;
  ScalarRangeExpressionType expr(std::string("x[0]"));
  DiffusionType diffusion_function("x", expr, 1);
  FunctionType function(diffusion_function.template as_grid_function<ElementType>());
  auto local_f = function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    const auto actual_order = local_f->order();
    EXPECT_EQ(expected_order, actual_order);
  }
}


TEST_F(ESV2007CutoffFunction_from_{{GRIDNAME}}, local_evaluate)
{
  ScalarRangeExpressionType expr(std::string("x[0]"));
  DiffusionType diffusion_function("x", expr, 1);
  FunctionType function(diffusion_function.template as_grid_function<ElementType>());
  auto local_f = function.local_function();
  auto local_diffusion = diffusion_function.template as_grid_function<ElementType>().local_function();
  double poincare_constant_ = 1.0 / (M_PIl * M_PIl);
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    local_diffusion->bind(element);
    double value_ = std::numeric_limits<double>::max();
    double value = std::numeric_limits<double>::max();
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), local_diffusion->order() + 2)) {
      const auto local_x = quadrature_point.position();
      value = local_f->evaluate(local_x);
      Common::FieldMatrix<double, d, d> diffusion_tensor_value(0.);
      diffusion_tensor_value[0][0] = local_diffusion->evaluate(local_x)[0];
      value_ = std::min(value_,
                       LA::make_eigen_solver(diffusion_tensor_value,
       {{'{{"type", LA::eigen_solver_types(diffusion_tensor_value).at(0)}, {"assert_positive_eigenvalues", "1e-15"}}'}})
                           .min_eigenvalues(1)
                           .at(0));
      }
    const auto hh = Grid::entity_diameter(element);
    value_ = (poincare_constant_ * hh * hh) / value_;
    EXPECT_EQ(value_, value);
  }
}


TEST_F(ESV2007CutoffFunction_from_{{GRIDNAME}}, local_jacobian)
{
  ScalarRangeExpressionType expr(std::string("x[0]"));
  DiffusionType diffusion_function("x", expr, 1);
  FunctionType function(diffusion_function.template as_grid_function<ElementType>());
  auto local_f = function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    DerivativeRangeReturnType expected_jacobian = DerivativeRangeReturnType();
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto value = local_f->jacobian(local_x);
      EXPECT_EQ(expected_jacobian, value);
    }
  }
}


{% endfor  %}
