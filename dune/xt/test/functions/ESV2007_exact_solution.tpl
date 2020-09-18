// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   René Fritze    (2019)
//   Tim Keil       (2018)
//   Tobias Leibner (2018)

#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/xt/test/main.hxx>

#include <dune/xt/grid/grids.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/ESV2007.hh>
#include <dune/xt/functions/grid-function.hh>

using namespace Dune::XT;

{% for GRIDNAME, GRID, r, rC in config['types'] %}


struct ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static constexpr size_t d = GridType::dimension;
  static constexpr size_t r = {{r}};
  static constexpr size_t rC = {{rC}};

  using FunctionType = Functions::ESV2007::Testcase1ExactSolution<d, r, rC>;

  using RangeReturnType = typename FunctionType::RangeReturnType;
  using DomainType = typename FunctionType::DomainType;
  using DerivativeRangeReturnType = typename FunctionType::DerivativeRangeReturnType;

  ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}()
    : grid_(Dune::XT::Grid::make_cube_grid<GridType>(-1., 1., 4))
  {
    grid_.visualize("grid");
  }

  const Grid::GridProvider<GridType> grid_;
};


TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_constructible)
{
  FunctionType function(3);
}

TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, has_default_config)
{
  auto cfg = FunctionType::defaults();
  EXPECT_EQ(cfg.get<std::string>("name"), FunctionType::static_id());
}


TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_creatable)
{
  // auto default_function = FunctionType::create();
  // EXPECT_EQ(default_function->order(), 3);
}


TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_visualizable)
{
  const auto leaf_view = grid_.leaf_view();
  FunctionType function(3);
  function.visualize(leaf_view, "test__ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");
}


TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_order)
{
  const int expected_order = 3;
  FunctionType default_function(3);
  const auto actual_order = default_function.order();
  EXPECT_EQ(expected_order, actual_order);
}


TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_evaluate)
{
  FunctionType function(3);
  for (auto point : {0.25, 0.5, 0.75}) {
    const DomainType xx(point);
    const RangeReturnType expected_value(cos(M_PI_2 * point) * cos(M_PI_2 * point));
    const auto actual_value = function.evaluate(xx);
    EXPECT_EQ(expected_value, actual_value);
  }
}

TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_jacobian)
{
  FunctionType function(3);
  for (auto point : {0.25, 0.5, 0.75}) {
    const DomainType xx(point);
    DerivativeRangeReturnType ret(0.);
    const double pre = -0.5 * M_PI;
    const double x_arg = M_PI_2 * xx[0];
    const double y_arg = M_PI_2 * xx[1];
    ret[0][0] = pre * sin(x_arg) * cos(y_arg);
    ret[0][1] = pre * cos(x_arg) * sin(y_arg);
    const auto actual_value = function.jacobian(xx);
    EXPECT_EQ(ret, actual_value);
  }
}

TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  FunctionType default_function(3);
  auto localizable_function = Functions::make_grid_function<ElementType>(default_function);
  auto local_f = localizable_function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
  }
}


TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_order)
{
  const int expected_order = 3;
  FunctionType default_function(3);
  auto localizable_function = Functions::make_grid_function<ElementType>(default_function);
  auto local_f = localizable_function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    const auto actual_order = local_f->order();
    EXPECT_EQ(expected_order, actual_order);
  }
}


TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_evaluate)
{
  FunctionType function(3);
  auto localizable_function = Functions::make_grid_function<ElementType>(function);
  auto local_f = localizable_function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto value = local_f->evaluate(local_x);
      const auto point = element.geometry().global(local_x);
      const RangeReturnType expected_value(cos(M_PI_2 * point[0]) * cos(M_PI_2 * point[1]));
      EXPECT_EQ(expected_value, value);
    }
  }
}


TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_jacobian)
{
  FunctionType function(3);
  auto localizable_function = Functions::make_grid_function<ElementType>(function);
  auto local_f = localizable_function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto value = local_f->jacobian(local_x);
      const auto point = element.geometry().global(local_x);
      DerivativeRangeReturnType ret(0.);
      const double pre = -0.5 * M_PI;
      const double x_arg = M_PI_2 * point[0];
      const double y_arg = M_PI_2 * point[1];
      ret[0][0] = pre * sin(x_arg) * cos(y_arg);
      ret[0][1] = pre * cos(x_arg) * sin(y_arg);
      EXPECT_EQ(ret, value);
    }
  }
}


{% endfor  %}
