// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
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

using namespace Dune::XT;

{% for GRIDNAME, GRID, r, rC in config['types'] %}


struct ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static const constexpr size_t d = GridType::dimension;
  static const size_t r = {{r}};
  static const size_t rC = {{rC}};

  using FunctionType = Functions::ESV2007::Testcase1ExactSolution<d, r, rC>;

  using RangeType = typename FunctionType::RangeType;
  using DomainType = typename FunctionType::DomainType;
  using DerivativeRangeType = typename FunctionType::DerivativeRangeType;

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
  auto cfg = FunctionType::default_config();
  EXPECT_EQ(cfg.get<std::string>("type"), FunctionType::static_id());
}


TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_creatable)
{
  auto default_function = FunctionType::create();
  EXPECT_EQ(default_function->order(), 3);
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
  auto default_function = FunctionType::create();
  const auto actual_order = default_function->order();
  EXPECT_EQ(expected_order, actual_order);
}


TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_evaluate)
{
  FunctionType function(3);
  for (auto point : {0.25, 0.5, 0.75}) {
    const DomainType xx(point);
    const RangeType expected_value(cos(M_PI_2l * point) * cos(M_PI_2l * point));
    const auto actual_value = function.evaluate(xx);
    EXPECT_EQ(expected_value, actual_value);
  }
}

TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_jacobian)
{
  FunctionType function(3);
  for (auto point : {0.25, 0.5, 0.75}) {
    const DomainType xx(point);
    DerivativeRangeType ret(0.);
    const double pre = -0.5 * M_PIl;
    const double x_arg = M_PI_2l * xx[0];
    const double y_arg = M_PI_2l * xx[1];
    ret[0][0] = pre * sin(x_arg) * cos(y_arg);
    ret[0][1] = pre * cos(x_arg) * sin(y_arg);
    const auto actual_value = function.jacobian(xx);
    EXPECT_EQ(ret, actual_value);
  }
}



TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_localizable)
{
  auto default_function = FunctionType::create();
  const auto& localizable_function = default_function->template as_grid_function<ElementType>();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    const auto local_f = localizable_function.local_function(element);
  }
}

TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  auto default_function = FunctionType::create();
  const auto& localizable_function = default_function->template as_grid_function<ElementType>();
  auto local_f = localizable_function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
  }
}


TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_order)
{
  const int expected_order = 3;
  auto default_function = FunctionType::create();
  const auto& localizable_function = default_function->template as_grid_function<ElementType>();
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
  const auto& localizable_function = function.template as_grid_function<ElementType>();
  auto local_f = localizable_function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto value = local_f->evaluate(local_x);
      const auto point = element.geometry().global(local_x);
      const RangeType expected_value(cos(M_PI_2l * point[0]) * cos(M_PI_2l * point[1]));
      EXPECT_EQ(expected_value, value);
    }
  }
}


TEST_F(ESV2007ExactSolutionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_jacobian)
{
  FunctionType function(3);
  const auto& localizable_function = function.template as_grid_function<ElementType>();
  auto local_f = localizable_function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto value = local_f->jacobian(local_x);
      const auto point = element.geometry().global(local_x);
      DerivativeRangeType ret(0.);
      const double pre = -0.5 * M_PIl;
      const double x_arg = M_PI_2l * point[0];
      const double y_arg = M_PI_2l * point[1];
      ret[0][0] = pre * sin(x_arg) * cos(y_arg);
      ret[0][1] = pre * cos(x_arg) * sin(y_arg);
      EXPECT_EQ(ret, value);
    }
  }
}


{% endfor  %}
