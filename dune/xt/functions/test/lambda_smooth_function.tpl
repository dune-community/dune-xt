// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014 - 2016, 2018)
//   Tobias Leibner  (2014 - 2017)

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/geometry/quadraturerules.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/lambda/smooth-function.hh>

using namespace Dune::XT;

{% for GRIDNAME, GRID, r, rC in config['types'] %}

struct GlobalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static const constexpr size_t d = GridType::dimension;
  static const size_t r = {{r}};
  static const size_t rC = {{rC}};

  using LambdaType = Functions::SmoothLambdaFunction<d, r, rC>;

  using RangeType = typename LambdaType::RangeType;
  using DomainType = typename LambdaType::DomainType;
  using DerivativeRangeType = typename LambdaType::DerivativeRangeType;

  GlobalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}()
    : grid_(Grid::make_cube_grid<GridType>(-1., 1., 4))
  {
  }

  const Grid::GridProvider<GridType> grid_;
};


TEST_F(GlobalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_constructible)
{
  LambdaType lambda(1, [](const auto& /*xx*/, const auto& /*param*/) { RangeType ret(1.); return ret;});
}

TEST_F(GlobalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_visualizable)
{
  const auto leaf_view = grid_.leaf_view();
  LambdaType function(1, [](const auto& /*xx*/, const auto& /*param*/) { RangeType ret(1.); return ret;});
  function.visualize(leaf_view, "test__GlobalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");
}

TEST_F(GlobalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_order)
{
  for (auto vv : {1, 2, 5}) {
    LambdaType function(vv);
    const auto actual_order = function.order();
    EXPECT_EQ(vv, actual_order);
  }
}

TEST_F(GlobalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_evaluate)
{
  LambdaType function(4,
                      [](const auto& xx, const auto& param) { RangeType ret(std::pow(xx[0], param.get("power").at(0))); return ret;},
                      "x_power_p",
                      Common::ParameterType("power", 1));
  for (auto point : {-100., -10., 0., 10., 100.}) {
    const DomainType xx(point);
    const RangeType expected_value(std::pow(xx[0], 2.));
    const auto actual_value = function.evaluate(xx, {"power", 2.});
    EXPECT_EQ(expected_value, actual_value);
  }
}

TEST_F(GlobalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_jacobian)
{
  DerivativeRangeType expected_jacobian;
  LambdaType function(1,
                      [](const auto& /*xx*/, const auto& /*param*/) { RangeType ret(1.); return ret;},
                      "jacobian.function",
                      {},
                      [](const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeType();});
  for (auto point : {-100., -10., 0., 10., 100.}) {
    const DomainType xx(point);
    const auto actual_jacobian = function.jacobian(xx);
    EXPECT_EQ(expected_jacobian, actual_jacobian);
  }
}

TEST_F(GlobalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_localizable)
{
  LambdaType function(1.);
  const auto& localizable_function = function.template as_localizable<ElementType>();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    const auto local_f = localizable_function.local_function(element);
  }
}

TEST_F(GlobalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  LambdaType function(1.);
  const auto& localizable_function = function.template as_localizable<ElementType>();
  auto local_f = localizable_function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
  }
}

TEST_F(GlobalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_order)
{
  for (auto vv : {1, 3, 5}) {
    const int expected_order = vv;
    LambdaType function(vv);
    const auto& localizable_function = function.template as_localizable<ElementType>();
    auto local_f = localizable_function.local_function();
    const auto leaf_view = grid_.leaf_view();
    for (auto&& element : Dune::elements(leaf_view)) {
      local_f->bind(element);
      const auto actual_order = local_f->order();
      EXPECT_EQ(expected_order, actual_order);
    }
  }
}

TEST_F(GlobalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_evaluate)
{
  LambdaType function(4,
                      [](const auto& xx, const auto& param) { RangeType ret(std::pow(xx[0], param.get("power").at(0))); return ret;},
                      "x_power_p",
                      Common::ParameterType("power", 1));
  const auto& localizable_function = function.template as_localizable<ElementType>();
  auto local_f = localizable_function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    const auto geometry = element.geometry();
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto global_x = geometry.global(local_x);
      const auto actual_value = local_f->evaluate(local_x, {"power", 2});
      const RangeType expected_value(std::pow(global_x[0], 2.));
      EXPECT_EQ(expected_value, actual_value);
    }
  }
}

TEST_F(GlobalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_jacobian)
{
  DerivativeRangeType expected_jacobian;
  LambdaType function(1,
                      [](const auto& /*xx*/, const auto& /*param*/) { RangeType ret(1.); return ret;},
                      "jacobian.function",
                      {},
                      [](const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeType();});
  const auto& localizable_function = function.template as_localizable<ElementType>();
  auto local_f = localizable_function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto actual_jacobian = local_f->jacobian(local_x);
      EXPECT_EQ(expected_jacobian, actual_jacobian);
    }
  }
}

{% endfor  %}
