// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   René Fritze    (2019)
//   Tim Keil       (2018)
//   Tobias Leibner (2018 - 2019)

#include <dune/xt/test/main.hxx>

#include <dune/xt/grid/grids.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/grid-function.hh>
#include <dune/xt/functions/visualization.hh>

using namespace Dune::XT;
using namespace Dune::XT::Functions;

{% for GRIDNAME, GRID, r, rC in config['types'] %}

struct ConstantFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static constexpr size_t d = GridType::dimension;
  static constexpr size_t r = {{r}};
  static constexpr size_t rC = {{rC}};

  using FunctionType = ConstantFunction<d, r, rC>;

  using RangeReturnType = typename FunctionType::RangeReturnType;
  using DomainType = typename FunctionType::DomainType;
  using DerivativeRangeReturnType = typename FunctionType::DerivativeRangeReturnType;

  ConstantFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}()
    : grid_(Grid::make_cube_grid<GridType>(-1., 1., 4))
  {
  }

  const Grid::GridProvider<GridType> grid_;
};


TEST_F(ConstantFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_constructible)
{
  // constructor taking FieldType
  FunctionType function1(0.);
  FunctionType function2(0., "name");
  // constructor taking RangeType
  FunctionType function3(RangeReturnType(0.));
  FunctionType function4(RangeReturnType(0.), "name");
}

TEST_F(ConstantFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, has_default_config)
{
  auto cfg = FunctionType::defaults();
  EXPECT_EQ(cfg.get<std::string>("name"), FunctionType::static_id());
}

TEST_F(ConstantFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_visualizable)
{
  const auto leaf_view = grid_.leaf_view();
  FunctionType function(1.);
  visualize(function, leaf_view, "test__ConstantFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");
}

TEST_F(ConstantFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_order)
{
  const int expected_order = 0;
  for (auto vv : {-10., 3., 17., 41.}) {
    const RangeReturnType value(vv);
    FunctionType function(value);
    const auto actual_order = function.order();
    EXPECT_EQ(expected_order, actual_order);
  }
}

TEST_F(ConstantFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_evaluate)
{
  for (auto value : {-10., 3., 17., 41.}) {
    const RangeReturnType expected_value(value);
    FunctionType function(expected_value);
    for (auto point : {-100., -10., 0., 10., 100.}) {
      const DomainType xx(point);
      const auto actual_value = function.evaluate(xx);
      EXPECT_EQ(expected_value, actual_value);
    }
  }
}

TEST_F(ConstantFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_jacobian)
{
  for (auto vv : {-10., 3., 17., 41.}) {
    const RangeReturnType value(vv);
    DerivativeRangeReturnType expected_jacobian;
    //expected_jacobian *= 0;
    FunctionType function(value);
    for (auto point : {-100., -10., 0., 10., 100.}) {
      const DomainType xx(point);
      const auto actual_jacobian = function.jacobian(xx);
      EXPECT_EQ(expected_jacobian, actual_jacobian);
    }
  }
}

TEST_F(ConstantFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  FunctionType function(1.);
  auto localizable_function = make_grid_function<ElementType>(function);
  auto local_f = localizable_function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
  }
}

TEST_F(ConstantFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_order)
{
  const int expected_order = 0;
  for (auto vv : {-10., 3., 17., 41.}) {
    const RangeReturnType value(vv);
    FunctionType function(value);
    auto localizable_function = make_grid_function<ElementType>(function);
    auto local_f = localizable_function.local_function();
    const auto leaf_view = grid_.leaf_view();
    for (auto&& element : Dune::elements(leaf_view)) {
      local_f->bind(element);
      const auto actual_order = local_f->order();
      EXPECT_EQ(expected_order, actual_order);
    }
  }
}

TEST_F(ConstantFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_evaluate)
{
  for (auto value : {-10., 3., 17., 41.}) {
    const RangeReturnType expected_value(value);
    FunctionType function(expected_value);
    auto localizable_function = make_grid_function<ElementType>(function);
    auto local_f = localizable_function.local_function();
    const auto leaf_view = grid_.leaf_view();
    for (auto&& element : Dune::elements(leaf_view)) {
      local_f->bind(element);
      for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
        const auto local_x = quadrature_point.position();
        const auto actual_value = local_f->evaluate(local_x);
        EXPECT_EQ(expected_value, actual_value);
      }
    }
  }
}

TEST_F(ConstantFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_jacobian)
{
  for (auto vv : {-10., 3., 17., 41.}) {
    RangeReturnType value(vv);
    DerivativeRangeReturnType expected_jacobian;
    //expected_jacobian *= 0;
    FunctionType function(value);
    auto localizable_function = make_grid_function<ElementType>(function);
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
}


struct ConstantGridFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static constexpr size_t d = GridType::dimension;
  static constexpr size_t r = {{r}};
  static constexpr size_t rC = {{rC}};

  using FunctionType = ConstantGridFunction<ElementType, r, rC>;

  using RangeReturnType = typename FunctionType::LocalFunctionType::RangeReturnType;
  using DomainType = typename FunctionType::LocalFunctionType::DomainType;
  using DerivativeRangeReturnType = typename FunctionType::LocalFunctionType::DerivativeRangeReturnType;

  ConstantGridFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}()
    : grid_(Grid::make_cube_grid<GridType>(-1., 1., 4))
  {
  }

  const Grid::GridProvider<GridType> grid_;
};


TEST_F(ConstantGridFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_constructible)
{
  // passing FieldType
  FunctionType function1(0.);
  FunctionType function2(0., "name");
  // passing RangeType
  FunctionType function3(RangeReturnType(0.));
  FunctionType function4(RangeReturnType(0.), "name");
}

TEST_F(ConstantGridFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_visualizable)
{
  const auto leaf_view = grid_.leaf_view();
  FunctionType function(1.);
  visualize(function, leaf_view, "test__ConstantFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");
}

TEST_F(ConstantGridFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  FunctionType function(1.);
  auto local_f = function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
  }
}

TEST_F(ConstantGridFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_order)
{
  const int expected_order = 0;
  for (auto vv : {-10., 3., 17., 41.}) {
    const RangeReturnType value(vv);
    FunctionType function(value);
    auto local_f = function.local_function();
    const auto leaf_view = grid_.leaf_view();
    for (auto&& element : Dune::elements(leaf_view)) {
      local_f->bind(element);
      const auto actual_order = local_f->order();
      EXPECT_EQ(expected_order, actual_order);
    }
  }
}

TEST_F(ConstantGridFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_evaluate)
{
  for (auto value : {-10., 3., 17., 41.}) {
    const RangeReturnType expected_value(value);
    FunctionType function(expected_value);
    auto local_f = function.local_function();
    const auto leaf_view = grid_.leaf_view();
    for (auto&& element : Dune::elements(leaf_view)) {
      local_f->bind(element);
      for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
        const auto local_x = quadrature_point.position();
        const auto actual_value = local_f->evaluate(local_x);
        EXPECT_EQ(expected_value, actual_value);
      }
    }
  }
}

TEST_F(ConstantGridFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_jacobian)
{
  for (auto vv : {-10., 3., 17., 41.}) {
    RangeReturnType value(vv);
    DerivativeRangeReturnType expected_jacobian;
    //expected_jacobian *= 0;
    FunctionType function(value);
    auto local_f = function.local_function();
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
}


struct ConstantFluxFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static constexpr size_t d = GridType::dimension;
  static constexpr size_t s = 2;
  static constexpr size_t r = {{r}};
  static constexpr size_t rC = {{rC}};

  using FunctionType = ConstantFluxFunction<ElementType, s, r, rC>;

  using RangeReturnType = typename FunctionType::LocalFunctionType::RangeReturnType;
  using DomainType = typename FunctionType::LocalFunctionType::DomainType;
  using StateType = typename FunctionType::LocalFunctionType::StateType;
  using JacobianRangeReturnType = typename FunctionType::LocalFunctionType::JacobianRangeReturnType;

  ConstantFluxFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}()
    : grid_(Grid::make_cube_grid<GridType>(-1., 1., 4))
  {
  }

  const Grid::GridProvider<GridType> grid_;
};


TEST_F(ConstantFluxFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_constructible)
{
  // passing FieldType
  FunctionType function1(0.);
  FunctionType function2(0., "name");
  // passing RangeType
  FunctionType function3(RangeReturnType(0.));
  FunctionType function4(RangeReturnType(0.), "name");
}

TEST_F(ConstantFluxFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  FunctionType function(1.);
  auto local_f = function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
  }
}

TEST_F(ConstantFluxFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_order)
{
  const int expected_order = 0;
  for (auto vv : {-10., 3., 17., 41.}) {
    const RangeReturnType value(vv);
    FunctionType function(value);
    auto local_f = function.local_function();
    const auto leaf_view = grid_.leaf_view();
    for (auto&& element : Dune::elements(leaf_view)) {
      local_f->bind(element);
      const auto actual_order = local_f->order();
      EXPECT_EQ(expected_order, actual_order);
    }
  }
}

TEST_F(ConstantFluxFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_evaluate)
{
  for (auto value : {-10., 3., 17., 41.}) {
    const RangeReturnType expected_value(value);
    FunctionType function(expected_value);
    auto local_f = function.local_function();
    const auto leaf_view = grid_.leaf_view();
    StateType u;
    for (auto&& element : Dune::elements(leaf_view)) {
      local_f->bind(element);
      for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
        const auto local_x = quadrature_point.position();
        const auto actual_value = local_f->evaluate(local_x, u);
        EXPECT_EQ(expected_value, actual_value);
      }
    }
  }
}

TEST_F(ConstantFluxFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_jacobian)
{
  for (auto vv : {-10., 3., 17., 41.}) {
    RangeReturnType value(vv);
    JacobianRangeReturnType expected_jacobian;
    //expected_jacobian *= 0;
    FunctionType function(value);
    auto local_f = function.local_function();
    const auto leaf_view = grid_.leaf_view();
    StateType u;
    for (auto&& element : Dune::elements(leaf_view)) {
      local_f->bind(element);
      for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
        const auto local_x = quadrature_point.position();
        const auto actual_jacobian = local_f->jacobian(local_x, u);
        EXPECT_EQ(expected_jacobian, actual_jacobian);
      }
    }
  }
}




{% endfor  %}
