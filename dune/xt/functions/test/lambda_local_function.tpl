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

#include <dune/xt/grid/grids.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/lambda/local-function.hh>

using namespace Dune::XT;

{% for GRIDNAME, GRID, r, rC in config['types'] %}

struct LocalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static const constexpr size_t d = GridType::dimension;
  static const size_t r = {{r}};
  static const size_t rC = {{rC}};

  using LocalLambdaType = Functions::LocalLambdaFunction<ElementType, r, rC>;

  using RangeType = typename LocalLambdaType::LocalFunctionType::RangeType;
  using DomainType = typename LocalLambdaType::LocalFunctionType::DomainType;
  using DerivativeRangeType = typename LocalLambdaType::LocalFunctionType::DerivativeRangeType;

  LocalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}()
    : grid_(Grid::make_cube_grid<GridType>(-1., 1., 4))
  {
  }

  const Grid::GridProvider<GridType> grid_;
};


TEST_F(LocalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_constructible)
{
  size_t element_index = 0;
  const auto leaf_view = grid_.leaf_view();

  LocalLambdaType  function(
      /*order=*/0,
      [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
      [&](const auto& /*xx*/, const auto& /*param*/) { return RangeType(element_index); },
      /*parameter=*/{},
      "element_index_",
      [](const auto& /*element*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeType(); },
      [](const auto& /*element*/, const auto& /*alpha*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeType(); });
}


TEST_F(LocalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_localizable)
{
  size_t element_index = 0;
  const auto leaf_view = grid_.leaf_view();

  LocalLambdaType  function(
      /*order=*/0,
      [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
      [&](const auto& /*xx*/, const auto& /*param*/) { return RangeType(element_index); },
      /*parameter=*/{},
      "element_index_",
      [](const auto& /*element*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeType(); },
      [](const auto& /*element*/, const auto& /*alpha*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeType(); });
  for (auto&& element : Dune::elements(leaf_view)) {
    const auto local_f = function.local_function(element);
  }
}

TEST_F(LocalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  size_t element_index = 0;
  const auto leaf_view = grid_.leaf_view();

  LocalLambdaType  function(
      /*order=*/0,
      [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
      [&](const auto& /*xx*/, const auto& /*param*/) { return RangeType(element_index); },
      /*parameter=*/{},
      "element_index_",
      [](const auto& /*element*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeType(); },
      [](const auto& /*element*/, const auto& /*alpha*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeType(); });
  auto local_f = function.local_function();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
  }
}

TEST_F(LocalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_order)
{
  size_t element_index = 0;
  const auto leaf_view = grid_.leaf_view();

  for (auto vv : {1, 3, 5}) {
    const int expected_order = vv;
    LocalLambdaType  function(
        /*order=*/vv,
        [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
        [&](const auto& /*xx*/, const auto& /*param*/) { return RangeType(element_index); },
        /*parameter=*/{},
        "element_index_",
        [](const auto& /*element*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeType(); },
        [](const auto& /*element*/, const auto& /*alpha*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeType(); });

    auto local_f = function.local_function();
    for (auto&& element : Dune::elements(leaf_view)) {
      local_f->bind(element);
      const auto actual_order = local_f->order();
      EXPECT_EQ(expected_order, actual_order);
    }
  }
}

TEST_F(LocalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_evaluate)
{
  size_t element_index = 0;
  const auto leaf_view = grid_.leaf_view();

  LocalLambdaType  function(
      /*order=*/0,
      [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
      [&](const auto& /*xx*/, const auto& /*param*/) { return RangeType(element_index); },
      /*parameter=*/{},
      "element_index_",
      [](const auto& /*element*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeType(); },
      [](const auto& /*element*/, const auto& /*alpha*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeType(); });
  auto local_f = function.local_function();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      Common::ParameterType param("power", 1);
      RangeType expected_value(leaf_view.indexSet().index(element));
      const auto actual_value = local_f->evaluate(local_x);
      EXPECT_EQ(expected_value, actual_value);
    }
  }
}

TEST_F(LocalLambdaFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_jacobian)
{
  size_t element_index = 0;
  const auto leaf_view = grid_.leaf_view();

  LocalLambdaType  function(
      /*order=*/0,
      [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
      [&](const auto& /*xx*/, const auto& /*param*/) { return RangeType(element_index); },
      /*parameter=*/{},
      "element_index_",
      [](const auto& /*element*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeType(); },
      [](const auto& /*element*/, const auto& /*alpha*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeType(); });
  auto local_f = function.local_function();
  DerivativeRangeType expected_jacobian;
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
