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

#include <dune/xt/functions/generic/flux-function.hh>

using namespace Dune::XT;

{% for GRIDNAME, GRID, s, r, rC in config['types'] %}

struct GenericFluxFunction_from_{{GRIDNAME}}_and_dim_{{s}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static constexpr size_t d = GridType::dimension;
  static constexpr size_t s = {{s}};
  static constexpr size_t r = {{r}};
  static constexpr size_t rC = {{rC}};

  using GenericType = Functions::GenericFluxFunction<ElementType, s, r, rC>;

  using RangeType = typename GenericType::LocalFunctionType::RangeType;
  using RangeReturnType = typename GenericType::LocalFunctionType::RangeReturnType;
  using DomainType = typename GenericType::LocalFunctionType::DomainType;
  using StateType = typename GenericType::LocalFunctionType::StateType;
  using JacobianRangeType = typename GenericType::LocalFunctionType::JacobianRangeType;
  using JacobianRangeReturnType = typename GenericType::LocalFunctionType::JacobianRangeReturnType;


  GenericFluxFunction_from_{{GRIDNAME}}_and_dim_{{s}}_to_{{r}}_times_{{rC}}()
    : grid_(Grid::make_cube_grid<GridType>(-1., 1., 4))
  {
  }

  const Grid::GridProvider<GridType> grid_;
};


TEST_F(GenericFluxFunction_from_{{GRIDNAME}}_and_dim_{{s}}_to_{{r}}_times_{{rC}}, is_constructible)
{
  size_t element_index = 0;
  const auto leaf_view = grid_.leaf_view();

  GenericType  function(
      /*order=*/0,
      [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
      [&](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return RangeReturnType(element_index); },
      /*parameter=*/{},
      "element_index_",
      [](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return JacobianRangeReturnType(); }
          );

  GenericType function2(
      /*order_func=*/[](const auto& /*param*/) { return 3; },
      [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
      [&](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return RangeReturnType(element_index); },
      /*parameter=*/{},
      "element_index_",
      [](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return JacobianRangeReturnType(); }
          );

  // check that defaults are compiling
  GenericType  function3(/*order=*/0);

  GenericType  function4(/*order_func=*/[](const auto& /*param*/) { return 3; });
}

TEST_F(GenericFluxFunction_from_{{GRIDNAME}}_and_dim_{{s}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  size_t element_index = 0;
  const auto leaf_view = grid_.leaf_view();

  GenericType  function(
      /*order=*/0,
      [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
      [&](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return RangeReturnType(element_index); },
      /*parameter=*/{},
      "element_index_",
      [](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return JacobianRangeReturnType(); }
          );
  auto local_f = function.local_function();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
  }
}

TEST_F(GenericFluxFunction_from_{{GRIDNAME}}_and_dim_{{s}}_to_{{r}}_times_{{rC}}, local_order)
{
  size_t element_index = 0;
  const auto leaf_view = grid_.leaf_view();

  for (auto vv : {1, 3, 5}) {
    const int expected_order = vv;
    GenericType  function(
        /*order=*/vv,
        [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
        [&](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return RangeReturnType(element_index); },
        /*parameter=*/{},
        "element_index_",
        [](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return JacobianRangeReturnType(); });

    auto local_f = function.local_function();
    for (auto&& element : Dune::elements(leaf_view)) {
      local_f->bind(element);
      const auto actual_order = local_f->order();
      EXPECT_EQ(expected_order, actual_order);
    }
  }
}

TEST_F(GenericFluxFunction_from_{{GRIDNAME}}_and_dim_{{s}}_to_{{r}}_times_{{rC}}, local_evaluate)
{
  size_t element_index = 0;
  const auto leaf_view = grid_.leaf_view();

  GenericType  function(
      /*order=*/0,
      [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
      [&](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return RangeReturnType(element_index); },
      /*parameter=*/{},
      "element_index_",
      [](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return JacobianRangeReturnType(); }
  );
  auto local_f = function.local_function();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    StateType u;
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      Common::ParameterType param("power", 1);
      RangeReturnType expected_value(leaf_view.indexSet().index(element));
      const auto actual_value = local_f->evaluate(local_x, u);
      EXPECT_EQ(expected_value, actual_value);
    }
  }
}

TEST_F(GenericFluxFunction_from_{{GRIDNAME}}_and_dim_{{s}}_to_{{r}}_times_{{rC}}, local_jacobian)
{
  size_t element_index = 0;
  const auto leaf_view = grid_.leaf_view();

  GenericType  function(
      /*order=*/0,
      [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
      [&](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return RangeReturnType(element_index); },
      /*parameter=*/{},
      "element_index_",
      [](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return JacobianRangeReturnType(); }
          );
  auto local_f = function.local_function();
  JacobianRangeReturnType expected_jacobian;
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

{% endfor  %}
