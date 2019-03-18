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

#include <dune/xt/functions/generic/flux-function.hh>

using namespace Dune::XT;

{% for GRIDNAME, GRID, s, r, rC in config['types'] %}

struct GenericFluxFunction_from_{{GRIDNAME}}_and_dim_{{s}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static const constexpr size_t d = GridType::dimension;
  static const size_t s = {{s}};
  static const size_t r = {{r}};
  static const size_t rC = {{rC}};

  using GenericType = Functions::GenericFluxFunction<ElementType, s, r, rC>;

  using RangeType = typename GenericType::LocalFunctionType::RangeType;
  using RangeReturnType = typename GenericType::LocalFunctionType::RangeReturnType;
  using DomainType = typename GenericType::LocalFunctionType::DomainType;
  using StateType = typename GenericType::LocalFunctionType::StateType;
  using PartialURangeType = typename GenericType::LocalFunctionType::PartialURangeType;
  using PartialURangeReturnType = typename GenericType::LocalFunctionType::PartialURangeReturnType;


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
      [](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return PartialURangeReturnType(); }
          );

  GenericType function2(
      /*order_func=*/[](const auto& /*param*/) { return 3; },
      [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
      [&](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return RangeReturnType(element_index); },
      /*parameter=*/{},
      "element_index_",
      [](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return PartialURangeReturnType(); }
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
      [](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return PartialURangeReturnType(); }
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
        [](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return PartialURangeReturnType(); });

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
      [](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return PartialURangeReturnType(); }
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

TEST_F(GenericFluxFunction_from_{{GRIDNAME}}_and_dim_{{s}}_to_{{r}}_times_{{rC}}, local_partial_u)
{
  size_t element_index = 0;
  const auto leaf_view = grid_.leaf_view();

  GenericType  function(
      /*order=*/0,
      [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
      [&](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return RangeReturnType(element_index); },
      /*parameter=*/{},
      "element_index_",
      [](const auto& /*xx*/, const auto& /*u*/, const auto& /*param*/) { return PartialURangeReturnType(); }
          );
  auto local_f = function.local_function();
  PartialURangeReturnType expected_partial_u;
  StateType u;
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto actual_partial_u = local_f->partial_u(local_x, u);
      EXPECT_EQ(expected_partial_u, actual_partial_u);
    }
  }
}

{% endfor  %}
