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

#include <dune/xt/functions/generic/grid-function.hh>

using namespace Dune::XT;

{% for GRIDNAME, GRID, r, rC in config['types'] %}

struct GenericFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static const constexpr size_t d = GridType::dimension;
  static const size_t r = {{r}};
  static const size_t rC = {{rC}};

  using GenericType = Functions::GenericGridFunction<ElementType, r, rC>;

  using RangeType = typename GenericType::LocalFunctionType::RangeType;
  using RangeReturnType = typename GenericType::LocalFunctionType::RangeReturnType;
  using DomainType = typename GenericType::LocalFunctionType::DomainType;
  using DerivativeRangeType = typename GenericType::LocalFunctionType::DerivativeRangeType;
  using DerivativeRangeReturnType = typename GenericType::LocalFunctionType::DerivativeRangeReturnType;


  GenericFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}()
    : grid_(Grid::make_cube_grid<GridType>(-1., 1., 4))
  {
  }

  template<size_t range_cols = rC, bool anything = true>
  struct JacobianHelper
  {
    template <class JacobianType = DerivativeRangeReturnType>
    static JacobianType local_jacobian(const DomainType& local_x)
    {
      JacobianType ret;
      for (size_t rr = 0; rr < r; ++rr)
        for (size_t ii = 0; ii < rC; ++ii)
          ret[rr][ii][0] = 2 * local_x[0];
      return ret;
    }

    template <class JacobianType = DerivativeRangeReturnType>
    static JacobianType expected_jacobian(const ElementType& element, const DomainType& local_x)
    {
      JacobianType ret;
      const auto J_inv_T = element.geometry().jacobianInverseTransposed(local_x);
      for (size_t rr = 0; rr < r; ++rr)
        for (size_t dd = 0; dd < d; ++dd)
          for (size_t ii = 0; ii < rC; ++ii)
            ret[rr][ii][dd] = J_inv_T[dd][0] * 2 * local_x[0];
      return ret;
    }
  };

  template<bool anything>
  struct JacobianHelper<1, anything>
  {
    template <class JacobianType = DerivativeRangeReturnType>
    static JacobianType local_jacobian(const DomainType& local_x)
    {
      JacobianType ret;
      for (size_t rr = 0; rr < r; ++rr)
        ret[rr][0] = 2 * local_x[0];
      return ret;
    }

    template <class JacobianType = DerivativeRangeReturnType>
    static JacobianType expected_jacobian(const ElementType& element, const DomainType& local_x)
    {
      JacobianType ret;
      const auto J_inv_T = element.geometry().jacobianInverseTransposed(local_x);
      for (size_t rr = 0; rr < r; ++rr)
        for (size_t dd = 0; dd < d; ++dd)
          ret[rr][dd] = J_inv_T[dd][0] * 2 * local_x[0];
      return ret;
    }
  };

  const Grid::GridProvider<GridType> grid_;
};


TEST_F(GenericFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_constructible)
{
  size_t element_index = 0;
  const auto leaf_view = grid_.leaf_view();

  GenericType  function(
      /*order=*/0,
      [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
      [&](const auto& /*xx*/, const auto& /*param*/) { return RangeReturnType(element_index); },
      /*parameter=*/{},
      "element_index_",
      [](const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeReturnType(); },
      [](const auto& /*alpha*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeReturnType(); });

  GenericType function2(
      /*order_func=*/[](const auto& /*param*/) { return 3; },
      [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
      [&](const auto& /*xx*/, const auto& /*param*/) { return RangeReturnType(element_index); },
      /*parameter=*/{},
      "element_index_",
      [](const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeReturnType(); },
      [](const auto& /*alpha*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeReturnType(); });

  // check that defaults are compiling
  GenericType  function3(/*order=*/0);

  GenericType  function4(/*order_func=*/[](const auto& /*param*/) { return 3; });
}

TEST_F(GenericFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  size_t element_index = 0;
  const auto leaf_view = grid_.leaf_view();

  GenericType  function(
      /*order=*/0,
      [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
      [&](const auto& /*xx*/, const auto& /*param*/) { return RangeReturnType(element_index); },
      /*parameter=*/{},
      "element_index_",
      [](const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeReturnType(); },
      [](const auto& /*alpha*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeReturnType(); });
  auto local_f = function.local_function();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
  }
}

TEST_F(GenericFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_order)
{
  size_t element_index = 0;
  const auto leaf_view = grid_.leaf_view();

  for (auto vv : {1, 3, 5}) {
    const int expected_order = vv;
    GenericType  function(
        /*order=*/vv,
        [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
        [&](const auto& /*xx*/, const auto& /*param*/) { return RangeReturnType(element_index); },
        /*parameter=*/{},
        "element_index_",
        [](const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeReturnType(); },
        [](const auto& /*alpha*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeReturnType(); });

    auto local_f = function.local_function();
    for (auto&& element : Dune::elements(leaf_view)) {
      local_f->bind(element);
      const auto actual_order = local_f->order();
      EXPECT_EQ(expected_order, actual_order);
    }
  }
}

TEST_F(GenericFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_evaluate)
{
  size_t element_index = 0;
  const auto leaf_view = grid_.leaf_view();

  GenericType  function(
      /*order=*/0,
      [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
      [&](const auto& /*xx*/, const auto& /*param*/) { return RangeReturnType(element_index); },
      /*parameter=*/{},
      "element_index_",
      [](const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeReturnType(); },
      [](const auto& /*alpha*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeReturnType(); });
  auto local_f = function.local_function();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      Common::ParameterType param("power", 1);
      RangeReturnType expected_value(leaf_view.indexSet().index(element));
      const auto actual_value = local_f->evaluate(local_x);
      EXPECT_EQ(expected_value, actual_value);
    }
  }
}

TEST_F(GenericFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_jacobian)
{
  const auto leaf_view = grid_.leaf_view();

  GenericType  function(
      /*order=*/2,
      [](const auto&) {},
      [](const auto& xx, const auto& /*param*/) { return RangeReturnType(xx[0] * xx[0]); },
      /*parameter=*/{},
      "x^2",
      [](const auto& xx, const auto& /*param*/) { return JacobianHelper<>::local_jacobian(xx); },
      [](const auto& /*alpha*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeReturnType(); });
  auto local_f = function.local_function();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto expected_jacobian = JacobianHelper<>::expected_jacobian(element, local_x);
      const auto actual_jacobian = local_f->jacobian(local_x);
      EXPECT_EQ(expected_jacobian, actual_jacobian);
    }
  }
}

{% endfor  %}
