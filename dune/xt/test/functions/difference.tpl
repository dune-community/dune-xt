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

#include <dune/xt/functions/base/combined-grid-functions.hh>
#include <dune/xt/functions/indicator.hh>

using namespace Dune::XT;

{% for GRIDNAME, GRID, r, rC in config['types'] %}

struct DifferenceFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static const constexpr size_t d = GridType::dimension;
  static const size_t r = {{r}};
  static const size_t rC = {{rC}};

  using IndicatorFunctionType = Dune::XT::Functions::IndicatorGridFunction<ElementType, r, rC>;

  using DifferenceFunctionType = Dune::XT::Functions::DifferenceGridFunction<IndicatorFunctionType, IndicatorFunctionType>;

  using RangeType = typename IndicatorFunctionType::RangeType;
  using DomainType = typename IndicatorFunctionType::DomainType;
  using DerivativeRangeType = typename IndicatorFunctionType::LocalFunctionType::DerivativeRangeType;

  DifferenceFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}()
    : grid_(Dune::XT::Grid::make_cube_grid<GridType>(-1., 1., 4))
  {
  }

  const Dune::XT::Grid::GridProvider<GridType> grid_;
};


TEST_F(DifferenceFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_constructible)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  RangeType first_value(10.);
  RangeType second_value(5.);
  RangeType negative_second_value(-5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  IndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  IndicatorFunctionType expected({{'{{domains_1, first_value}, {domains_2, negative_second_value}}'}});
  DifferenceFunctionType manual_difference(f, g);
}


TEST_F(DifferenceFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, operator_works)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  RangeType first_value(10.);
  RangeType second_value(5.);
  RangeType negative_second_value(-5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  IndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  const auto& difference DUNE_UNUSED = f - g;
}


TEST_F(DifferenceFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_visualizable)
{
  const auto leaf_view = grid_.leaf_view();
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  RangeType first_value(10.);
  RangeType second_value(5.);
  RangeType negative_second_value(-5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  IndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  DifferenceFunctionType manual_difference(f, g);
  const auto& difference = f - g;

  manual_difference.visualize(leaf_view, "test__DifferenceFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");
  difference.visualize(leaf_view, "test__DifferenceFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");
}

TEST_F(DifferenceFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  RangeType first_value(10.);
  RangeType second_value(5.);
  RangeType negative_second_value(-5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  IndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  DifferenceFunctionType manual_difference(f, g);
  const auto& difference = f - g;

  auto local_f = manual_difference.local_function();
  auto local_difference = difference.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    local_difference->bind(element);
  }
}

TEST_F(DifferenceFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_order)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  RangeType first_value(10.);
  RangeType second_value(5.);
  RangeType negative_second_value(-5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  IndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  DifferenceFunctionType manual_difference(f, g);
  const auto& difference = f - g;

  const auto leaf_view = grid_.leaf_view();
  const int expected_order = 0;
  auto local_f = manual_difference.local_function();
  auto local_difference = difference.local_function();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    local_difference->bind(element);
    const auto actual_order = local_f->order();
    const auto actual_order_difference = local_difference->order();
    EXPECT_EQ(expected_order, actual_order_difference);
    EXPECT_EQ(expected_order, actual_order);
  }
}


TEST_F(DifferenceFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_evaluate)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  RangeType first_value(10.);
  RangeType second_value(5.);
  RangeType negative_second_value(-5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  IndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  IndicatorFunctionType expected({{'{{domains_1, first_value}, {domains_2, negative_second_value}}'}});

  DifferenceFunctionType manual_difference(f, g);
  const auto& difference = f - g;

  const auto leaf_view = grid_.leaf_view();
  auto local_manual = manual_difference.local_function();
  auto local_difference = difference.local_function();
  auto local_expected = expected.local_function();
  for (auto&& element : elements(leaf_view)) {
    local_manual->bind(element);
    local_difference->bind(element);
    local_expected->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto value_manual = local_manual->evaluate(local_x);
      const auto value_difference = local_difference->evaluate(local_x);
      const auto value_expected = local_expected->evaluate(local_x);
      EXPECT_EQ(value_manual, value_expected);
      EXPECT_EQ(value_difference, value_expected);
    }
  }
}

TEST_F(DifferenceFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_jacobian)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  RangeType first_value(10.);
  RangeType second_value(5.);
  RangeType negative_second_value(-5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  IndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  IndicatorFunctionType expected({{'{{domains_1, first_value}, {domains_2, negative_second_value}}'}});

  DifferenceFunctionType manual_difference(f, g);
  const auto& difference = f - g;

  const auto leaf_view = grid_.leaf_view();
  auto local_manual = manual_difference.local_function();
  auto local_difference = difference.local_function();
  auto local_expected = expected.local_function();
  for (auto&& element : elements(leaf_view)) {
    local_manual->bind(element);
    local_difference->bind(element);
    local_expected->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto value_manual = local_manual->jacobian(local_x);
      const auto value_difference = local_difference->jacobian(local_x);
      const auto value_expected = local_expected->jacobian(local_x);
      EXPECT_EQ(value_manual, value_expected);
      EXPECT_EQ(value_difference, value_expected);
    }
  }
}

{% endfor  %}
