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

#include <dune/xt/functions/indicator.hh>
#include <dune/xt/functions/grid-function.hh>
#include <dune/xt/functions/visualization.hh>

using namespace Dune::XT;
using namespace Dune::XT::Functions;


{% for GRIDNAME, GRID, r, rC in config['types'] %}

struct IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static constexpr size_t d = GridType::dimension;
  static constexpr size_t r = {{r}};
  static constexpr size_t rC = {{rC}};

  using FunctionType = IndicatorFunction<d, r, rC>;

  using RangeType = typename FunctionType::RangeReturnType;
  using DomainType = typename FunctionType::DomainType;
  using DerivativeRangeType = typename FunctionType::DerivativeRangeReturnType;

  IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}()
    : grid_(Grid::make_cube_grid<GridType>(-1., 1., 4))
  {
  }

  const Grid::GridProvider<GridType> grid_;
};



TEST_F(IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_constructible)
{
  //first constructor
  DomainType lower_left(-1.);
  DomainType middle_1(-0.5);
  DomainType middle_2(-0.25);
  DomainType upper_right(0.5);
  RangeType first_value(1.);
  RangeType second_value(2.);

  FunctionType function_first_ctor({std::make_tuple(middle_1, upper_right, first_value)});

  FunctionType function_first_ctor_overlap({std::make_tuple(lower_left, middle_2, first_value), std::make_tuple(middle_1, upper_right, second_value)});

  //second constructor (the following functions are equivalent to the ones before (test see below))
  Common::FieldMatrix<double, d, 2> domains(0);
  Common::FieldMatrix<double, d, 2> overlap_domains_first(0);
  Common::FieldMatrix<double, d, 2> overlap_domains_second(0);

  for (size_t dd = 0; dd < d; ++dd) {
    domains[dd][0] = lower_left[dd];
    domains[dd][1] = upper_right[dd];
    overlap_domains_first[dd][0] = lower_left[dd];
    overlap_domains_first[dd][1] = middle_2[dd];
    overlap_domains_second[dd][0] = middle_1[dd];
    overlap_domains_second[dd][1] = upper_right[dd];
  }

  FunctionType function_second_ctor({std::make_pair(domains, first_value)});

  FunctionType function_second_ctor_overlap({std::make_pair(overlap_domains_first, first_value), std::make_pair(overlap_domains_second, second_value)});
}


TEST_F(IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, has_default_config)
{
  auto cfg = FunctionType::defaults();
  EXPECT_EQ(cfg.get<std::string>("name"), FunctionType::static_id());
}

TEST_F(IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_creatable)
{
  // auto default_function = FunctionType::create();

  // const auto leaf_view = grid_.leaf_view();
  // auto local_f = default_function->local_function();
  // for (auto&& element : elements(leaf_view)) {
  //   local_f->bind(element);
  //   const auto actual_order = local_f->order();
  //   EXPECT_EQ(0, actual_order);
  // }
}

TEST_F(IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_visualizable)
{
  // NOTE: Visualization of the indicators does only make sense for a suitable fine mesh !

  const auto leaf_view = grid_.leaf_view();

  DomainType lower_left(-1.);
  DomainType middle_1(-0.5);
  DomainType middle_2(-0.25);
  DomainType upper_right(0.5);
  RangeType first_value(1.);
  RangeType second_value(2.);

  FunctionType function({std::make_tuple(middle_1, upper_right, first_value)});
  visualize(function, leaf_view, "test__IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");

  FunctionType function_overlap({std::make_tuple(lower_left, middle_2, first_value), std::make_tuple(middle_1, upper_right, second_value)});
  function_overlap.visualize(leaf_view, "test__IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__with_overlap__is_visualizable");
}

TEST_F(IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_order)
{
  DomainType lower_left(-0.5);
  DomainType upper_right(0.5);
  RangeType first_value(1.);

  FunctionType function({std::make_tuple(lower_left, upper_right, first_value)});
  const int expected_order = 0;
  const auto actual_order = function.order();
  EXPECT_EQ(expected_order, actual_order);
}

TEST_F(IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_evaluate)
{
  const auto leaf_view = grid_.leaf_view();
  // single indicator
  for (auto ll : {-1., -0.75, -0.5, -0.25}) {
    DomainType lower_left(ll);
    for (auto ur : {1., 0.75, 0.5, 0.25}) {
      DomainType upper_right(ur);
      RangeType value(1.);
      FunctionType function({std::make_tuple(lower_left, upper_right, value)});
      for (auto&& element : elements(leaf_view)) {
        for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
          // expected
          RangeType expected_value(0.);
          DomainType point = element.geometry().global(quadrature_point.position());
          if (Common::FloatCmp::le(lower_left, point) && Common::FloatCmp::le(point, upper_right)) {
            expected_value = value;
          }
          const auto actual_value = function.evaluate(element.geometry().global(quadrature_point.position()));
          EXPECT_EQ(expected_value, actual_value);
        }
      }
    }
  }
}

TEST_F(IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_jacobian)
{
  const auto leaf_view = grid_.leaf_view();
  const DerivativeRangeType expected_jacobian;
  for (auto ll : {-1., -0.75, -0.5, -0.25}) {
    DomainType lower_left(ll);
    for (auto ur : {1., 0.75, 0.5, 0.25}) {
      DomainType upper_right(ur);
      for (auto vv : {1., 2., 3., 4.}) {
        RangeType value(vv);
        FunctionType function({std::make_tuple(lower_left, upper_right, value)});
        for (auto&& element : elements(leaf_view)) {
          for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
            const auto actual_jacobian = function.jacobian(element.geometry().global(quadrature_point.position()));
            EXPECT_EQ(expected_jacobian, actual_jacobian);
          }
        }
      }
    }
  }
}

TEST_F(IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  DomainType lower_left(-1.);
  DomainType upper_right(0.5);
  RangeType first_value(1.);

  FunctionType default_function({std::make_tuple(lower_left, upper_right, first_value)});
  auto localizable_function = make_grid_function<ElementType>(default_function);
  auto local_f = localizable_function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : elements(leaf_view)) {
    local_f->bind(element);
  }
}

TEST_F(IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_order)
{
  const auto leaf_view = grid_.leaf_view();
  const int expected_order = 0;
  for (auto ll : {-1., -0.75, -0.5, -0.25}) {
    DomainType lower_left(ll);
    for (auto ur : {1., 0.75, 0.5, 0.25}) {
      DomainType upper_right(ur);
      for (auto vv : {1., 2., 3., 4.}) {
        RangeType value(vv);
        FunctionType function({std::make_tuple(lower_left, upper_right, value)});
        auto localizable_function = make_grid_function<ElementType>(function);
        auto local_f = localizable_function.local_function();
        for (auto&& element : elements(leaf_view)) {
          local_f->bind(element);
          const auto actual_order = local_f->order();
          EXPECT_EQ(expected_order, actual_order);
        }
      }
    }
  }
}


TEST_F(IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_evaluate)
{
  const auto leaf_view = grid_.leaf_view();
  // single indicator
  for (auto ll : {-1., -0.75, -0.5, -0.25}) {
    DomainType lower_left(ll);
    for (auto ur : {1., 0.75, 0.5, 0.25}) {
      DomainType upper_right(ur);
      for (auto vv : {1., 2., 3., 4.}) {
        RangeType value(vv);
        FunctionType function({std::make_tuple(lower_left, upper_right, value)});
        auto localizable_function = make_grid_function<ElementType>(function);
        auto local_f = localizable_function.local_function();
        for (auto&& element : elements(leaf_view)) {
          local_f->bind(element);
          for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
            // expected
            RangeType expected_value(0.);
            DomainType point = element.geometry().global(quadrature_point.position());
            if (Common::FloatCmp::le(lower_left, point) && Common::FloatCmp::le(point, upper_right))
              expected_value = value;
            const auto local_x = quadrature_point.position();
            const auto actual_value = local_f->evaluate(local_x);
            EXPECT_EQ(expected_value, actual_value);
          }
        }
      }
    }
  }
  // multiple indicators
  DomainType lower_left(-1.);
  DomainType middle(-0.25);
  DomainType upper_right(0.5);
  RangeType first_value(1.);
  RangeType second_value(2.);
  FunctionType function_multiple({std::make_tuple(lower_left, middle, first_value), std::make_tuple(middle, upper_right, second_value)});
  auto localizable_function_mult = make_grid_function<ElementType>(function_multiple);
  auto local_f_mult = localizable_function_mult.local_function();

  for (auto&& element : elements(leaf_view)) {
    local_f_mult->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      // expected
      RangeType expected_value(0.);
      DomainType point = element.geometry().global(quadrature_point.position());
      if (Common::FloatCmp::le(lower_left, point) && Common::FloatCmp::le(point, middle))
        expected_value += first_value;
      if (Common::FloatCmp::le(middle, point) && Common::FloatCmp::le(point, upper_right))
        expected_value += second_value;
      const auto local_x = quadrature_point.position();
      const auto actual_value = local_f_mult->evaluate(local_x);
      EXPECT_EQ(expected_value, actual_value);
    }
  }
  // overlapping indicators
  DomainType lower_left_ol(-0.5);
  DomainType upper_right_ol(0);
  FunctionType function_overlap({std::make_tuple(lower_left, middle, first_value),
                                 std::make_tuple(middle, upper_right, second_value),
                                 std::make_tuple(lower_left_ol, upper_right_ol, first_value)});
  auto localizable_function_ol = make_grid_function<ElementType>(function_overlap);
  auto local_f_ol = localizable_function_ol.local_function();
  for (auto&& element : elements(leaf_view)) {
    // actual
    local_f_ol->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      // expected
      RangeType expected_value(0.);
      DomainType point = element.geometry().global(quadrature_point.position());
      if (Common::FloatCmp::le(lower_left, point) && Common::FloatCmp::le(point, middle))
        expected_value += first_value;
      if (Common::FloatCmp::le(middle, point) && Common::FloatCmp::le(point, upper_right))
        expected_value += second_value;
      if (Common::FloatCmp::le(lower_left_ol, point) && Common::FloatCmp::le(point, upper_right_ol))
        expected_value += first_value; // overlapping indicator
      const auto local_x = quadrature_point.position();
      const auto actual_value = local_f_ol->evaluate(local_x);
      EXPECT_EQ(expected_value, actual_value);
    }
  }
}



TEST_F(IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_jacobian)
{
  const auto leaf_view = grid_.leaf_view();
  const DerivativeRangeType expected_jacobian;
  for (auto ll : {-1., -0.75, -0.5, -0.25}) {
    DomainType lower_left(ll);
    for (auto ur : {1., 0.75, 0.5, 0.25}) {
      DomainType upper_right(ur);
      for (auto vv : {1., 2., 3., 4.}) {
        RangeType value(vv);
        FunctionType function({std::make_tuple(lower_left, upper_right, value)});
        auto localizable_function = make_grid_function<ElementType>(function);
        auto local_f = localizable_function.local_function();
        for (auto&& element : elements(leaf_view)) {
          local_f->bind(element);
          for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
            const auto local_x = quadrature_point.position();
            const auto actual_jacobian = local_f->jacobian(local_x);
            EXPECT_EQ(expected_jacobian, actual_jacobian);
          }
        }
      }
    }
  }
}

{% endfor  %}
