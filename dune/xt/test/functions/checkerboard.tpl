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

#include <dune/xt/functions/checkerboard.hh>

using namespace Dune::XT;


{% for GRIDNAME, GRID, r, rC in config['types'] %}

struct CheckerboardFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static constexpr size_t d = GridType::dimension;
  static constexpr size_t r = {{r}};
  static constexpr size_t rC = {{rC}};

  using FunctionType = Dune::XT::Functions::CheckerboardFunction<ElementType, r, rC>;

  using RangeType = typename FunctionType::RangeType;
  using DomainType = typename FunctionType::DomainType;
  using DerivativeRangeType = typename FunctionType::LocalFunctionType::DerivativeRangeType;

  CheckerboardFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}()
    : grid_(Dune::XT::Grid::make_cube_grid<GridType>(-1., 1., 4))
  {
  }

  const Dune::XT::Grid::GridProvider<GridType> grid_;
};


TEST_F(CheckerboardFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_constructible)
{
  DomainType lower_left(0.);
  DomainType upper_right(1.);
  Common::FieldVector<size_t, d> num_elements(2.);
  size_t num_squares = 1;
  std::vector<RangeType> values;
  for (size_t dd = 0; dd < d; ++dd)
      num_squares *= num_elements[dd];
  for (size_t ii = 0; ii < num_squares; ++ii) {
      RangeType entry(ii+1);
      values.emplace_back(entry);
  }
  FunctionType function(lower_left, upper_right, num_elements, values);
}


TEST_F(CheckerboardFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, has_default_config)
{
  auto cfg = FunctionType::defaults();
  EXPECT_EQ(cfg.get<std::string>("name"), FunctionType::static_id());
}

TEST_F(CheckerboardFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_creatable)
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

TEST_F(CheckerboardFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_visualizable)
{
  const auto leaf_view = grid_.leaf_view();

  DomainType lower_left_all(-1.);
  DomainType lower_left(0.);
  DomainType upper_right(1.);
  Common::FieldVector<size_t, d> num_elements(2.);
  size_t num_squares = 1;
  std::vector<RangeType> values;
  for (size_t dd = 0; dd < d; ++dd)
      num_squares *= num_elements[dd];
  for (size_t ii = 0; ii < num_squares; ++ii) {
      RangeType entry(ii+1);
      values.emplace_back(entry);
  }
  // part of Domain
  FunctionType function(lower_left, upper_right, num_elements, values);
  // entire Domain
  FunctionType function_all(lower_left_all, upper_right, num_elements, values);

  function.visualize(leaf_view, "test__CheckerboardFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");
  function_all.visualize(leaf_view, "test__AllDomain__CheckerboardFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");
}

TEST_F(CheckerboardFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  DomainType lower_left(0.);
  DomainType upper_right(1.);
  Common::FieldVector<size_t, d> num_elements(2.);
  size_t num_squares = 1;
  std::vector<RangeType> values;
  for (size_t dd = 0; dd < d; ++dd)
      num_squares *= num_elements[dd];
  for (size_t ii = 0; ii < num_squares; ++ii) {
      RangeType entry(ii+1);
      values.emplace_back(entry);
  }
  FunctionType default_function(lower_left, upper_right, num_elements, values);

  const auto leaf_view = grid_.leaf_view();

  auto local_f = default_function.local_function();

  for (auto&& element : Dune::elements(leaf_view))
    local_f->bind(element);
}

TEST_F(CheckerboardFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_order)
{
  const auto leaf_view = grid_.leaf_view();
  const int expected_order = 0;
  Common::FieldVector<size_t, d> num_elements(2.);
  size_t num_squares = 1;
  std::vector<RangeType> values;
  for (size_t dd = 0; dd < d; ++dd)
      num_squares *= num_elements[dd];
  for (size_t ii = 0; ii < num_squares; ++ii) {
      RangeType entry(ii+1);
      values.emplace_back(entry);
  }
  for (auto ll : {-1., -0.75, -0.5, -0.25}) {
    DomainType lower_left(ll);
    for (auto ur : {1., 0.75, 0.5, 0.25}) {
      DomainType upper_right(ur);
      FunctionType function(lower_left, upper_right, num_elements, values);
      auto local_f = function.local_function();
      for (auto&& element : Dune::elements(leaf_view)) {
        local_f->bind(element);
        const auto actual_order = local_f->order();
        EXPECT_EQ(expected_order, actual_order);
      }
    }
  }
}

TEST_F(CheckerboardFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_evaluate)
{
  const auto leaf_view = grid_.leaf_view();
  Common::FieldVector<size_t, d> num_elements(2.);
  size_t num_squares = 1;
  std::vector<RangeType> values;
  for (size_t dd = 0; dd < d; ++dd)
      num_squares *= num_elements[dd];
  for (size_t ii = 0; ii < num_squares; ++ii) {
      RangeType entry(ii+1);
      values.emplace_back(entry);
  }
  for (auto ll : {-1., -0.75, -0.5, -0.25}) {
    DomainType lower_left(ll);
    for (auto ur : {1., 0.75, 0.5, 0.25}) {
      DomainType upper_right(ur);
      FunctionType function(lower_left, upper_right, num_elements, values);
      auto local_f = function.local_function();
      for (auto&& element : Dune::elements(leaf_view)) {
        RangeType expected_value(0.);
        const auto center = element.geometry().center();
        std::vector<size_t> which_partition(d, 0);
        if (Dune::XT::Common::FloatCmp::le(lower_left, center) && Dune::XT::Common::FloatCmp::lt(center, upper_right)) {
          for (size_t dd = 0; dd < d; ++dd) {
            which_partition[dd] =
                std::min(size_t(std::floor(num_elements[dd]
                                           * ((center[dd] - lower_left[dd]) / (upper_right[dd] - lower_left[dd])))),
                         num_elements[dd] - 1);
          }
          size_t subdomain = 0;
          if (d == 1)
            subdomain = which_partition[0];
          else if (d == 2)
            subdomain = which_partition[0] + which_partition[1] * num_elements[0];
          else
            subdomain = which_partition[0] + which_partition[1] * num_elements[0] + which_partition[2] * num_elements[1] * num_elements[0];
          expected_value = values[subdomain];
        }
        local_f->bind(element);
        for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
          const auto local_x = quadrature_point.position();
          const auto actual_value = local_f->evaluate(local_x);
          EXPECT_EQ(expected_value, actual_value);
        }
      }
    }
  }
}


TEST_F(CheckerboardFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_jacobian)
{
  const auto leaf_view = grid_.leaf_view();
  Common::FieldVector<size_t, d> num_elements(2.);
  size_t num_squares = 1;
  std::vector<RangeType> values;
  for (size_t dd = 0; dd < d; ++dd)
      num_squares *= num_elements[dd];
  for (size_t ii = 0; ii < num_squares; ++ii) {
      RangeType entry(ii+1);
      values.emplace_back(entry);
  }
  const DerivativeRangeType expected_jacobian = DerivativeRangeType();
  for (auto ll : {-1., -0.75, -0.5, -0.25}) {
    DomainType lower_left(ll);
    for (auto ur : {1., 0.75, 0.5, 0.25}) {
      DomainType upper_right(ur);
      FunctionType function(lower_left, upper_right, num_elements, values);
      auto local_f = function.local_function();
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
}

{% endfor  %}
