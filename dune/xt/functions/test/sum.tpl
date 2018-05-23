#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/grid/grids.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/base/combined.hh>
#include <dune/xt/functions/indicator.hh>

using namespace Dune::XT;

{% for GRIDNAME, GRID, r, rC in config['types'] %}

struct SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static const constexpr size_t d = GridType::dimension;
  static const size_t r = {{r}};
  static const size_t rC = {{rC}};

  using IndicatorFunctionType = Dune::XT::Functions::IndicatorFunction<ElementType, r, rC>;

  using SumFunctionType = Dune::XT::Functions::SumFunction<IndicatorFunctionType, IndicatorFunctionType>;

  using RangeType = typename IndicatorFunctionType::RangeType;
  using DomainType = typename IndicatorFunctionType::DomainType;
  using DerivativeRangeType = typename IndicatorFunctionType::LocalFunctionType::DerivativeRangeType;

  SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}()
    : grid_(Dune::XT::Grid::make_cube_grid<GridType>(-1., 1., 4))
  {
  }

  const Dune::XT::Grid::GridProvider<GridType> grid_;
};


TEST_F(SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_constructible)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  RangeType first_value(10.);
  RangeType second_value(5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  IndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  IndicatorFunctionType expected({{'{{domains_1, first_value}, {domains_2, second_value}}'}});
  SumFunctionType manual_sum(f, g);
}


TEST_F(SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, operator_works)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  RangeType first_value(10.);
  RangeType second_value(5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  IndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  const auto& sum = f + g;
}


TEST_F(SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_visualizable)
{
  const auto leaf_view = grid_.leaf_view();
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  RangeType first_value(10.);
  RangeType second_value(5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  IndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  SumFunctionType manual_sum(f, g);
  const auto& sum = f + g;

  manual_sum.visualize(leaf_view, "test__SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");
  sum.visualize(leaf_view, "test__SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");
}

TEST_F(SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_localizable)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  RangeType first_value(10.);
  RangeType second_value(5.);


  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  IndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  SumFunctionType manual_sum(f, g);
  const auto& sum = f + g;
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    const auto local_f = manual_sum.local_function(element);
    const auto local_sum = sum.local_function(element);
  }
}

TEST_F(SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  RangeType first_value(10.);
  RangeType second_value(5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  IndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  SumFunctionType manual_sum(f, g);
  const auto& sum = f + g;

  auto local_f = manual_sum.local_function();
  auto local_sum = sum.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    local_sum->bind(element);
  }
}

TEST_F(SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_order)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  RangeType first_value(10.);
  RangeType second_value(5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  IndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  SumFunctionType manual_sum(f, g);
  const auto& sum = f + g;

  const auto leaf_view = grid_.leaf_view();
  const int expected_order = 0;
  auto local_f = manual_sum.local_function();
  auto local_sum = sum.local_function();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    local_sum->bind(element);
    const auto actual_order = local_f->order();
    const auto actual_order_sum = local_sum->order();
    EXPECT_EQ(expected_order, actual_order_sum);
    EXPECT_EQ(expected_order, actual_order);
  }
}


TEST_F(SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_evaluate)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  RangeType first_value(10.);
  RangeType second_value(5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  IndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  IndicatorFunctionType expected({{'{{domains_1, first_value}, {domains_2, second_value}}'}});

  SumFunctionType manual_sum(f, g);
  const auto& sum = f + g;

  const auto leaf_view = grid_.leaf_view();
  auto local_manual = manual_sum.local_function();
  auto local_sum = sum.local_function();
  auto local_expected = expected.local_function();
  for (auto&& element : elements(leaf_view)) {
    local_manual->bind(element);
    local_sum->bind(element);
    local_expected->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto value_manual = local_manual->evaluate(local_x);
      const auto value_sum = local_sum->evaluate(local_x);
      const auto value_expected = local_expected->evaluate(local_x);
      EXPECT_EQ(value_manual, value_expected);
      EXPECT_EQ(value_sum, value_expected);
    }
  }
}

TEST_F(SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_jacobian)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  RangeType first_value(10.);
  RangeType second_value(5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  IndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  IndicatorFunctionType expected({{'{{domains_1, first_value}, {domains_2, second_value}}'}});

  SumFunctionType manual_sum(f, g);
  const auto& sum = f + g;

  const auto leaf_view = grid_.leaf_view();
  auto local_manual = manual_sum.local_function();
  auto local_sum = sum.local_function();
  auto local_expected = expected.local_function();
  for (auto&& element : elements(leaf_view)) {
    local_manual->bind(element);
    local_sum->bind(element);
    local_expected->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto value_manual = local_manual->jacobian(local_x);
      const auto value_sum = local_sum->jacobian(local_x);
      const auto value_expected = local_expected->jacobian(local_x);
      EXPECT_EQ(value_manual, value_expected);
      EXPECT_EQ(value_sum, value_expected);
    }
  }
}

{% endfor  %}
