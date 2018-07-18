#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/grid/grids.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/base/combined_functions.hh>
#include <dune/xt/functions/constant.hh>

using namespace Dune::XT;

{% for GRIDNAME, GRID, r, rC in config['types'] %}

struct SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static const constexpr size_t d = GridType::dimension;
  static const size_t r = {{r}};
  static const size_t rC = {{rC}};

  using ConstantFunctionType = Dune::XT::Functions::ConstantFunction<d, r, rC>;

  using SumFunctionType = Dune::XT::Functions::SumFunction<ConstantFunctionType, ConstantFunctionType>;

  using RangeReturnType = typename ConstantFunctionType::RangeReturnType;
  using DomainType = typename ConstantFunctionType::DomainType;
  using DerivativeRangeReturnType = typename ConstantFunctionType::DerivativeRangeReturnType;

  SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}()
    : grid_(Dune::XT::Grid::make_cube_grid<GridType>(-1., 1., 4))
  {
  }

  const Dune::XT::Grid::GridProvider<GridType> grid_;
};


TEST_F(SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_constructible)
{
  ConstantFunctionType f(3.);
  ConstantFunctionType g(2.);

  ConstantFunctionType expected(5.);
  SumFunctionType manual_sum(f, g);
}


TEST_F(SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, operator_works)
{
  ConstantFunctionType f(3.);
  ConstantFunctionType g(2.);
  ConstantFunctionType h(8.);

  ConstantFunctionType expected(5.);
  SumFunctionType manual_sum(f, g);

  const auto& sum = f + g;
  const auto& minus = h - f;

  for (auto point : {-100., -10., 0., 10., 100.}) {
    const DomainType xx(point);
    const auto manual_value = manual_sum.evaluate(xx);
    const auto sum_value = sum.evaluate(xx);
    const auto minus_value = minus.evaluate(xx);
    EXPECT_EQ(sum_value, manual_value);
    EXPECT_EQ(minus_value, manual_value);
  }
}


//TEST_F(SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_visualizable)
//{
//  const auto leaf_view = grid_.leaf_view();
//  Common::FieldMatrix<double, d, 2> domains_1(0);
//  Common::FieldMatrix<double, d, 2> domains_2(0);

//  RangeReturnType first_value(10.);
//  RangeReturnType second_value(5.);

//  for (size_t dd = 0; dd < d; ++dd) {
//    domains_1[dd][0] = -1;
//    domains_1[dd][1] = -0.5;
//    domains_2[dd][0] = 0;
//    domains_2[dd][1] = 1;
//  }

//  ConstantFunctionType f({{'{{domains_1, first_value}}'}});
//  ConstantFunctionType g({{'{{domains_2, second_value}}'}});

//  SumFunctionType manual_sum(f, g);
//  const auto& sum = f + g;

//  manual_sum.visualize(leaf_view, "test__SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");
//  sum.visualize(leaf_view, "test__SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");
//}

//TEST_F(SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
//{
//  Common::FieldMatrix<double, d, 2> domains_1(0);
//  Common::FieldMatrix<double, d, 2> domains_2(0);

//  RangeReturnType first_value(10.);
//  RangeReturnType second_value(5.);

//  for (size_t dd = 0; dd < d; ++dd) {
//    domains_1[dd][0] = -1;
//    domains_1[dd][1] = -0.5;
//    domains_2[dd][0] = 0;
//    domains_2[dd][1] = 1;
//  }

//  ConstantFunctionType f({{'{{domains_1, first_value}}'}});
//  ConstantFunctionType g({{'{{domains_2, second_value}}'}});

//  SumFunctionType manual_sum(f, g);
//  const auto& sum = f + g;

//  auto local_f = manual_sum.local_function();
//  auto local_sum = sum.local_function();
//  const auto leaf_view = grid_.leaf_view();
//  for (auto&& element : Dune::elements(leaf_view)) {
//    local_f->bind(element);
//    local_sum->bind(element);
//  }
//}

//TEST_F(SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_order)
//{
//  Common::FieldMatrix<double, d, 2> domains_1(0);
//  Common::FieldMatrix<double, d, 2> domains_2(0);

//  RangeReturnType first_value(10.);
//  RangeReturnType second_value(5.);

//  for (size_t dd = 0; dd < d; ++dd) {
//    domains_1[dd][0] = -1;
//    domains_1[dd][1] = -0.5;
//    domains_2[dd][0] = 0;
//    domains_2[dd][1] = 1;
//  }

//  ConstantFunctionType f({{'{{domains_1, first_value}}'}});
//  ConstantFunctionType g({{'{{domains_2, second_value}}'}});

//  SumFunctionType manual_sum(f, g);
//  const auto& sum = f + g;

//  const auto leaf_view = grid_.leaf_view();
//  const int expected_order = 0;
//  auto local_f = manual_sum.local_function();
//  auto local_sum = sum.local_function();
//  for (auto&& element : Dune::elements(leaf_view)) {
//    local_f->bind(element);
//    local_sum->bind(element);
//    const auto actual_order = local_f->order();
//    const auto actual_order_sum = local_sum->order();
//    EXPECT_EQ(expected_order, actual_order_sum);
//    EXPECT_EQ(expected_order, actual_order);
//  }
//}


//TEST_F(SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_evaluate)
//{
//  Common::FieldMatrix<double, d, 2> domains_1(0);
//  Common::FieldMatrix<double, d, 2> domains_2(0);

//  RangeReturnType first_value(10.);
//  RangeReturnType second_value(5.);

//  for (size_t dd = 0; dd < d; ++dd) {
//    domains_1[dd][0] = -1;
//    domains_1[dd][1] = -0.5;
//    domains_2[dd][0] = 0;
//    domains_2[dd][1] = 1;
//  }

//  ConstantFunctionType f({{'{{domains_1, first_value}}'}});
//  ConstantFunctionType g({{'{{domains_2, second_value}}'}});

//  ConstantFunctionType expected({{'{{domains_1, first_value}, {domains_2, second_value}}'}});

//  SumFunctionType manual_sum(f, g);
//  const auto& sum = f + g;

//  const auto leaf_view = grid_.leaf_view();
//  auto local_manual = manual_sum.local_function();
//  auto local_sum = sum.local_function();
//  auto local_expected = expected.local_function();
//  for (auto&& element : elements(leaf_view)) {
//    local_manual->bind(element);
//    local_sum->bind(element);
//    local_expected->bind(element);
//    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
//      const auto local_x = quadrature_point.position();
//      const auto value_manual = local_manual->evaluate(local_x);
//      const auto value_sum = local_sum->evaluate(local_x);
//      const auto value_expected = local_expected->evaluate(local_x);
//      EXPECT_EQ(value_manual, value_expected);
//      EXPECT_EQ(value_sum, value_expected);
//    }
//  }
//}

//TEST_F(SumFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_jacobian)
//{
//  Common::FieldMatrix<double, d, 2> domains_1(0);
//  Common::FieldMatrix<double, d, 2> domains_2(0);

//  RangeReturnType first_value(10.);
//  RangeReturnType second_value(5.);

//  for (size_t dd = 0; dd < d; ++dd) {
//    domains_1[dd][0] = -1;
//    domains_1[dd][1] = -0.5;
//    domains_2[dd][0] = 0;
//    domains_2[dd][1] = 1;
//  }

//  ConstantFunctionType f({{'{{domains_1, first_value}}'}});
//  ConstantFunctionType g({{'{{domains_2, second_value}}'}});

//  ConstantFunctionType expected({{'{{domains_1, first_value}, {domains_2, second_value}}'}});

//  SumFunctionType manual_sum(f, g);
//  const auto& sum = f + g;

//  const auto leaf_view = grid_.leaf_view();
//  auto local_manual = manual_sum.local_function();
//  auto local_sum = sum.local_function();
//  auto local_expected = expected.local_function();
//  for (auto&& element : elements(leaf_view)) {
//    local_manual->bind(element);
//    local_sum->bind(element);
//    local_expected->bind(element);
//    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
//      const auto local_x = quadrature_point.position();
//      const auto value_manual = local_manual->jacobian(local_x);
//      const auto value_sum = local_sum->jacobian(local_x);
//      const auto value_expected = local_expected->jacobian(local_x);
//      EXPECT_EQ(value_manual, value_expected);
//      EXPECT_EQ(value_sum, value_expected);
//    }
//  }
//}

{% endfor  %}
