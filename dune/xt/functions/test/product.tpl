#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/grid/grids.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/base/combined.hh>
#include <dune/xt/functions/indicator.hh>

using namespace Dune::XT;

{% for GRIDNAME, GRID, r, rC in config['types'] %}

struct ProductFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static const constexpr size_t d = GridType::dimension;
  static const size_t r = {{r}};
  static const size_t rC = {{rC}};

  using IndicatorFunctionType = Dune::XT::Functions::IndicatorFunction<ElementType, r, rC>;
  using ScalarIndicatorFunctionType = Dune::XT::Functions::IndicatorFunction<ElementType, 1, 1>;

  using ProductFunctionType = Dune::XT::Functions::ProductFunction<ScalarIndicatorFunctionType, IndicatorFunctionType>;

  using RangeType = typename IndicatorFunctionType::RangeType;
  using ScalarRangeType = typename ScalarIndicatorFunctionType::RangeType;
  using DomainType = typename IndicatorFunctionType::DomainType;
  using DerivativeRangeType = typename IndicatorFunctionType::LocalFunctionType::DerivativeRangeType;

  ProductFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}()
    : grid_(Dune::XT::Grid::make_cube_grid<GridType>(-1., 1., 4))
  {
  }

  const Dune::XT::Grid::GridProvider<GridType> grid_;
};


TEST_F(ProductFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_constructible)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);
  Common::FieldMatrix<double, d, 2> domains_product(0);

  ScalarRangeType first_value(10.);
  RangeType second_value(5.);
  RangeType product_value(50.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = 0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
    domains_product[dd][0] = 0;
    domains_product[dd][1] = 0.5;
  }

  ScalarIndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  IndicatorFunctionType expected({{'{{domains_product, product_value}}'}});
  ProductFunctionType manual_product(f, g);
}


TEST_F(ProductFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, operator_works)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  ScalarRangeType first_value(10.);
  RangeType second_value(5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  ScalarIndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  const auto& product = f * g;
}


TEST_F(ProductFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_visualizable)
{
  const auto leaf_view = grid_.leaf_view();
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  ScalarRangeType first_value(10.);
  RangeType second_value(5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  ScalarIndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  ProductFunctionType manual_product(f, g);
  const auto& product = f * g;

  manual_product.visualize(leaf_view, "test__ProductFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");
  product.visualize(leaf_view, "test__ProductFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");
}

TEST_F(ProductFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_localizable)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  ScalarRangeType first_value(10.);
  RangeType second_value(5.);


  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  ScalarIndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  ProductFunctionType manual_product(f, g);
  const auto& product = f * g;
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    const auto local_f = manual_product.local_function(element);
    const auto local_product = product.local_function(element);
  }
}

TEST_F(ProductFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  ScalarRangeType first_value(10.);
  RangeType second_value(5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  ScalarIndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  ProductFunctionType manual_product(f, g);
  const auto& product = f * g;

  auto local_f = manual_product.local_function();
  auto local_product = product.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    local_product->bind(element);
  }
}

TEST_F(ProductFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_order)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);

  ScalarRangeType first_value(10.);
  RangeType second_value(5.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = -0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
  }

  ScalarIndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  ProductFunctionType manual_product(f, g);
  const auto& product = f * g;

  const auto leaf_view = grid_.leaf_view();
  const int expected_order = 0;
  auto local_f = manual_product.local_function();
  auto local_product = product.local_function();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    local_product->bind(element);
    const auto actual_order = local_f->order();
    const auto actual_order_product = local_product->order();
    EXPECT_EQ(expected_order, actual_order_product);
    EXPECT_EQ(expected_order, actual_order);
  }
}


TEST_F(ProductFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_evaluate)
{
  Common::FieldMatrix<double, d, 2> domains_1(0);
  Common::FieldMatrix<double, d, 2> domains_2(0);
  Common::FieldMatrix<double, d, 2> domains_product(0);

  ScalarRangeType first_value(10.);
  RangeType second_value(5.);
  RangeType product_value(50.);

  for (size_t dd = 0; dd < d; ++dd) {
    domains_1[dd][0] = -1;
    domains_1[dd][1] = 0.5;
    domains_2[dd][0] = 0;
    domains_2[dd][1] = 1;
    domains_product[dd][0] = 0;
    domains_product[dd][1] = 0.5;
  }

  ScalarIndicatorFunctionType f({{'{{domains_1, first_value}}'}});
  IndicatorFunctionType g({{'{{domains_2, second_value}}'}});

  IndicatorFunctionType expected({{'{{domains_product, product_value}}'}});

  ProductFunctionType manual_product(f, g);
  const auto& product = f * g;

  const auto leaf_view = grid_.leaf_view();
  auto local_manual = manual_product.local_function();
  auto local_product = product.local_function();
  auto local_expected = expected.local_function();
  for (auto&& element : elements(leaf_view)) {
    local_manual->bind(element);
    local_product->bind(element);
    local_expected->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto value_manual = local_manual->evaluate(local_x);
      const auto value_product = local_product->evaluate(local_x);
      const auto value_expected = local_expected->evaluate(local_x);
      EXPECT_EQ(value_manual, value_expected);
      EXPECT_EQ(value_product, value_expected);
    }
  }
}


{% endfor  %}
