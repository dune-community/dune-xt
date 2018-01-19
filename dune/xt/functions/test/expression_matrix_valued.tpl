#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/geometry/quadraturerules.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/expression.hh>

using namespace Dune::XT;

{% for GRIDNAME, GRID, r, rC in config['types'] %}

struct ExpressionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static const constexpr size_t d = GridType::dimension;
  static const size_t r = {{r}};
  static const size_t rC = {{rC}};

  using FunctionType = Dune::XT::Functions::ExpressionFunction<d, r, rC>;

  using RangeType = typename FunctionType::RangeType;
  using DomainType = typename FunctionType::DomainType;
  using DerivativeRangeType = typename FunctionType::DerivativeRangeType;

  using RangeExpressionType = typename Dune::XT::Common::FieldMatrix<std::string, r, rC>;
  using DerivativeRangeExpressionType =
      typename Dune::XT::Common::FieldVector<Dune::XT::Common::FieldMatrix<std::string, rC, d>, r>;
  using DerivativeRangeExpressionSingleType = typename Dune::XT::Common::FieldMatrix<std::string, rC, d>;

  ExpressionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}()
    : grid_(Dune::XT::Grid::make_cube_grid<GridType>(-1., 1., 4))
  {
  }

  const Dune::XT::Grid::GridProvider<GridType> grid_;
};


TEST_F(ExpressionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_constructible)
{
  // construct a function without providing its gradient
  RangeExpressionType expr_1(std::string("x[0]*x[0]"));
  FunctionType function("x", expr_1, 2);

  // construct function with its gradient
  DerivativeRangeExpressionSingleType grad_1_single(std::string("0"));
  DerivativeRangeExpressionSingleType grad_1_single_specified(std::string("2*x[0]"));
  DerivativeRangeExpressionType grad_1(grad_1_single);
  grad_1[0] = grad_1_single_specified;
  FunctionType function_grad("x", expr_1, grad_1, 2);

  // construct a second function
  RangeExpressionType expr_2(std::string(""));
  DerivativeRangeExpressionSingleType grad_2_single(std::string("0"));
  DerivativeRangeExpressionType grad_2(grad_2_single);
  for (size_t rr = 0; rr < r; ++rr) {
      for (size_t rc = 0; rc < rC; ++rc) {
          expr_2[rr][rc] = "exp(x[0])+sin(x[0])+x[0]";
          grad_2[rr][rc][0] = "exp(x[0])-cos(x[0])+1";
      }
  }

  FunctionType second_function("x", expr_2, 4);
  FunctionType second_function_grad("x", expr_2, grad_2, 4);
}


TEST_F(ExpressionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, has_default_config)
{
  auto cfg = FunctionType::default_config();
  EXPECT_EQ(cfg.get<std::string>("type"), FunctionType::static_id());
}

TEST_F(ExpressionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_creatable)
{
  auto default_function = FunctionType::create();
  EXPECT_EQ(3, default_function->order());
}

TEST_F(ExpressionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_visualizable)
{
  const auto leaf_view = grid_.leaf_view();

  RangeExpressionType expr_1(std::string("x[0]*x[0]"));
  FunctionType function("x", expr_1, 2);

  RangeExpressionType expr_2(std::string(""));
  for (size_t rr = 0; rr < r; ++rr) {
      for (size_t rc = 0; rc < rC; ++rc)
        expr_2[rr][rc] = "exp(x[0])+sin(x[0])+x[0]";
  }
  FunctionType second_function("x", expr_2, 4);

  function.visualize(leaf_view, "test__ExpressionFunction_1_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__function");
  second_function.visualize(leaf_view, "test__ExpressionFunction_2_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__id");
}

TEST_F(ExpressionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_order)
{
  const int expected_order = 2;

  RangeExpressionType expr(std::string("x[0]*x[0]"));
  FunctionType function("x", expr, 2);

  const auto actual_order_function = function.order();

  EXPECT_EQ(expected_order, actual_order_function);
}

TEST_F(ExpressionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_evaluate)
{
  // Constant functions
  for (auto value : {-10., 3., 17., 41.}) {
    const RangeType expected_value(value);
    const RangeExpressionType constant_expr(Common::to_string(value));
    FunctionType function("x", constant_expr, 0);
    for (auto point : {-1., -0.5, 0., 0.5, 1.}) {
      const DomainType xx(point);
      const auto actual_value = function.evaluate(xx);
      EXPECT_EQ(expected_value, actual_value);
    }
  }

  // Further functions
  RangeExpressionType expr_2(std::string(""));
  DerivativeRangeExpressionSingleType grad_2_single(std::string("0"));
  DerivativeRangeExpressionType grad_2(grad_2_single);
  for (size_t rr = 0; rr < r; ++rr) {
      for (size_t rc = 0; rc < rC; ++rc)
        expr_2[rr][rc] = "exp(x[0])+sin(x[0])+x[0]";
  }
  FunctionType second_function("x", expr_2, 4);

  for (auto point : {-1., -0.5, 0., 0.5, 1.}) {
      const DomainType xx(point);
      RangeType expected_value(exp(point) + sin(point) + point);
      const auto actual_value = second_function.evaluate(xx);
      EXPECT_EQ(expected_value, actual_value);
    }
}


TEST_F(ExpressionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_jacobian)
{
  // Constant functions
  for (auto value : {-10., 3., 17., 41.}) {
    const RangeExpressionType constant_expr(Common::to_string(value));
    DerivativeRangeExpressionSingleType constant_grad_single;
    DerivativeRangeExpressionType constant_grad(constant_grad_single);
    for (size_t rr = 0; rr < r; ++rr) {
      for (size_t dd = 0; dd < d; ++dd)
        constant_grad[rr][dd] = "0";
    }
    DerivativeRangeType expected_jacobian(Common::FieldMatrix<double, rC, d>(0.));
    FunctionType function("x", constant_expr, constant_grad, 0);
    for (auto point : {-1., -0.5, 0., 0.5, 1.}) {
      const DomainType xx(point);
      const auto actual_jacobian = function.jacobian(xx);
      EXPECT_EQ(expected_jacobian, actual_jacobian);
    }
  }

  // construct a second function
  RangeExpressionType expr_2(std::string(""));
  DerivativeRangeExpressionSingleType grad_2_single(std::string("0"));
  DerivativeRangeExpressionType grad_2(grad_2_single);
  for (size_t rr = 0; rr < r; ++rr) {
      for (size_t rc = 0; rc < rC; ++rc) {
          expr_2[rr][rc] = "exp(x[0])+sin(x[0])+x[0]";
          grad_2[rr][rc][0] = "exp(x[0])-cos(x[0])+1";
      }
  }
  FunctionType second_function("x", expr_2, grad_2, 4);
  for (auto point : {-1., -0.5, 0., 0.5, 1.}) {
    DerivativeRangeType expected_jacobian(Common::FieldMatrix<double, rC, d>(0.));
    const DomainType xx(point);
    for (size_t rr = 0; rr < r; ++rr) {
        for (size_t rc = 0; rc < rC; ++rc)
          expected_jacobian[rr][rc][0] = exp(point)-cos(point)+1;
    }
    const auto actual_jacobian = second_function.jacobian(xx);
    EXPECT_EQ(expected_jacobian, actual_jacobian);
  }
}


TEST_F(ExpressionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_localizable)
{
  auto default_function = FunctionType::create();
  const auto leaf_view = grid_.leaf_view();

  const auto& localizable_function = default_function->template as_localizable<ElementType>();

  for (auto&& element : Dune::elements(leaf_view)) {
    const auto local_f = localizable_function.local_function(element);
  }
}

TEST_F(ExpressionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  auto default_function = FunctionType::create();

  const auto& localizable_function = default_function->template as_localizable<ElementType>();

  auto local_f = localizable_function.local_function();

  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
  }
}

TEST_F(ExpressionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_order)
{
  const int expected_order = 2;

  RangeExpressionType expr(std::string("x[0]*x[0]"));
  FunctionType function("x", expr, 2);

  const auto leaf_view = grid_.leaf_view();

  const auto& localizable_function = function.template as_localizable<ElementType>();
  auto local_f = localizable_function.local_function();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    const auto actual_order = local_f->order();
    EXPECT_EQ(expected_order, actual_order);
  }
}

TEST_F(ExpressionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_evaluate)
{
  const auto leaf_view = grid_.leaf_view();

  // Constant functions
  for (auto value : {-10., 3., 17., 41.}) {
    const RangeType expected_value(value);
    const RangeExpressionType constant_expr(Common::to_string(value));
    FunctionType function("x", constant_expr, 0);
    const auto& localizable_function = function.template as_localizable<ElementType>();
    auto local_f = localizable_function.local_function();
    for (auto&& element : Dune::elements(leaf_view)) {
      local_f->bind(element);
      for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
        const auto local_x = quadrature_point.position();
        const auto actual_value = local_f->evaluate(local_x);
        EXPECT_EQ(expected_value, actual_value);
      }
    }
  }

}


TEST_F(ExpressionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_jacobian)
{
  const auto leaf_view = grid_.leaf_view();

  for (auto value : {-10., 3., 17., 41.}) {
    const RangeExpressionType constant_expr(Common::to_string(value));
    DerivativeRangeExpressionSingleType constant_grad_single;
    DerivativeRangeExpressionType constant_grad(constant_grad_single);
    for (size_t rr = 0; rr < r; ++rr) {
      for (size_t rc = 0; rc < rC; ++rc) {
        for (size_t dd = 0; dd < d; ++dd)
          constant_grad[rr][rc][dd] = "0";
      }
    }
    DerivativeRangeType expected_jacobian(Common::FieldMatrix<double, rC, d>(0.));
    FunctionType function("x", constant_expr, constant_grad, 0);
    const auto& localizable_function = function.template as_localizable<ElementType>();
    auto local_f = localizable_function.local_function();
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

{% endfor  %}

