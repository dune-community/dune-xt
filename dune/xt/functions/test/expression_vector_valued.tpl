#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/grid/grids.hh>
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

  using RangeExpressionType = typename Dune::XT::Common::FieldVector<std::string, r>;
  using DerivativeRangeExpressionType = typename Dune::XT::Common::FieldMatrix<std::string, r, d>;

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
  DerivativeRangeExpressionType grad_1(std::string("0"));
  grad_1[0][0] = "2*x[0]";
  FunctionType function_grad("x", expr_1, grad_1, 2);

  // construct a second function
  RangeExpressionType expr_2(std::string(""));
  DerivativeRangeExpressionType grad_2(std::string("0"));
  for (size_t rr = 0; rr < r; ++rr) {
      expr_2[rr] = "exp(x[0])+sin(x[0])+x[0]";
      grad_2[rr][0] = "exp(x[0])+cos(x[0])+1";
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
  for (size_t rr = 0; rr < r; ++rr)
      expr_2[rr] = "exp(x[0])+sin(x[0])+x[0]";
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

  // Identity function
  RangeExpressionType id_expr(std::string("0"));
  for (size_t rr = 0; rr < d; ++rr)
      id_expr[rr] = "x[" + Common::to_string(rr) + "]";
  FunctionType function_id_1("x", id_expr, 1);
  for (auto point : {-1., -0.5, 0., 0.5, 1.}) {
    const DomainType xx(point);
    RangeType expected(0.);
    for (size_t rr = 0; rr < d; ++rr)
        expected[rr] = point;
    const auto actual_value = function_id_1.evaluate(xx);
    EXPECT_EQ(expected, actual_value);
  }

  // Further functions
  RangeExpressionType expr_2(std::string(""));
  DerivativeRangeExpressionType grad_2(std::string("0"));
  for (size_t rr = 0; rr < r; ++rr) {
      expr_2[rr] = "exp(x[0])+sin(x[0])+x[0]";
      grad_2[rr][0] = "exp(x[0])+cos(x[0])+1";
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
    DerivativeRangeExpressionType constant_grad;
    for (size_t rr = 0; rr < r; ++rr) {
      for (size_t dd = 0; dd < d; ++dd)
        constant_grad[rr][dd] = "0";
    }
    const DerivativeRangeType expected_jacobian;
    FunctionType function("x", constant_expr, constant_grad, 0);
    for (auto point : {-1., -0.5, 0., 0.5, 1.}) {
      const DomainType xx(point);
      const auto actual_jacobian = function.jacobian(xx);
      EXPECT_EQ(expected_jacobian, actual_jacobian);
    }
  }

  // Further functions
  RangeExpressionType expr_2(std::string(""));
  DerivativeRangeExpressionType grad_2(std::string("0"));
  for (size_t rr = 0; rr < r; ++rr) {
      expr_2[rr] = "exp(x[0])+sin(x[0])+x[0]";
      grad_2[rr][0] = "exp(x[0])+cos(x[0])+1";
  }
  FunctionType second_function("x", expr_2, grad_2, 4);

  for (auto point : {-1., -0.5, 0., 0.5, 1.}) {
    const DomainType xx(point);
    DerivativeRangeType expected_jacobian_1;
    for (size_t rr = 0; rr < r; ++rr)
      expected_jacobian_1[rr][0] = exp(point) + cos(point) + 1;
    const auto actual_jacobian = second_function.jacobian(xx);
    EXPECT_EQ(expected_jacobian_1, actual_jacobian);
  }
}

TEST_F(ExpressionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_localizable)
{
  auto default_function = FunctionType::create();
  const auto leaf_view = grid_.leaf_view();

  const auto& localizable_function = default_function->template as_grid_function<ElementType>();

  for (auto&& element : Dune::elements(leaf_view)) {
    const auto local_f = localizable_function.local_function(element);
  }
}

TEST_F(ExpressionFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  auto default_function = FunctionType::create();

  const auto& localizable_function = default_function->template as_grid_function<ElementType>();

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

  const auto& localizable_function = function.template as_grid_function<ElementType>();
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
    const auto& localizable_function = function.template as_grid_function<ElementType>();
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
    DerivativeRangeExpressionType constant_grad;
    for (size_t rr = 0; rr < r; ++rr) {
      for (size_t dd = 0; dd < d; ++dd)
        constant_grad[rr][dd] = "0";
    }
    const DerivativeRangeType expected_jacobian;
    FunctionType function("x", constant_expr, constant_grad, 0);
    const auto& localizable_function = function.template as_grid_function<ElementType>();
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
