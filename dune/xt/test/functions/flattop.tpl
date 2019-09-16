#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/grid/grids.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/flattop.hh>

using namespace Dune::XT;

{% for GRIDNAME, GRID, r, rC in config['types'] %}


struct FlattopFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static const constexpr size_t d = GridType::dimension;
  static const size_t r = {{r}};
  static const size_t rC = {{rC}};

  using FunctionType = Functions::FlatTopFunction<d, r, rC>;

  using RangeReturnType = typename FunctionType::RangeReturnType;
  using DomainType = typename FunctionType::DomainType;
  using DerivativeRangeReturnType = typename FunctionType::DerivativeRangeReturnType;

  FlattopFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}()
    : grid_(Dune::XT::Grid::make_cube_grid<GridType>(-2., 2., 4))
  {
    grid_.visualize("grid");
  }

  const Grid::GridProvider<GridType> grid_;
};


TEST_F(FlattopFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_constructible)
{
  const DomainType left(0);
  const DomainType right(1);
  const DomainType delta(1e-6);
  const double top_value = 20;
  FunctionType function(left, right, delta, top_value);
}


TEST_F(FlattopFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, has_default_config)
{
  auto cfg = FunctionType::defaults();
  EXPECT_EQ(cfg.get<std::string>("name"), FunctionType::static_id());
}


TEST_F(FlattopFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_creatable)
{
  // auto default_function = FunctionType::create();
  // EXPECT_EQ(default_function->order(), 3 * d);
}


TEST_F(FlattopFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_visualizable)
{
  const auto leaf_view = grid_.leaf_view();
  const DomainType left(0);
  const DomainType right(1);
  const DomainType delta(1e-6);
  const double top_value = 20;
  FunctionType function(left, right, delta, top_value);
  function.visualize(leaf_view, "test__FlattopFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");
}


TEST_F(FlattopFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_order)
{
  const int expected_order = 3 * d;
  const DomainType left(0);
  const DomainType right(1);
  const DomainType delta(1e-6);
  const double top_value = 20;
  FunctionType function(left, right, delta, top_value);
  const auto actual_order = function.order();
  EXPECT_EQ(expected_order, actual_order);
}


TEST_F(FlattopFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, global_evaluate)
{
  for (auto value : {-10., 3., 17., 41.}) {
    const DomainType left(0);
    const DomainType right(1);
    const DomainType delta(1e-6);
    const RangeReturnType expected_value(value);
    FunctionType function(left, right, delta, value);
    for (auto point : {0.25, 0.5, 0.75}) {
      const DomainType xx(point);
      const auto actual_value = function.evaluate(xx);
      EXPECT_EQ(expected_value, actual_value);
    }
  }
}

TEST_F(FlattopFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  const DomainType left(0);
  const DomainType right(1);
  const DomainType delta(1e-6);
  const double top_value = 20;
  FunctionType function(left, right, delta, top_value);
  const auto& localizable_function = function.template as_grid_function<ElementType>();
  auto local_f = localizable_function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
  }
}


TEST_F(FlattopFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_order)
{
  const int expected_order = 3 * d;
  const DomainType left(0);
  const DomainType right(1);
  const DomainType delta(1e-6);
  const double top_value = 20;
  FunctionType function(left, right, delta, top_value);
  const auto& localizable_function = function.template as_grid_function<ElementType>();
  auto local_f = localizable_function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    const auto actual_order = local_f->order();
    EXPECT_EQ(expected_order, actual_order);
  }
}


TEST_F(FlattopFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_evaluate)
{
  const DomainType left(1);
  const DomainType right(2);
  const DomainType delta(1e-6);
  const double top_value = 20;
  FunctionType func(left, right, delta, top_value);
  const auto& localizable_function = func.template as_grid_function<ElementType>();
  auto local_f = localizable_function.local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto value = local_f->evaluate(local_x);
      const auto point = element.geometry().global(local_x);
      // check
      if (Dune::XT::Common::FloatCmp::lt(point, left - delta) || Dune::XT::Common::FloatCmp::gt(point, right + delta)) {
        // outside
        EXPECT_EQ(0.0, value) << point;
      } else if (Dune::XT::Common::FloatCmp::ge(point, left + delta)
                 && Dune::XT::Common::FloatCmp::le(point, right - delta)) {
        // inside
        EXPECT_EQ(top_value, value) << point;
      } else {
        // boundary layer
        if (top_value > 0.0) {
          EXPECT_GE(top_value, value) << point[0];
          EXPECT_LE(0.0, value) << point[0];
        } else {
          EXPECT_LE(top_value, value) << point[0];
          EXPECT_GE(0.0, value) << point[0];
        }
      }
    }
  }
}

{% endfor  %}
