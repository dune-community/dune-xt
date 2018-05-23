#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/grid/grids.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/spe10/model1.hh>

using namespace Dune::XT;

{% for GRIDNAME, GRID, r, rC in config['types'] %}


struct Spe10Model1Function_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static const constexpr size_t d = GridType::dimension;
  static const size_t r = {{r}};
  static const size_t rC = {{rC}};

  using FunctionType = Functions::Spe10::Model1Function<ElementType, r, rC>;

  using DerivativeRangeType = typename FunctionType::LocalFunctionType::DerivativeRangeType;

  Spe10Model1Function_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}()
    : grid_(Dune::XT::Grid::make_cube_grid<GridType>({0, 0}, {762.0, 152.4}, {100, 20}))
  {
    grid_.visualize("grid");
  }

  const Dune::XT::Grid::GridProvider<GridType> grid_;
};


TEST_F(Spe10Model1Function_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_constructible)
{
  auto filename = Dune::XT::Functions::Spe10::internal::model1_filename;
  FunctionType function(
      filename,
      {0, 0},
      {Dune::XT::Functions::Spe10::internal::model_1_length_x, Dune::XT::Functions::Spe10::internal::model_1_length_z});
}

TEST_F(Spe10Model1Function_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, has_default_config)
{
  auto cfg = FunctionType::default_config();
  EXPECT_EQ(cfg.get<std::string>("type"), FunctionType::static_id());
}

TEST_F(Spe10Model1Function_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_creatable)
{
  auto default_function = FunctionType::create();

  const auto leaf_view = grid_.leaf_view();
  auto local_f = default_function->local_function();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    const auto actual_order = local_f->order();
    EXPECT_EQ(0, actual_order);
  }
}

TEST_F(Spe10Model1Function_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_visualizable)
{
  auto default_function = FunctionType::create();
  const auto leaf_view = grid_.leaf_view();
  default_function->visualize(leaf_view, "test__Spe10Model1Function_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");
}

TEST_F(Spe10Model1Function_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_localizable)
{
  auto default_function = FunctionType::create();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    const auto local_f = default_function->local_function(element);
  }
}

TEST_F(Spe10Model1Function_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  auto default_function = FunctionType::create();
  auto local_f = default_function->local_function();
  const auto leaf_view = grid_.leaf_view();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
  }
}

TEST_F(Spe10Model1Function_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_order)
{
  auto default_function = FunctionType::create();
  const auto leaf_view = grid_.leaf_view();
  const int expected_order = 0;
  auto local_f = default_function->local_function();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    const auto actual_order = local_f->order();
    EXPECT_EQ(expected_order, actual_order);
  }
}


TEST_F(Spe10Model1Function_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_evaluate)
{
  const auto leaf_view = grid_.leaf_view();
  auto default_function = FunctionType::create();
  auto local_f = default_function->local_function();
  double lower = Dune::XT::Functions::Spe10::internal::model1_min_value;
  double upper = Dune::XT::Functions::Spe10::internal::model1_max_value + 1;

  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto actual_value = local_f->evaluate(local_x);
      EXPECT_LE(actual_value, upper);
      EXPECT_LE(lower, actual_value);
    }
  }
}


TEST_F(Spe10Model1Function_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, local_jacobian)
{
  const auto leaf_view = grid_.leaf_view();
  const DerivativeRangeType expected_jacobian = DerivativeRangeType();
  auto default_function = FunctionType::create();
  auto local_f = default_function->local_function();
  for (auto&& element : Dune::elements(leaf_view)) {
    local_f->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto actual_jacobian = local_f->jacobian(local_x);
      EXPECT_EQ(expected_jacobian, actual_jacobian);
    }
  }
}

{% endfor  %}
