#include <dune/xt/common/test/main.hxx>

#include <dune/xt/grid/grids.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/indicator.hh>

using namespace Dune::XT;


{% for GRIDNAME, GRID, r, rC in config['types'] %}

struct IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static const constexpr size_t d = GridType::dimension;
  static const size_t r = {{r}};
  static const size_t rC = {{rC}};

  using FunctionType = Functions::IndicatorGridFunction<ElementType, r, rC>;

  using RangeType = typename FunctionType::RangeType;
  using DomainType = typename FunctionType::DomainType;
  using DerivativeRangeType = typename FunctionType::LocalFunctionType::DerivativeRangeType;

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

  FunctionType function_first_ctor({{'{{lower_left, upper_right, first_value}}'}});

  FunctionType function_first_ctor_overlap({{'{{lower_left, middle_2, first_value}, {middle_1, upper_right, second_value}}'}});

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

  FunctionType function_second_ctor({{'{{domains, first_value}}'}});

  FunctionType function_second_ctor_overlap({{'{{overlap_domains_first, first_value}, {overlap_domains_second, second_value}}'}});
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
  const auto leaf_view = grid_.leaf_view();

  DomainType lower_left(-1.);
  DomainType middle_1(-0.5);
  DomainType middle_2(-0.25);
  DomainType upper_right(0.5);
  RangeType first_value(1.);
  RangeType second_value(2.);

  FunctionType function({{'{{lower_left, upper_right, first_value}}'}});
  function.visualize(leaf_view, "test__IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__is_visualizable");

  FunctionType function_overlap({{'{{lower_left, middle_2, first_value}, {middle_1, upper_right, second_value}}'}});
  function_overlap.visualize(leaf_view, "test__IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}__with_overlap__is_visualizable");
}

TEST_F(IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, is_bindable)
{
  DomainType lower_left(-1.);
  DomainType upper_right(0.5);
  RangeType first_value(1.);

  FunctionType default_function({{'{{lower_left, upper_right, first_value}}'}});
  auto local_f = default_function.local_function();
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
        FunctionType function({{'{{lower_left, upper_right, value}}'}});
        auto local_f = function.local_function();
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
        FunctionType function({{'{{lower_left, upper_right, value}}'}});
        auto local_f = function.local_function();
        for (auto&& element : elements(leaf_view)) {
          // expected
          RangeType expected_value(0.);
          const auto center = element.geometry().center();
          if (Common::FloatCmp::le(lower_left, center) && Common::FloatCmp::lt(center, upper_right))
            expected_value = value;
          // actual
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
  // multiple indicators
  DomainType lower_left(-1.);
  DomainType middle(-0.25);
  DomainType upper_right(0.5);
  RangeType first_value(1.);
  RangeType second_value(2.);
  FunctionType function_multiple({{'{{lower_left, middle, first_value}, {middle, upper_right, second_value}}'}});
  auto local_f_mult = function_multiple.local_function();
  for (auto&& element : elements(leaf_view)) {
    // expected
    RangeType expected_value(0.);
    const auto center = element.geometry().center();
    if (Common::FloatCmp::le(lower_left, center) && Common::FloatCmp::lt(center, middle))
      expected_value += first_value;
    if (Common::FloatCmp::le(middle, center) && Common::FloatCmp::lt(center, upper_right))
      expected_value += second_value;
    // actual
    local_f_mult->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto actual_value = local_f_mult->evaluate(local_x);
      EXPECT_EQ(expected_value, actual_value);
    }
  }
  // overlapping indicators
  DomainType lower_left_ol(-0.5);
  DomainType upper_right_ol(0);
  FunctionType function_overlap({{'{{lower_left, middle, first_value},{middle, upper_right, second_value},{lower_left_ol, upper_right_ol, first_value}}'}});
  auto local_f_ol = function_overlap.local_function();
  for (auto&& element : elements(leaf_view)) {
    // expected
    RangeType expected_value(0.);
    const auto center = element.geometry().center();
    if (Common::FloatCmp::le(lower_left, center) && Common::FloatCmp::lt(center, middle))
      expected_value += first_value;
    if (Common::FloatCmp::le(middle, center) && Common::FloatCmp::lt(center, upper_right))
      expected_value += second_value;
    if (Common::FloatCmp::le(lower_left_ol, center) && Common::FloatCmp::lt(center, upper_right_ol))
      expected_value += first_value; // overlapping indicator
    // actual
    local_f_ol->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto actual_value = local_f_ol->evaluate(local_x);
      EXPECT_EQ(expected_value, actual_value);
    }
  }
}


TEST_F(IndicatorFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, ctors_are_equivalent)
{
  const auto leaf_view = grid_.leaf_view();

  DomainType lower_left(-1.);
  DomainType upper_right(0.5);
  RangeType first_value(1.);

  FunctionType function_first_ctor({{'{{lower_left, upper_right, first_value}}'}});

  //second constructor
  Common::FieldMatrix<double, d, 2> domains(0);

  for (size_t dd = 0; dd < d; ++dd) {
    domains[dd][0] = lower_left[dd];
    domains[dd][1] = upper_right[dd];
  }

  FunctionType function_second_ctor({{'{{domains, first_value}}'}});

  auto local_f_1 = function_first_ctor.local_function();
  auto local_f_2 = function_second_ctor.local_function();
  for (auto&& element : elements(leaf_view)) {
    local_f_1->bind(element);
    local_f_2->bind(element);
    for (const auto& quadrature_point : Dune::QuadratureRules<double, d>::rule(element.type(), 3)) {
      const auto local_x = quadrature_point.position();
      const auto value_1 = local_f_1->evaluate(local_x);
      const auto value_2 = local_f_2->evaluate(local_x);
      EXPECT_EQ(value_1, value_2);
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
        FunctionType function({{'{{lower_left, upper_right, value}}'}});
        auto local_f = function.local_function();
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
