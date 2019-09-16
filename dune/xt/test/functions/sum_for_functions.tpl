#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/xt/test/main.hxx>

#include <dune/xt/grid/grids.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/base/combined-functions.hh>
#include <dune/xt/functions/constant.hh>

/**
 * todo: Tests for combined_functions
 * currently testing only constant functions with minus and plus
 */

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


{% endfor  %}
