// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   René Fritze    (2019)
//   Tobias Leibner (2019)

#include <dune/xt/test/main.hxx>

#include <dune/xt/grid/grids.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/generic/grid-function.hh>
#include <dune/xt/functions/visualization.hh>

using namespace Dune::XT;
using namespace Dune::XT::Functions;

{% for GRIDNAME, GRID, r, rC in config['types'] %}

struct VisualizeGenericGridFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}} : public ::testing::Test
{
  using GridType = {{GRID}};
  using ElementType = typename GridType::template Codim<0>::Entity;
  static constexpr size_t d = GridType::dimension;
  static constexpr size_t r = {{r}};
  static constexpr size_t rC = {{rC}};

  using GenericType = GenericGridFunction<ElementType, r, rC>;

  using RangeType = typename GenericType::LocalFunctionType::RangeType;
  using RangeReturnType = typename GenericType::LocalFunctionType::RangeReturnType;
  using DomainType = typename GenericType::LocalFunctionType::DomainType;
  using DerivativeRangeType = typename GenericType::LocalFunctionType::DerivativeRangeType;
  using DerivativeRangeReturnType = typename GenericType::LocalFunctionType::DerivativeRangeReturnType;


  VisualizeGenericGridFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}()
    : grid_(Grid::make_cube_grid<GridType>(-1., 1., 4))
  {
  }

  const Grid::GridProvider<GridType> grid_;
};


TEST_F(VisualizeGenericGridFunction_from_{{GRIDNAME}}_to_{{r}}_times_{{rC}}, visualize)
{
  size_t element_index = 0;
  const auto leaf_view = grid_.leaf_view();

  GenericType function(
      /*order=*/0,
      [&](const auto& element) { element_index = leaf_view.indexSet().index(element); },
      [&](const auto& /*xx*/, const auto& /*param*/) { return RangeReturnType(element_index); },
      /*parameter=*/{},
      "element_index_",
      [](const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeReturnType(); },
      [](const auto& /*alpha*/, const auto& /*xx*/, const auto& /*param*/) { return DerivativeRangeReturnType(); });

  visualize(function, leaf_view, "default");
  SumVisualizer<r, rC> sum_visualizer;
  visualize(function, leaf_view, "sum", true, Dune::VTK::appendedraw, {}, sum_visualizer);
  ComponentVisualizer<r, rC> first_comp_visualizer(0);
  visualize(function, leaf_view, "component", false, Dune::VTK::appendedraw, {}, first_comp_visualizer);
  GenericVisualizer<r, rC> sin_abs_visualizer(1, [](const int, const RangeType& val) { return std::sin(val.two_norm()); });
  visualize(function, leaf_view, "sin_abs", false, Dune::VTK::appendedraw, {}, sin_abs_visualizer);

  if (r == 1)
    visualize_gradient(function, leaf_view, "gradient");
}

{% endfor  %}
