// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2020)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/functional.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/provider.hh>


PYBIND11_MODULE(_test_grid_intersection, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;
  using namespace Dune::XT;

  using G = GRID_2D;
  using I = Grid::extract_intersection_t<typename G::LeafGridView>;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");

  m.def("call_on_each_intersection", [](Grid::GridProvider<G>& grid, std::function<void(const I&)> generic_function) {
    auto gv = grid.leaf_view();
    for (auto&& element : elements(gv))
      for (auto&& intersection : intersections(gv, element))
        generic_function(intersection);
  });
}
