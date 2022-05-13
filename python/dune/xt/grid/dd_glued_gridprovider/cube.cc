// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tim Keil (2020 - 2021)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/dd/glued.hh>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/traits.hh>

using namespace Dune;
using namespace Dune::XT::Grid::bindings;


template <class G_M, class G_L, class local_element_type>
struct make_cube_dd_grid
{
  static void bind(pybind11::module& m)
  {
    using namespace pybind11::literals;

    m.def(
        "make_cube_dd_grid",
        [](XT::Grid::GridProvider<G_M>& macro_grid,
           const local_element_type&,
           const unsigned int num_refinements) {
          return std::make_unique<XT::Grid::DD::Glued<G_M, G_L, XT::Grid::Layers::leaf>>(
              macro_grid, num_refinements, false, true);
        },
        "macro_grid"_a,
        "local_element_type"_a,
        "num_refinements"_a = 0);
  } // ... bind(...)
}; // struct make_cube_dd_grid<...>


PYBIND11_MODULE(_grid_dd_glued_gridprovider_cube, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid._grid_dd_glued_gridprovider_provider");
  py::module::import("dune.xt.grid._grid_traits");

#if HAVE_DUNE_GRID_GLUE
  make_cube_dd_grid<YASP_2D_EQUIDISTANT_OFFSET, YASP_2D_EQUIDISTANT_OFFSET, Cube>::bind(m);
#  if HAVE_DUNE_ALUGRID
  make_cube_dd_grid<ALU_2D_SIMPLEX_CONFORMING, ALU_2D_SIMPLEX_CONFORMING, Simplex>::bind(m);
  make_cube_dd_grid<YASP_2D_EQUIDISTANT_OFFSET, ALU_2D_SIMPLEX_CONFORMING, Simplex>::bind(m);

  make_cube_dd_grid<ALU_2D_CUBE, ALU_2D_CUBE, Cube>::bind(m);
  make_cube_dd_grid<ALU_2D_CUBE, ALU_2D_SIMPLEX_CONFORMING, Simplex>::bind(m);
#  endif
#endif
}
