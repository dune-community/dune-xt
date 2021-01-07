// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:

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


template <class G, class element_type>
struct make_cube_dd_grid
{
  static void bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;
    using D = typename G::ctype;
    static const size_t d = G::dimension;

    m.def(
        "make_cube_dd_grid",
        [](XT::Grid::GridProvider<G>& macro_grid, const unsigned int num_refinements) {
          return std::make_unique<XT::Grid::DD::Glued<G, G, XT::Grid::Layers::leaf>>(
              macro_grid, num_refinements, false, true);
        },
        "macro_grid"_a,
        "num_refinements"_a = 0,
        py::keep_alive<0, 1>());
  } // ... bind(...)
}; // struct make_cube_dd_grid<...>


PYBIND11_MODULE(_grid_dd_glued_gridprovider_cube, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid._grid_dd_glued_gridprovider_provider");
  py::module::import("dune.xt.grid._grid_traits");

  //  make_cube_dd_grid<ONED_1D, void>::bind(m);  // so far only 2d grid
  make_cube_dd_grid<YASP_2D_EQUIDISTANT_OFFSET, Cube>::bind(m);
//  make_cube_dd_grid<YASP_3D_EQUIDISTANT_OFFSET, Cube>::bind(m);
#if HAVE_DUNE_ALUGRID
  make_cube_dd_grid<ALU_2D_SIMPLEX_CONFORMING, Simplex>::bind(m);
  make_cube_dd_grid<ALU_2D_CUBE, Cube>::bind(m);
//  make_cube_dd_grid<ALU_3D_SIMPLEX_CONFORMING, Simplex>::bind(m);
//  make_cube_dd_grid<ALU_3D_CUBE, Cube>::bind(m);
#endif
}
