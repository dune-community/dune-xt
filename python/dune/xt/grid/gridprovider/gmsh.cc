// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/gridprovider/gmsh.hh>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/traits.hh>

using namespace Dune;
using namespace Dune::XT::Grid::bindings;


template <class G, class element_type>
struct make_gmsh_grid
{
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def(
        "make_gmsh_grid",
        [](const std::string& filename, const Dimension<d>&, const element_type&) {
          return XT::Grid::make_gmsh_grid<G>(filename);
        },
        "filename"_a,
        "dim"_a,
        "element_type"_a);
  } // ... bind(...)
}; // struct make_gmsh_grid


PYBIND11_MODULE(_grid_gridprovider_gmsh, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid._grid_gridprovider_provider");
  py::module::import("dune.xt.grid._grid_traits");

#if HAVE_DUNE_ALUGRID
  make_gmsh_grid<ALU_2D_SIMPLEX_CONFORMING, Simplex>::bind(m);
  make_gmsh_grid<ALU_2D_CUBE, Cube>::bind(m);
  make_gmsh_grid<ALU_3D_SIMPLEX_CONFORMING, Simplex>::bind(m);
  make_gmsh_grid<ALU_3D_CUBE, Cube>::bind(m);
#endif
}
