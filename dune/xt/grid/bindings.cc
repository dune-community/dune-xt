// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <string>
#include <vector>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/common/configuration.pbh>
#include <dune/xt/common/fvector.pbh>

#include "grids.hh"
#include "layers.hh"
#include "intersection.hh"
#include "type_traits.hh"

#include "boundaryinfo.bindings.hh"
#include "gridprovider.pbh"
#include "walker.pbh"
#include "walker/apply-on.bindings.hh"

namespace py = pybind11;
using namespace pybind11::literals;


template <class G>
void addbind_for_Grid(py::module& m, const std::string& grid_id)
{
  using namespace Dune::XT;
  using namespace Dune::XT::Grid;

  typedef typename Layer<G, Layers::level, Backends::view>::type LevelView;
  typedef typename Layer<G, Layers::leaf, Backends::view>::type LeafView;
#if HAVE_DUNE_FEM
  typedef typename Layer<G, Layers::level, Backends::part>::type LevelPart;
  typedef typename Layer<G, Layers::leaf, Backends::part>::type LeafPart;
#endif

  bind_GridProvider<G>(m, grid_id);
  bind_make_cube_grid<G>(m, grid_id);

  bind_DdSubdomainsGridProvider<G>(m, grid_id);
  bind_make_cube_dd_subdomains_grid<G>(m, grid_id);

#define BIND(V, id) bind_Walker<V>(m, grid_id + "_" + id);

  BIND(LevelView, "level_view");
  BIND(LeafView, "leaf_view");
#if HAVE_DUNE_FEM
  BIND(LevelPart, "level_part");
  BIND(LeafPart, "leaf_part");
#endif
#undef BIND
} // ... addbind_for_Grid(...)


PYBIND11_PLUGIN(_grid)
{
  py::module m("_grid", "dune-xt-grid");

  py::module::import("dune.xt.common");

  addbind_for_Grid<Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>>(m, "2d_cube_yaspgrid");
#if HAVE_DUNE_ALUGRID
  addbind_for_Grid<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>(m, "2d_simplex_aluconform");
#endif
  //#if HAVE_DUNE_UGGRID || HAVE_UG
  //  addbind_for_Grid<Dune::UGGrid<2>>(m, "2d_simplex_uggrid");
  //#endif
  //#if HAVE_ALBERTA
  //  addbind_for_Grid<Dune::AlbertaGrid<2, 2>>(m, "2d_simplex_albertagrid");
  //#endif

  DUNE_XT_GRID_BOUNDARYINFO_BIND(m);
  DUNE_XT_GRID_WALKER_APPLYON_BIND(m);

  m.def("init_logger",
        [](const ssize_t max_info_level,
           const ssize_t max_debug_level,
           const bool enable_warnings,
           const bool enable_colors,
           const std::string& info_color,
           const std::string& debug_color,
           const std::string& warning_color) {
          Dune::XT::Common::TimedLogger().create(
              max_info_level, max_debug_level, enable_warnings, enable_colors, info_color, debug_color, warning_color);
        },
        "max_info_level"_a = -1,
        "max_debug_level"_a = -1,
        "enable_warnings"_a = true,
        "enable_colors"_a = true,
        "info_color"_a = "blue",
        "debug_color"_a = "darkgray",
        "warning_color"_a = "red");

  return m.ptr();
}

#endif // HAVE_DUNE_PYBINDXI
