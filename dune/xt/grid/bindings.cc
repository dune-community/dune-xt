// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
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

#include "boundaryinfo.pbh"
#include "gridprovider.pbh"
#include "walker.pbh"
#include "walker/apply-on.pbh"

namespace py = pybind11;
using namespace pybind11::literals;


template <class I, bool required = true>
struct for_Grid_and_Intersection
{
  template <class G>
  static void addbind(py::module& m, const std::string& grid_id, const std::string& id)
  {
    using namespace Dune::XT::Grid;
    pybind11::class_<BoundaryInfo<I>>(
        m, std::string("BoundaryInfo__" + grid_id + id).c_str(), std::string("BoundaryInfo__" + grid_id + id).c_str());
    bind_BoundaryInfo<AllDirichletBoundaryInfo<I>>(m, "AllDirichletBoundaryInfo__" + grid_id + id);
    bind_BoundaryInfo<AllNeumannBoundaryInfo<I>>(m, "AllNeumannBoundaryInfo__" + grid_id + id);
    bind_BoundaryInfo<NormalBasedBoundaryInfo<I>>(m, "NormalBasedBoundaryInfo__" + grid_id + id);
  }
};

template <class I>
struct for_Grid_and_Intersection<I, false>
{
  template <class G>
  static void addbind(py::module& /*m*/, const std::string& /*grid_id*/, const std::string& /*id*/)
  {
  }
};


template <class G>
void addbind_for_Grid(py::module& m, const std::string& grid_id)
{
  using namespace Dune::XT;
  using namespace Dune::XT::Grid;

  typedef typename Layer<G, Layers::level, Backends::view>::type LevelView;
  typedef typename Layer<G, Layers::leaf, Backends::view>::type LeafView;
  typedef typename Intersection<LeafView>::type FVI;
  typedef typename Intersection<LevelView>::type LVI;
#if HAVE_DUNE_FEM
  typedef typename Layer<G, Layers::level, Backends::part>::type LevelPart;
  typedef typename Layer<G, Layers::leaf, Backends::part>::type LeafPart;
  typedef typename Intersection<LeafPart>::type FPI;
  typedef typename Intersection<LevelPart>::type LPI;
#endif

  bind_GridProvider<G>(m, grid_id);
  bind_make_cube_grid<G>(m, grid_id);

  for_Grid_and_Intersection<FVI, true>::template addbind<G>(m,
                                                            grid_id,
                                                            (std::is_same<FVI, LVI>::value
#if HAVE_DUNE_FEM
                                                             && std::is_same<FVI, FPI>::value
                                                             && std::is_same<FVI, LPI>::value
#endif
                                                             )
                                                                ? ""
                                                                : "_leaf_view");
  for_Grid_and_Intersection<LVI, !(std::is_same<LVI, FVI>::value)>::template addbind<G>(m, grid_id, "_level_view");
#if HAVE_DUNE_FEM
  for_Grid_and_Intersection<FPI,
                            !(std::is_same<FPI, FVI>::value
                              || std::is_same<FPI, LVI>::value)>::template addbind<G>(m, grid_id, "_leaf_part");
  for_Grid_and_Intersection<LPI,
                            !(std::is_same<LPI, FVI>::value || std::is_same<LPI, LVI>::value
                              || std::is_same<LPI, FPI>::value)>::template addbind<G>(m, grid_id, "_level_part");
#endif // HAVE_DUNE_FEM


#define BIND(V, id)                                                                                                    \
  bind_Walker<V>(m, grid_id + "_" + id);                                                                               \
                                                                                                                       \
  m.def(std::string(std::string("make_all_dirichlet_boundary_info__") + id).c_str(),                                   \
        [](const GridProvider<typename extract_grid<V>::type>& /*grid_provider*/) {                                    \
          return AllDirichletBoundaryInfo<typename Intersection<V>::Type>();                                           \
        },                                                                                                             \
        "grid_provider"_a);                                                                                            \
  m.def(std::string(std::string("make_all_neumann_boundary_info__") + id).c_str(),                                     \
        [](const GridProvider<typename extract_grid<V>::type>& /*grid_provider*/) {                                    \
          return AllNeumannBoundaryInfo<typename Intersection<V>::Type>();                                             \
        },                                                                                                             \
        "grid_provider"_a);                                                                                            \
  m.def(std::string(std::string("make_normalbased_boundary_info__") + id).c_str(),                                     \
        [](const GridProvider<typename extract_grid<V>::type>& /*grid_provider*/, const Common::Configuration& cfg) {  \
          return *NormalBasedBoundaryInfo<typename Intersection<V>::Type>::create(cfg);                                \
        },                                                                                                             \
        "grid_provider"_a,                                                                                             \
        "cfg"_a = normalbased_boundaryinfo_default_config());                                                          \
  addbind_WhichIntersection<V>(m, grid_id, id); /*                                                           end BIND  \
                                                   */

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
#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
  addbind_for_Grid<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>(m, "2d_simplex_aluconform");
#endif
#if HAVE_UG
  addbind_for_Grid<Dune::UGGrid<2>>(m, "2d_simplex_uggrid");
#endif

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
