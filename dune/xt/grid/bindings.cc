#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <string>
#include <vector>

#include <dune/pybindxi/pybind11.h>
//#include <dune/pybindxi/stl_bind.h> // <- see dune/xt/common/bindings.cc

#include <dune/xt/common/configuration.pbh>
#include <dune/xt/common/fvector.pbh>

#include "grids.hh"
#include "layers.hh"
#include "intersection.hh"
#include "type_traits.hh"

#include "boundaryinfo.pbh"
#include "gridprovider.pbh"
#include "walker.pbh"

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

#define BIND(V, id)                                                                                                    \
  bind_Walker<V>(m, grid_id + "_" + id);                                                                               \
  bind_BoundaryInfo<AllDirichletBoundaryInfo<typename Intersection<V>::Type>>(                                         \
      m, "AllDirichletBoundaryInfo__" + grid_id + "_" + id);                                                           \
  bind_BoundaryInfo<AllNeumannBoundaryInfo<typename Intersection<V>::Type>>(                                           \
      m, "AllNeumannBoundaryInfo__" + grid_id + "_" + id);                                                             \
  bind_BoundaryInfo<NormalBasedBoundaryInfo<typename Intersection<V>::Type>>(                                          \
      m, "NormalBasedBoundaryInfo__" + grid_id + "_" + id);                                                            \
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
  m.def(std::string(std::string("make_all_normalbased_boundary_info__") + id).c_str(),                                 \
        [](const GridProvider<typename extract_grid<V>::type>& /*grid_provider*/, const Common::Configuration& cfg) {  \
          return *NormalBasedBoundaryInfo<typename Intersection<V>::Type>::create(cfg);                                \
        },                                                                                                             \
        "grid_provider"_a,                                                                                             \
        "cfg"_a = normalbased_boundaryinfo_default_config()) /* end BIND */

  BIND(LevelView, "level_view");
  BIND(LeafView, "leaf_view");
#if HAVE_DUNE_FEM
  BIND(LevelPart, "level_part");
  BIND(LeafPart, "leaf_part");
#endif
#undef BIND
} // ... addbind_for_Grid(...)


PYBIND11_PLUGIN(grid)
{
  py::module m("grid", "dune-xt-grid");

  py::module::import("common");

  addbind_for_Grid<Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>>(m, "2d_cube_yaspgrid");
  addbind_for_Grid<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>(m, "2d_simplex_aluconform");

  return m.ptr();
}

#endif // HAVE_DUNE_PYBINDXI
