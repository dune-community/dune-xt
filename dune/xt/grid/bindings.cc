#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <string>
#include <vector>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

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


template <class I, bool required = true>
struct for_Grid_and_Intersection
{
  template <class G>
  static void addbind(py::module& m, const std::string& grid_id, const std::string& id)
  {
    using namespace Dune::XT::Grid;
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
