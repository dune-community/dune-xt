#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <string>
#include <vector>

#include <dune/pybindxi/pybind11.h>
//#include <dune/pybindxi/stl_bind.h> // <- see dune/xt/common/bindings.cc

#include <dune/xt/common/configuration.pbh>
#include <dune/xt/common/fvector.pbh>

#include "grids.hh"

#include "boundaryinfo.pbh"
#include "gridprovider.pbh"
#include "walker.pbh"

namespace py = pybind11;
using namespace pybind11::literals;


PYBIND11_PLUGIN(grid)
{
  py::module m("grid", "dune-xt-grid");

  py::module::import("common");

  typedef Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>> YaspGrid2dType;
  const std::string yasp_grid_2d_id = "2d_cube_yaspgrid";

  Dune::XT::Grid::bind_GridProvider<YaspGrid2dType>(m, yasp_grid_2d_id);
  Dune::XT::Grid::bind_make_cube_grid<YaspGrid2dType>(m, yasp_grid_2d_id);

  Dune::XT::Grid::bind_BoundaryTypeInterface(m);
  Dune::XT::Grid::bind_BoundaryType<Dune::XT::Grid::NoBoundary>(m, "NoBoundary");
  Dune::XT::Grid::bind_BoundaryType<Dune::XT::Grid::UnknownBoundary>(m, "UnknownBoundary");
  Dune::XT::Grid::bind_BoundaryType<Dune::XT::Grid::DirichletBoundary>(m, "DirichletBoundary");
  Dune::XT::Grid::bind_BoundaryType<Dune::XT::Grid::NeumannBoundary>(m, "NeumannBoundary");
  Dune::XT::Grid::bind_BoundaryType<Dune::XT::Grid::RobinBoundary>(m, "RobinBoundary");

  Dune::XT::Grid::bind_Walker<typename YaspGrid2dType::LevelGridView>(m, yasp_grid_2d_id + "_level_view");
  Dune::XT::Grid::bind_Walker<typename YaspGrid2dType::LeafGridView>(m, yasp_grid_2d_id + "_leaf_view");

  return m.ptr();
}

#endif // HAVE_DUNE_PYBINDXI
