// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)
//   René Fritze     (2020)
//   Tim Keil        (2021)
//   Tobias Leibner  (2021)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/boundaryinfo/allreflecting.hh>
#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/coupling.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/view/coupling.hh>

#include <python/dune/xt/grid/grids.bindings.hh>

using namespace Dune::XT;
using namespace Dune::XT::Grid::bindings;

namespace Dune::XT::Grid::bindings {

template <class GV>
class AllReflectingBoundaryInfo
{
  using G = Dune::XT::Grid::extract_grid_t<GV>;
  using I = Dune::XT::Grid::extract_intersection_t<GV>;

public:
  using type = Dune::XT::Grid::AllReflectingBoundaryInfo<I>;
  using base_type = Dune::XT::Grid::BoundaryInfo<I>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id = grid_name<G>::value(),
                         const std::string& layer_id = "",
                         const std::string& class_id = "all_reflecting_boundary_info")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init([]() { return std::make_unique<type>(); }));
    return c;
  }

  static void bind_leaf_factory(pybind11::module& m, const std::string& class_id = "all_reflecting_boundary_info")
  {
    using namespace pybind11::literals;
    m.def(
        Common::to_camel_case(class_id).c_str(), [](Grid::GridProvider<G>&) { return new type(); }, "grid_provider"_a);
  } // ... bind_leaf_factory(...)

  static void bind_coupling_factory(pybind11::module& m, const std::string& class_id = "all_reflecting_boundary_info")
  {
    using namespace pybind11::literals;
    m.def(
        Common::to_camel_case(class_id).c_str(),
        [](CouplingGridProvider<GV>&) { return new type(); },
        "coupling_grid_provider"_a);
  } // ... bind_coupling_factory(...)
}; // AllReflectingBoundaryInfo

} // namespace Dune::XT::Grid::bindings

template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct AllReflectingBoundaryInfo_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using LGV = typename G::LeafGridView;
  static const size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::XT::Grid::bindings::grid_name;
    Dune::XT::Grid::bindings::AllReflectingBoundaryInfo<LGV>::bind(m, grid_name<G>::value(), "leaf");
    Dune::XT::Grid::bindings::AllReflectingBoundaryInfo<LGV>::bind_leaf_factory(m);
#if HAVE_DUNE_GRID_GLUE
    if constexpr (d < 3) {
      using GridGlueType = Dune::XT::Grid::DD::Glued<G, G, Dune::XT::Grid::Layers::leaf>;
      using CGV = Dune::XT::Grid::CouplingGridView<GridGlueType>;
      Dune::XT::Grid::bindings::AllReflectingBoundaryInfo<CGV>::bind(m, grid_name<G>::value(), "coupling");
      Dune::XT::Grid::bindings::AllReflectingBoundaryInfo<CGV>::bind_coupling_factory(m);
    }
#endif
    AllReflectingBoundaryInfo_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct AllReflectingBoundaryInfo_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_grid_boundaryinfo_allreflecting, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.grid._grid_gridprovider_provider");
  py::module::import("dune.xt.grid._grid_boundaryinfo_interfaces");

  AllReflectingBoundaryInfo_for_all_grids<>::bind(m);
}
