// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)
//   René Fritze     (2020)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/grid/functors/refinement.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/view/coupling.hh>

#include "interfaces.hh"


namespace Dune::XT::Grid::bindings {


template <class GV>
class MaximumEntityVolumeRefineFunctor
{
  using G = typename GV::Grid;
  static_assert(is_grid<G>::value);
  using I = extract_intersection_t<GV>;

public:
  using type = Grid::MaximumEntityVolumeRefineFunctor<GV>;
  using base_type = Grid::ElementFunctor<GV>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id = grid_name<G>::value(),
                         const std::string& layer_id = "",
                         const std::string& class_id = "maximum_element_volume_refine_functor")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto ClassId = Common::to_camel_case(class_id);
    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
    if (!layer_id.empty())
      ClassName += "_" + layer_id;
    const std::string doc{ClassId + "( " + grid_id + " variant)"};
    bound_type c(m, ClassName.c_str(), doc.c_str());
    c.def(py::init([](GridProvider<G>& grid_provider, const double& volume) {
            return std::make_unique<type>(grid_provider.grid(), volume, 1.);
          }),
          "grid_provider"_a,
          "volume"_a,
          py::keep_alive<0, 1>());
    c.def("__repr__", [ClassId](type&) { return ClassId + "(grid_provider=\?\?\?, volume=\?\?\?)"; });
    // there's no result member in the functor??
    // c.def_property_readonly("result", [](const type& self) { return self.result(); });

    m.def(
        ClassId.c_str(),
        [](GridProvider<G>& grid_provider, const double& volume) {
          return std::make_unique<type>(grid_provider.grid(), volume, 1.);
        },
        "grid_provider"_a,
        "volume"_a,
        py::keep_alive<0, 1>());

    return c;
  } // ... bind(...)

  static void bind_leaf_factory(pybind11::module& m,
                                const std::string& class_id = "maximum_element_volume_refine_functor")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;
    m.def(
        Common::to_camel_case(class_id).c_str(),
        [](GridProvider<G>& grid_provider, const double& volume) {
          return std::make_unique<type>(grid_provider.grid(), volume, 1.);
        },
        "grid_provider"_a,
        "volume"_a,
        py::keep_alive<0, 1>());
  } // ... bind_leaf_factory(...)

  static void bind_coupling_factory(pybind11::module& m,
                                    const std::string& class_id = "maximum_element_volume_refine_functor")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;
    m.def(
        Common::to_camel_case(class_id).c_str(),
        [](CouplingGridProvider<G>& coupling_grid_provider, const double& volume) {
          return std::make_unique<type>(coupling_grid_provider.grid(), volume, 1.);
        },
        "coupling_grid_provider"_a,
        "volume"_a,
        py::keep_alive<0, 1>());
  }
}; // class MaximumEntityVolumeRefineFunctor


} // namespace Dune::XT::Grid::bindings


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct MaximumEntityVolumeRefineFunctor_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using LGV = typename G::LeafGridView;
  static const size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::XT::Grid::bindings::grid_name;
    Dune::XT::Grid::bindings::MaximumEntityVolumeRefineFunctor<LGV>::bind(m, grid_name<G>::value(), "leaf");
    Dune::XT::Grid::bindings::MaximumEntityVolumeRefineFunctor<LGV>::bind_leaf_factory(m);
    if constexpr (d == 2) {
      using GridGlueType = Dune::XT::Grid::DD::Glued<G, G, Dune::XT::Grid::Layers::leaf>;
      using CGV = Dune::XT::Grid::CouplingGridView<GridGlueType>;
      Dune::XT::Grid::bindings::MaximumEntityVolumeRefineFunctor<CGV>::bind(m, grid_name<G>::value(), "coupling");
      Dune::XT::Grid::bindings::MaximumEntityVolumeRefineFunctor<CGV>::bind_coupling_factory(m);
    }
    MaximumEntityVolumeRefineFunctor_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct MaximumEntityVolumeRefineFunctor_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_grid_functors_refinement, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.grid._grid_boundaryinfo_interfaces");
  py::module::import("dune.xt.grid._grid_boundaryinfo_types");
  py::module::import("dune.xt.grid._grid_gridprovider_provider");
  py::module::import("dune.xt.grid._grid_functors_interfaces");

  MaximumEntityVolumeRefineFunctor_for_all_grids<>::bind(m);
}
