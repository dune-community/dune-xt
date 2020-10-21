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
#include <dune/xt/grid/functors/refinement.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/type_traits.hh>

#include "interfaces.hh"


namespace Dune::XT::Grid::bindings {


template <class G>
class MaximumEntityVolumeRefineFunctor
{
  static_assert(is_grid<G>::value);
  using GV = typename G::LeafGridView;
  using I = extract_intersection_t<GV>;

public:
  using type = Grid::MaximumEntityVolumeRefineFunctor<GV>;
  using base_type = Grid::ElementFunctor<GV>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "maximum_element_volume_refine_functor",
                         const std::string& grid_id = grid_name<G>::value())
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto ClassId = Common::to_camel_case(class_id);
    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
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
}; // class MaximumEntityVolumeRefineFunctor


} // namespace Dune::XT::Grid::bindings


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct MaximumEntityVolumeRefineFunctor_for_all_grids
{
  static void bind(pybind11::module& m)
  {
    Dune::XT::Grid::bindings::MaximumEntityVolumeRefineFunctor<Dune::XT::Common::tuple_head_t<GridTypes>>::bind(m);
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
