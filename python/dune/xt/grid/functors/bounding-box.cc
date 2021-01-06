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
#include <dune/xt/grid/functors/bounding-box.hh>
#include <dune/xt/grid/view/coupling.hh>
#include <python/dune/xt/common/fvector.hh>

#include "interfaces.hh"


namespace Dune::XT::Grid::bindings {


template <class GV>
class MinMaxCoordinateFunctor
{
  using G = typename GV::Grid;
  static_assert(is_grid<G>::value);
  using I = extract_intersection_t<GV>;

public:
  using type = Grid::MinMaxCoordinateFunctor<GV>;
  using base_type = Grid::ElementFunctor<GV>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id = grid_name<G>::value(),
                         const std::string& layer_id = "",
                         const std::string& class_id = "bounding_box_functor",
                         )
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto ClassId = Common::to_camel_case(class_id);
    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
    if (!layer_id.empty())
      ClassName += "_" + layer_id;
    bound_type c(m, ClassName.c_str(), std::string(ClassId + "( " + grid_id + " variant)").c_str());
    c.def(py::init());
    c.def("__repr__", [ClassId](type&) { return ClassId + "()"; });
    c.def_property_readonly("result", [](const type& self) { return self.result(); });

    return c;
  } // ... bind(...)

  static bound_type bind_leaf_factory(pybind11::module& m,
                                      const std::string& class_id = "bounding_box_functor")
  {
      using namespace pybind11::literals;
      // factories
      m.def(
          Common::to_camel_case(class_id).c_str(),
          [](GridProvider<G>&) { return new type(); }, "grid_provider"_a);
  }

  static bound_type bind_coupling_factory(pybind11::module& m,
                                          const std::string& class_id = "bounding_box_functor")
  {
      using namespace pybind11::literals;
      // factories
      m.def(
          Common::to_camel_case(class_id).c_str(),
          [](CouplingGridProvider<GV>&) { return new type(); }, "coupling_grid_provider"_a);
  }

}; // class MinMaxCoordinateFunctor


} // namespace Dune::XT::Grid::bindings


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct MinMaxCoordinateFunctor_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;

  static void bind(pybind11::module& m)
  {
    using Dune::XT::Grid::bindings::grid_name;
    Dune::XT::Grid::bindings::MinMaxCoordinateFunctor<GV>::bind(m, grid_name<G>::value(), "leaf");
//    Dune::XT::Grid::bindings::MinMaxCoordinateFunctor<GV>::bind_leaf_factory(m);
    MinMaxCoordinateFunctor_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct MinMaxCoordinateFunctor_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


template <class GridTypes = Dune::XT::Grid::Available2dGridTypes>
struct MinMaxCoordinateFunctor_for_all_coupling_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GridGlueType = Dune::XT::Grid::DD::Glued<G,G,Dune::XT::Grid::Layers::leaf>;
  using CGV = Dune::XT::Grid::CouplingGridView<GridGlueType>;

  static void bind(pybind11::module& m)
  {
    using Dune::XT::Grid::bindings::grid_name;
    Dune::XT::Grid::bindings::MinMaxCoordinateFunctor<CGV>::bind(m, grid_name<G>::value(), "coupling");
//    Dune::XT::Grid::bindings::MinMaxCoordinateFunctor<CGV>::bind_coupling_factory(m);
    MinMaxCoordinateFunctor_for_all_coupling_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct MinMaxCoordinateFunctor_for_all_coupling_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_grid_functors_bounding_box, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.grid._grid_gridprovider_provider");
  py::module::import("dune.xt.grid._grid_functors_interfaces");

  MinMaxCoordinateFunctor_for_all_grids<>::bind(m);
  MinMaxCoordinateFunctor_for_all_coupling_grids<>::bind(m);
}
