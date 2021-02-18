// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)
//   René Fritze     (2020)
//   Tobias Leibner  (2020)
//
// This file is part of the dune-gdt project:

#include "config.h"

#include <sstream>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/operators.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/grid/print.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/view/coupling.hh>

#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>


namespace Dune::XT::Grid::bindings {


template <class GV>
class Intersection
{
  using G = XT::Grid::extract_grid_t<GV>;
  using I = XT::Grid::extract_intersection_t<GV>;
  using D = typename G::ctype;
  static constexpr size_t d = G::dimension;

public:
  using type = I;
  using bound_type = pybind11::class_<I>;

  using LocalCoordinateType = FieldVector<D, d - 1>;
  using GlobalCoordinateType = FieldVector<D, d>;

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id = grid_name<G>::value(),
                         const std::string& layer_id = "",
                         const std::string& class_id = "intersection")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());

    // properties
    c.def_property_readonly("boundary", [](type& self) { return self.boundary(); });
    c.def_property_readonly("boundary_segment_index", [](type& self) { return self.boundarySegmentIndex(); });
    c.def_property_readonly("neighbor", [](type& self) { return self.neighbor(); });
    c.def_property_readonly("conforming", [](type& self) { return self.conforming(); });
    c.def_property_readonly("index_in_inside", [](type& self) { return self.indexInInside(); });
    c.def_property_readonly("index_in_outside", [](type& self) { return self.indexInOutside(); });
    c.def_property_readonly("center_unit_outer_normal", [](type& self) { return self.centerUnitOuterNormal(); });

    // special methods/operators
    c.def("__str__", [](type& self) {
      std::stringstream ss;
      ss << print(self);
      return ss.str();
    });
    c.def("__repr__", [](type& self) {
      std::stringstream ss;
      ss << repr(self);
      return ss.str();
    });
    c.def(
        "__eq__", [](type& self, const type& other) { return self == other; }, py::is_operator());
    c.def(
        "__neq__", [](type& self, const type& other) { return self != other; }, py::is_operator());

    // methods
    c.def(
        "to_global",
        [](type& self, const LocalCoordinateType& x) { return self.geometry().global(x); },
        "point_in_reference_intersection"_a);
    c.def(
        "to_local",
        [](type& self, const GlobalCoordinateType& x) { return self.geometry().local(x); },
        "point_in_physical_space"_a);
    c.def(
        "outer_normal",
        [](type& self, const LocalCoordinateType& x) { return self.outerNormal(x); },
        "point_in_reference_intersection"_a);
    c.def(
        "unit_outer_normal",
        [](type& self, const LocalCoordinateType& x) { return self.unitOuterNormal(x); },
        "point_in_reference_intersection"_a);

    return c;
  } // ... bind(...)
}; // class Intersection


} // namespace Dune::XT::Grid::bindings


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct Intersection_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using LGV = typename G::LeafGridView;
  static const size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::XT::Grid::bindings::grid_name;
    Dune::XT::Grid::bindings::Intersection<LGV>::bind(m, grid_name<G>::value(), "leaf");
#if HAVE_DUNE_GRID_GLUE
    if constexpr (d == 2) {
      using GridGlueType = Dune::XT::Grid::DD::Glued<G, G, Dune::XT::Grid::Layers::leaf>;
      using CGV = Dune::XT::Grid::CouplingGridView<GridGlueType>;
      Dune::XT::Grid::bindings::Intersection<CGV>::bind(m, grid_name<G>::value(), "coupling");
    }
#endif
    Intersection_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct Intersection_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_grid_intersection, m)
{
  namespace py = pybind11;
  using namespace Dune;

  py::module::import("dune.xt.common");

  Intersection_for_all_grids<>::bind(m);
}
