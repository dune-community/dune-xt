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
#include <dune/xt/grid/functors/boundary-detector.hh>

#include <python/dune/xt/grid/grids.bindings.hh>

#include "interfaces.hh"


namespace Dune::XT::Grid::bindings {


template <class G>
class BoundaryDetectorFunctor
{
  static_assert(is_grid<G>::value, "");
  using GV = typename G::LeafGridView;
  using I = extract_intersection_t<GV>;

public:
  using type = Grid::BoundaryDetectorFunctor<GV>;
  using base_type = Grid::IntersectionFunctor<GV>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "boundary_detector_functor",
                         const std::string& grid_id = grid_name<G>::value())
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto ClassId = Common::to_camel_case(class_id);
    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
    bound_type c(m, ClassName.c_str(), std::string(ClassId + "( " + grid_id + " variant)").c_str());
    c.def(py::init([](const BoundaryInfo<I>& boundary_info,
                      const BoundaryType& boundary_type,
                      const std::string& logging_prefix) {
            return std::make_unique<type>(boundary_info, boundary_type, logging_prefix);
          }),
          "boundary_info"_a,
          "boundary_type"_a,
          "logging_prefix"_a = "",
          py::keep_alive<0, 1>());
    c.def("__repr__", [ClassId](type&) { return ClassId + "(boundary_info=, boundary_type=)"; });
    c.def_property_readonly("result", [](const type& self) { return self.result(); });

    m.def(
        ClassId.c_str(),
        [](const GridProvider<G>&,
           const BoundaryInfo<I>& boundary_info,
           const BoundaryType& boundary_type,
           const std::string& logging_prefix) {
          return std::make_unique<type>(boundary_info, boundary_type, logging_prefix);
        },
        "grid_provider"_a,
        "boundary_info"_a,
        "boundary_type"_a,
        "logging_prefix"_a = "",
        py::keep_alive<0, 2>());

    return c;
  } // ... bind(...)
}; // class BoundaryDetectorFunctor


} // namespace Dune::XT::Grid::bindings


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct BoundaryDetectorFunctor_for_all_grids
{
  static void bind(pybind11::module& m)
  {
    Dune::XT::Grid::bindings::BoundaryDetectorFunctor<Dune::XT::Common::tuple_head_t<GridTypes>>::bind(m);
    BoundaryDetectorFunctor_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct BoundaryDetectorFunctor_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_grid_functors_boundary_detector, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.grid._grid_boundaryinfo_interfaces");
  py::module::import("dune.xt.grid._grid_boundaryinfo_types");
  py::module::import("dune.xt.grid._grid_gridprovider_provider");
  py::module::import("dune.xt.grid._grid_functors_interfaces");

  BoundaryDetectorFunctor_for_all_grids<>::bind(m);
}
