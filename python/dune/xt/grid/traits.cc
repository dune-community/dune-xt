// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)
//   Tim Keil        (2021)
//   Tobias Leibner  (2021)

#include "config.h"

#include "traits.hh"

//#include <dune/xt/grid/type_traits.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct CouplingIntersection_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using LGV = typename G::LeafGridView;
  static const size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
#if HAVE_DUNE_GRID_GLUE
    if constexpr (d < 3) {
      using Dune::XT::Grid::bindings::grid_name;
      Dune::XT::Grid::bindings::CouplingIntersectionBinder<G, G>::bind(m, grid_name<G>::value());
    }
#endif
    CouplingIntersection_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct CouplingIntersection_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_grid_traits, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT::Grid::bindings;

  py::module::import("dune.xt.common");

  pybind11::class_<Simplex>(m, "Simplex", "tag: Simplex").def(py::init()).def("__repr__", [](const Simplex&) {
    return "Simplex()";
  });
  pybind11::class_<Cube>(m, "Cube", "tag: Cube").def(py::init()).def("__repr__", [](const Cube&) { return "Cube()"; });
  pybind11::class_<Pyramid>(m, "Pyramid", "tag: Pyramid").def(py::init()).def("__repr__", [](const Pyramid&) {
    return "Pyramid()";
  });
  pybind11::class_<Prism>(m, "Prism", "tag: Prism").def(py::init()).def("__repr__", [](const Prism&) {
    return "Prism()";
  });

  pybind11::class_<Dimension<0>>(m, "Dimension0", "tag: Dimension0")
      .def(py::init())
      .def("__repr__", [](const Dimension<0>&) { return "Dim(0)"; });
  pybind11::class_<Dimension<1>>(m, "Dimension1", "tag: Dimension1")
      .def(py::init())
      .def("__repr__", [](const Dimension<1>&) { return "Dim(1)"; });
  pybind11::class_<Dimension<2>>(m, "Dimension2", "tag: Dimension2")
      .def(py::init())
      .def("__repr__", [](const Dimension<2>&) { return "Dim(2)"; });
  pybind11::class_<Dimension<3>>(m, "Dimension3", "tag: Dimension3")
      .def(py::init())
      .def("__repr__", [](const Dimension<3>&) { return "Dim(3)"; });
  pybind11::class_<Dimension<4>>(m, "Dimension4", "tag: Dimension4")
      .def(py::init())
      .def("__repr__", [](const Dimension<4>&) { return "Dim(4)"; });
  pybind11::class_<Dimension<5>>(m, "Dimension5", "tag: Dimension5")
      .def(py::init())
      .def("__repr__", [](const Dimension<5>&) { return "Dim(5)"; });
  pybind11::class_<Dimension<6>>(m, "Dimension6", "tag: Dimension6")
      .def(py::init())
      .def("__repr__", [](const Dimension<6>&) { return "Dim(6)"; });
  pybind11::class_<Dimension<7>>(m, "Dimension7", "tag: Dimension7")
      .def(py::init())
      .def("__repr__", [](const Dimension<7>&) { return "Dim(7)"; });
  pybind11::class_<Dimension<8>>(m, "Dimension8", "tag: Dimension8")
      .def(py::init())
      .def("__repr__", [](const Dimension<8>&) { return "Dim(8)"; });
  pybind11::class_<Dimension<9>>(m, "Dimension9", "tag: Dimension9")
      .def(py::init())
      .def("__repr__", [](const Dimension<9>&) { return "Dim(9)"; });

  pybind11::class_<LeafIntersection>(m, "LeafIntersection", "tag: LeafIntersection")
      .def(py::init())
      .def("__repr__", [](const LeafIntersection&) { return "LeafIntersection"; });

  CouplingIntersection_for_all_grids<>::bind(m);
}
