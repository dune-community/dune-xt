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
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/boundaryinfo/allreflecting.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/type_traits.hh>

#include <python/dune/xt/grid/grids.bindings.hh>

using namespace Dune::XT;
using namespace Dune::XT::Grid::bindings;


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct AllReflectingBoundaryInfo_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  using I = Dune::XT::Grid::extract_intersection_t<GV>;
  using type = Dune::XT::Grid::AllReflectingBoundaryInfo<I>;
  using base_type = Dune::XT::Grid::BoundaryInfo<I>;

  static void bind(pybind11::module& m,
                   const std::string& class_id = "all_reflecting_boundary_info",
                   const std::string& grid_id = grid_name<G>::value())
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
    pybind11::class_<type, base_type> c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init([]() { return std::make_unique<type>(); }));

    m.def(Common::to_camel_case(class_id).c_str(),
          [](const Grid::GridProvider<G>&) { return type(); },
          "grid_provider"_a);

    AllReflectingBoundaryInfo_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct AllReflectingBoundaryInfo_for_all_grids<boost::tuples::null_type>
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
