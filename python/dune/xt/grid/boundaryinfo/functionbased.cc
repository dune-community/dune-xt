// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2021)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/boundaryinfo/functionbased.hh>
#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/coupling.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/view/coupling.hh>

#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

using namespace Dune::XT;
using namespace Dune::XT::Grid::bindings;

namespace Dune::XT::Grid::bindings {

template <class GV>
struct FunctionBasedBoundaryInfo
{
  using G = Dune::XT::Grid::extract_grid_t<GV>;
  using I = Dune::XT::Grid::extract_intersection_t<GV>;

public:
  using type = Dune::XT::Grid::FunctionBasedBoundaryInfo<I>;
  using base_type = Dune::XT::Grid::BoundaryInfo<I>;
  using bound_type = pybind11::class_<type, base_type>;
  using D = typename type::DomainFieldType;

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id = grid_name<G>::value(),
                         const std::string& layer_id = "",
                         const std::string& class_id = "function_based_boundary_info")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    const auto ClassName = XT::Common::to_camel_case(class_name);

    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(
        py::init([](const Grid::BoundaryType& default_boundary_type, const D& tol, const std::string& logging_prefix) {
          return std::make_unique<type>(default_boundary_type, tol, logging_prefix);
        }),
        "default_boundary_type"_a,
        "tolerance"_a = 1e-10,
        "logging_prefix"_a = "");
    c.def(
        "register_new_function",
        [](type& self, const typename type::FunctionType& function, const Grid::BoundaryType& boundary_type) {
          self.register_new_function(function, boundary_type);
        },
        "function"_a,
        "boundary_type"_a);
    return c;
  }

  static void bind_leaf_factory(pybind11::module& m, const std::string& class_id = "function_based_boundary_info")
  {
    using namespace pybind11::literals;
    m.def(
        Common::to_camel_case(class_id).c_str(),
        [](const Grid::GridProvider<G>&,
           const Grid::BoundaryType& default_boundary_type,
           const D& tol,
           const std::string& logging_prefix) {
          return std::make_unique<type>(default_boundary_type, tol, logging_prefix);
        },
        "grid_provider"_a,
        "default_boundary_type"_a,
        "tolerance"_a = 1e-10,
        "logging_prefix"_a = "");
  } // ... bind_leaf_factory(...)

  static void bind_coupling_factory(pybind11::module& m, const std::string& class_id = "function_based_boundary_info")
  {
    using namespace pybind11::literals;
    m.def(
        Common::to_camel_case(class_id).c_str(),
        [](const CouplingGridProvider<GV>&,
           const Grid::BoundaryType& default_boundary_type,
           const D& tol,
           const std::string& logging_prefix) {
          return std::make_unique<type>(default_boundary_type, tol, logging_prefix);
        },
        "grid_provider"_a,
        "default_boundary_type"_a,
        "tolerance"_a = 1e-10,
        "logging_prefix"_a = "");
  } // ... bind_coupling_factory(...)


}; // struct FunctionBasedBoundaryInfo

} // namespace Dune::XT::Grid::bindings


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct FunctionBasedBoundaryInfo_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using LGV = typename G::LeafGridView;
  static const size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::XT::Grid::bindings::grid_name;
    Dune::XT::Grid::bindings::FunctionBasedBoundaryInfo<LGV>::bind(m, grid_name<G>::value(), "leaf");
    Dune::XT::Grid::bindings::FunctionBasedBoundaryInfo<LGV>::bind_leaf_factory(m);
#if HAVE_DUNE_GRID_GLUE
    if constexpr (d < 3) {
      using GridGlueType = Dune::XT::Grid::DD::Glued<G, G, Dune::XT::Grid::Layers::leaf>;
      using CGV = Dune::XT::Grid::CouplingGridView<GridGlueType>;
      Dune::XT::Grid::bindings::FunctionBasedBoundaryInfo<CGV>::bind(m, grid_name<G>::value(), "coupling");
      Dune::XT::Grid::bindings::FunctionBasedBoundaryInfo<CGV>::bind_coupling_factory(m);
    }
#endif
    FunctionBasedBoundaryInfo_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct FunctionBasedBoundaryInfo_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_grid_boundaryinfo_functionbased, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.grid._grid_gridprovider_provider");
  py::module::import("dune.xt.grid._grid_boundaryinfo_interfaces");
  py::module::import("dune.xt.grid._grid_boundaryinfo_types");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_1d");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_2d");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_3d");

  FunctionBasedBoundaryInfo_for_all_grids<>::bind(m);
}
