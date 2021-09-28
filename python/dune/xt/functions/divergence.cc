// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019 - 2020)
//   René Fritze     (2019 - 2020)
//   Tobias Leibner  (2019 - 2021)

#include "config.h"

#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/functional.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/base/derivatives-of-grid-functions.hh>

#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/fmatrix.hh>
#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/grid/traits.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune::XT::Functions::bindings {


template <class G, class E>
class DivergenceGridFunction
{
  using GP = XT::Grid::GridProvider<G>;
  static const size_t d = G::dimension;
  using GF = Functions::GridFunction<E, d>;

public:
  using type = Functions::DivergenceGridFunction<GF>;
  using base_type = Functions::GridFunctionInterface<E>;
  using bound_type = pybind11::class_<type, base_type>;


  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id = Grid::bindings::grid_name<G>::value(),
                         const std::string& layer_id = "",
                         const std::string& class_id = "divergence_grid_function")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    class_name += "_to_1d";
    const auto ClassName = Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), Common::to_camel_case(class_id).c_str());
    c.def(py::init<GF, const std::string&>(), "grid_function"_a, "name"_a = "DivergenceGridFunction");

    m.def(
        "divergence",
        [](const GridFunctionInterface<E, d>& grid_function, const std::string& name) {
          return new type(grid_function, name);
        },
        "grid_function"_a,
        "name"_a = "DivergenceGridFunction");
    m.def(
        "divergence",
        [](GF grid_function, const std::string& name) { return new type(grid_function, name); },
        "grid_function"_a,
        "name"_a = "DivergenceGridFunction");

    return c;
  }
}; // class GridFunction


} // namespace Dune::XT::Functions::bindings


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct DivergenceGridFunction_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;
  using E = Dune::XT::Grid::extract_entity_t<GV>;

  static void bind(pybind11::module& m)
  {
    using Dune::XT::Functions::bindings::DivergenceGridFunction;
    using Dune::XT::Grid::bindings::grid_name;

    DivergenceGridFunction<G, E>::bind(m, grid_name<G>::value());

    DivergenceGridFunction_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct DivergenceGridFunction_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_functions_divergence, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_1d");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_2d");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_3d");
  py::module::import("dune.xt.functions._functions_gridfunction");

  DivergenceGridFunction_for_all_grids<>::bind(m);
} // PYBIND11_MODULE(...)
