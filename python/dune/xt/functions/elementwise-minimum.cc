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
#include <dune/xt/functions/elementwise-minimum.hh>

#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/fmatrix.hh>
#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/grid/traits.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune::XT::Functions::bindings {


template <class G, class E, size_t r = 1, size_t rC = 1>
class ElementwiseMinimumFunction
{
  using GP = XT::Grid::GridProvider<G>;
  using GF = Functions::GridFunction<E, r, rC>;
  static const size_t d = G::dimension;

public:
  using type = Functions::ElementwiseMinimumFunction<Functions::GridFunctionInterface<E, r, rC>>;
  using base_type = Functions::GridFunctionInterface<E>;
  using bound_type = pybind11::class_<type, base_type>;


  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id = Grid::bindings::grid_name<G>::value(),
                         const std::string& layer_id = "",
                         const std::string& class_id = "elementwise_minimum_function")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    class_name += "_to_" + Common::to_string(r);
    if (rC > 1)
      class_name += "x" + Common::to_string(rC);
    class_name += "d";
    const auto ClassName = Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), Common::to_camel_case(class_id).c_str());
    c.def(py::init<GF, const int, const std::string&>(),
          "grid_function"_a,
          "search_quadrature_order"_a,
          "name"_a = "ElementwiseMinimumFunction");

    const auto FactoryName = Common::to_camel_case(class_id);
    m.def(
        FactoryName.c_str(),
        [](GF grid_function, const int search_quadrature_order, const std::string& name) {
          return new type(grid_function, search_quadrature_order, name);
        },
        "grid_function"_a,
        "search_quadrature_order"_a,
        "name"_a = "ElementwiseMinimumFunction");
    m.def(
        FactoryName.c_str(),
        [](const Functions::GridFunction<E, r, rC>& grid_function,
           const int search_quadrature_order,
           const std::string& name) { return new type(grid_function, search_quadrature_order, name); },
        "grid_function"_a,
        "search_quadrature_order"_a,
        "name"_a = "ElementwiseMinimumFunction");

    return c;
  }
}; // class GridFunction


} // namespace Dune::XT::Functions::bindings


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct ElementwiseMinimumFunction_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;
  using E = Dune::XT::Grid::extract_entity_t<GV>;

  static void bind(pybind11::module& m)
  {
    using Dune::XT::Functions::bindings::ElementwiseMinimumFunction;
    using Dune::XT::Grid::bindings::grid_name;

    ElementwiseMinimumFunction<G, E, 1, 1>::bind(m, grid_name<G>::value());
    ElementwiseMinimumFunction<G, E, 2, 2>::bind(m, grid_name<G>::value());
    ElementwiseMinimumFunction<G, E, 3, 3>::bind(m, grid_name<G>::value());

    ElementwiseMinimumFunction_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct ElementwiseMinimumFunction_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_functions_elementwise_minimum, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_1d");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_2d");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_3d");
  py::module::import("dune.xt.functions._functions_gridfunction");

  ElementwiseMinimumFunction_for_all_grids<>::bind(m);
} // PYBIND11_MODULE(...)
