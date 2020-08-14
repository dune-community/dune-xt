// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

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
#include <dune/xt/functions/elementwise-diameter.hh>

#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/fmatrix.hh>
#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/grid/traits.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace bindings {


template <class G, class E>
class ElementwiseDiameterFunction
{
  using GP = XT::Grid::GridProvider<G>;
  static const size_t d = G::dimension;

public:
  using type = Functions::ElementwiseDiameterFunction<E>;
  using base_type = Functions::GridFunctionInterface<E>;
  using bound_type = pybind11::class_<type, base_type>;

public:
  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id = Grid::bindings::grid_name<G>::value(),
                         const std::string& layer_id = "",
                         const std::string& class_id = "elementwise_diameter_function")
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
    c.def(py::init<const std::string&>(), "name"_a = "ElementwiseDiameterFunction");

    const auto FactoryName = Common::to_camel_case(class_id);
    m.def(
        FactoryName.c_str(),
        [](const GP&, const std::string& name) { return new type(name); },
        "grid"_a,
        "name"_a = "ElementwiseDiameterFunction",
        py::keep_alive<0, 1>());

    return c;
  }
}; // class GridFunction


} // namespace bindings
} // namespace Functions
} // namespace XT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct ElementwiseDiameterFunction_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  using E = Dune::XT::Grid::extract_entity_t<GV>;

  static void bind(pybind11::module& m)
  {
    using Dune::XT::Functions::bindings::ElementwiseDiameterFunction;
    using Dune::XT::Grid::bindings::grid_name;

    ElementwiseDiameterFunction<G, E>::bind(m, grid_name<G>::value());

    ElementwiseDiameterFunction_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct ElementwiseDiameterFunction_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_functions_elementwise_diameter, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_1d");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_2d");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_3d");

  ElementwiseDiameterFunction_for_all_grids<>::bind(m);
} // PYBIND11_MODULE(...)
