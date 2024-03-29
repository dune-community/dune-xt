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
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/checkerboard.hh>

#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/fmatrix.hh>
#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/grid/traits.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune::XT::Functions::bindings {


template <class G, class E, size_t r = 1, size_t rC = 1, class R = double>
class CheckerboardFunction
{
  using GP = XT::Grid::GridProvider<G>;
  static constexpr size_t d = G::dimension;

public:
  using type = Functions::CheckerboardFunction<E, r, rC, R>;
  using base_type = Functions::GridFunctionInterface<E, r, rC, R>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id = Grid::bindings::grid_name<G>::value(),
                         const std::string& layer_id = "",
                         const std::string& class_id = "checkerboard_function")
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
    if (!std::is_same<R, double>::value)
      class_name += "_" + Common::Typename<R>::value(/*fail_wo_typeid=*/true);
    const auto ClassName = Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), Common::to_camel_case(class_id).c_str());
    c.def(py::init<const typename type::DomainType&,
                   const typename type::DomainType&,
                   const FieldVector<size_t, d>&,
                   const std::vector<typename type::RangeType>&,
                   const std::string>(),
          "lower_left"_a,
          "upper_right"_a,
          "num_elements"_a,
          "values"_a,
          "name"_a = type::static_id());

    m.def(
        Common::to_camel_case(class_id).c_str(),
        [](const GP& /*grid*/,
           const typename type::DomainType& lower_left,
           const typename type::DomainType& upper_right,
           const FieldVector<size_t, d>& num_elements,
           const std::vector<typename type::RangeType>& values,
           const std::string& name = "checkerboard") {
          return new type(lower_left, upper_right, num_elements, values, name);
        },
        "grid"_a,
        "lower_left"_a,
        "upper_right"_a,
        "num_elements"_a,
        "values"_a,
        "name"_a = type::static_id());
    return c;
  }
}; // class CheckerboardFunction


} // namespace Dune::XT::Functions::bindings


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct CheckerboardFunction_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;
  using E = Dune::XT::Grid::extract_entity_t<GV>;
  static constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::XT::Functions::bindings::CheckerboardFunction;
    using Dune::XT::Grid::bindings::grid_name;

    CheckerboardFunction<G, E, 1, 1>::bind(m);
    CheckerboardFunction<G, E, 1, 2>::bind(m);
    CheckerboardFunction<G, E, 1, 3>::bind(m);
    CheckerboardFunction<G, E, 2, 1>::bind(m);
    CheckerboardFunction<G, E, 2, 2>::bind(m);
    CheckerboardFunction<G, E, 2, 3>::bind(m);
    CheckerboardFunction<G, E, 3, 1>::bind(m);
    CheckerboardFunction<G, E, 3, 2>::bind(m);
    CheckerboardFunction<G, E, 3, 3>::bind(m);

    CheckerboardFunction_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct CheckerboardFunction_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_functions_checkerboard, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_1d");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_2d");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_3d");

  CheckerboardFunction_for_all_grids<>::bind(m);
} // PYBIND11_MODULE(...)
