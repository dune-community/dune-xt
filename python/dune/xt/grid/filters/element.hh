// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)
//   René Fritze     (2020)

#ifndef PYTHON_DUNE_XT_GRID_FILTERS_ELEMENT_HH
#define PYTHON_DUNE_XT_GRID_FILTERS_ELEMENT_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/filters/element.hh>
#include <dune/xt/grid/gridprovider/coupling.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace bindings {


template <template <class> class Filter, class GV>
class InitlessElementFilter
{
  using G = typename GV::Grid;
  static_assert(is_grid<G>::value);

public:
  using type = Filter<GV>;
  using base_type = Grid::ElementFilter<GV>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id,
                         const std::string& layer_id = "",
                         const std::string& grid_id = grid_name<G>::value())
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto ClassId = Common::to_camel_case(class_id);
    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
    if (!layer_id.empty())
      ClassName += "_" + layer_id;
    ClassId += "_" + layer_id;
    bound_type c(m, ClassName.c_str(), std::string(ClassId + "( " + grid_id + " variant)").c_str());
    c.def(py::init([]() { return std::make_unique<type>(); }));
    c.def("__repr__", [ClassId](type&) { return ClassId + "()"; });

    return c;
  } // ... bind(...)

  static void bind_leaf_factory(pybind11::module& m, const std::string& class_id)
  {
    m.def(Common::to_camel_case(class_id).c_str(), [](const Grid::GridProvider<G>&) { return new type(); });
  } // ... bind_leaf_factory(...)

  static void bind_coupling_factory(pybind11::module& m, const std::string& class_id = "Walker")
  {
    m.def(Common::to_camel_case(class_id).c_str(), [](const CouplingGridProvider<GV>&) { return new type(); });
  }

}; // class InitlessElementFilter


} // namespace bindings
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // PYTHON_DUNE_XT_GRID_FILTERS_ELEMENT_HH
