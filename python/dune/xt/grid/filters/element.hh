// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_XT_GRID_FILTERS_ELEMENT_HH
#define PYTHON_DUNE_XT_GRID_FILTERS_ELEMENT_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/filters/element.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace bindings {


template <template <class> class Filter, class G>
class InitlessElementFilter
{
  static_assert(is_grid<G>::value);
  using GV = typename G::LeafGridView;

public:
  using type = Filter<GV>;
  using base_type = Grid::ElementFilter<GV>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type
  bind(pybind11::module& m, const std::string& class_id, const std::string& grid_id = grid_name<G>::value())
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto ClassId = Common::to_camel_case(class_id);
    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
    bound_type c(m, ClassName.c_str(), std::string(ClassId + "( " + grid_id + " variant)").c_str());
    c.def(py::init([]() { return std::make_unique<type>(); }));
    c.def("__repr__", [ClassId](type&) { return ClassId + "()"; });

    m.def(ClassId.c_str(), [](const Grid::GridProvider<G>&) { return new type(); });

    return c;
  } // ... bind(...)
}; // class InitlessElementFilter


} // namespace bindings
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // PYTHON_DUNE_XT_GRID_FILTERS_ELEMENT_HH
