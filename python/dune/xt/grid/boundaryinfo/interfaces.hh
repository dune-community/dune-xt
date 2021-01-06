// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)
//   René Fritze     (2020)

#ifndef PYTHON_DUNE_XT_GRID_BOUNDARYINFO_INTERFACES_HH
#define PYTHON_DUNE_XT_GRID_BOUNDARYINFO_INTERFACES_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/boundaryinfo/interfaces.hh>
#include <dune/xt/grid/type_traits.hh>
#include <python/dune/xt/common/timedlogging.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune::XT::Grid::bindings {


template <class GV>
class BoundaryInfo
{
  using G = typename GV::Grid;
  static_assert(is_grid<G>::value);
  using I = extract_intersection_t<GV>;

public:
  using type = Grid::BoundaryInfo<I>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id = grid_name<G>::value(),
                         const std::string& layer_id = "",
                         const std::string& class_id = "boundary_info")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def_readonly("logger", &type::logger);
    c.def("__repr__", [](type& self) { return self.str(); });

    return c;
  } // ... bind(...)
}; // class BoundaryInfo


} // namespace Dune::XT::Grid::bindings

#endif // PYTHON_DUNE_XT_GRID_BOUNDARYINFO_INTERFACES_HH
