// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_XT_GRID_BOUNDARYINFO_INTERFACES_HH
#define PYTHON_DUNE_XT_GRID_BOUNDARYINFO_INTERFACES_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/boundaryinfo/interfaces.hh>
#include <python/dune/xt/common/timedlogging.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace bindings {


template <class G>
class BoundaryInfo
{
  static_assert(is_grid<G>::value, "");
  using GV = typename G::LeafGridView;
  using I = extract_intersection_t<GV>;

public:
  using type = Grid::BoundaryInfo<I>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "boundary_info",
                         const std::string& grid_id = grid_name<G>::value())
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def_readonly("logger", &type::logger);
    c.def("__repr__", [](type& self) { return self.str(); });

    return c;
  } // ... bind(...)
}; // class BoundaryInfo


} // namespace bindings
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // PYTHON_DUNE_XT_GRID_BOUNDARYINFO_INTERFACES_HH
