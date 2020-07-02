// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017, 2019)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2018, 2020)

#ifndef PYTHON_DUNE_XT_GRID_GRIDPROVIDER_HH
#define PYTHON_DUNE_XT_GRID_GRIDPROVIDER_HH

#include <dune/geometry/type.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/xt/common/parallel/mpi_comm_wrapper.hh>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/exceptions.hh>
#include <dune/xt/grid/gridprovider/provider.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace bindings {


template <class G>
class GridProvider
{
public:
  using type = Grid::GridProvider<G>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "grid_provider",
                         const std::string& grid_id = grid_name<G>::value())
  {
    namespace py = pybind11;
    using namespace pybind11::literals;
    const int dim = type::GridType::dimension;

    const std::string class_name = class_id + "_" + grid_id;
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), (XT::Common::to_camel_case(class_id) + " (" + grid_id + " variant)").c_str());
    c.def_property_readonly("dimension", [dim](type&){return dim;});
    c.def_property_readonly("max_level", &type::max_level);
    c.def("size", [dim](type& self, const int codim){
      DUNE_THROW_IF(codim < 0 || codim > dim,
                    Exceptions::wrong_codimension,
                    "dim = " << dim << "\n   codim = " << codim);
      auto grid_view = self.leaf_view();
      MultipleCodimMultipleGeomTypeMapper<decltype(grid_view)> mapper(grid_view, [codim](GeometryType gt, int dimgrid) {
        return dimgrid - Common::numeric_cast<int>(gt.dim()) == codim;
      });
      return mapper.size();
    }, "codim"_a);
    c.def("visualize", [](type& self, const std::string& filename, const Common::Configuration& boundary_info_cfg){
      self.visualize(filename, boundary_info_cfg);
    },
    "filename"_a,
        "boundary_info_cfg"_a = Common::Configuration());
    c.def("global_refine", [](type& self, const int count){self.global_refine(count);}, "count"_a = 1);
    return c;
  } // ... bind(...)
}; // class GridProvider


} // namespace bindings
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // PYTHON_DUNE_XT_GRID_GRIDPROVIDER_HH
