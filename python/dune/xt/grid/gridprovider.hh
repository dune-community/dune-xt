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
#include <dune/xt/grid/gridprovider/dgf.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/functions/generic/grid-function.hh>

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
    c.def_property_readonly("dimension", [dim](type&) { return dim; });
    c.def_property_readonly("max_level", &type::max_level);
    c.def(
        "size",
        [dim](type& self, const int codim) {
          DUNE_THROW_IF(
              codim < 0 || codim > dim, Exceptions::wrong_codimension, "dim = " << dim << "\n   codim = " << codim);
          auto grid_view = self.leaf_view();
          MultipleCodimMultipleGeomTypeMapper<decltype(grid_view)> mapper(
              grid_view,
              [codim](GeometryType gt, int dimgrid) { return dimgrid - Common::numeric_cast<int>(gt.dim()) == codim; });
          return mapper.size();
        },
        "codim"_a);
    c.def(
        "centers",
        [dim](type& self, const int codim) {
          DUNE_THROW_IF(codim != 0, NotImplemented, "Only for codim 0 at the moment!");
          DUNE_THROW_IF(
              codim < 0 || codim > dim, Exceptions::wrong_codimension, "dim = " << dim << "\n   codim = " << codim);
          auto grid_view = self.leaf_view();
          MultipleCodimMultipleGeomTypeMapper<decltype(grid_view)> mapper(
              grid_view,
              [codim](GeometryType gt, int dimgrid) { return dimgrid - Common::numeric_cast<int>(gt.dim()) == codim; });
          XT::LA::CommonDenseMatrix<double> centers(mapper.size(), Common::numeric_cast<size_t>(dim), 0.);
          for (auto&& element : elements(grid_view)) {
            auto index = mapper.index(element);
            auto center = element.geometry().center();
            for (size_t jj = 0; jj < Common::numeric_cast<size_t>(dim); ++jj)
              centers.set_entry(index, jj, center[jj]);
          }
          return centers;
        },
        "codim"_a = 0,
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "visualize",
        [](type& self, const std::string& filename, const std::string& layer) {
          DUNE_THROW_IF(layer != "leaf", NotImplemented, "Visualization of level views not implemented yet!");
          auto grid_view = self.leaf_view();
          using GV = decltype(grid_view);
          const MultipleCodimMultipleGeomTypeMapper<GV> mapper(
              grid_view, [](GeometryType gt, int dimgrid) { return dimgrid == Common::numeric_cast<int>(gt.dim()); });
          double element_index = 0;
          Functions::GenericGridFunction<extract_entity_t<GV>> element_index_function(
              /*order=*/[](const auto&) { return 0; },
              /*post_bind=*/[&mapper, &element_index](const auto& element) { element_index = mapper.index(element); },
              /*evaluate=*/[&element_index](const auto&, const auto&) { return element_index; },
              /*param_type=*/{},
              /*name=*/"Element index");
          element_index_function.visualize(grid_view, filename, /*subsampling=*/false);
        },
        "filename"_a,
        "layer"_a = "leaf",
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "global_refine",
        [](type& self, const int count) { self.global_refine(count); },
        "count"_a = 1,
        py::call_guard<py::gil_scoped_release>());
    c.def("refine_steps_for_half", [](type& /*self*/) { return DGFGridInfo<G>::refineStepsForHalf(); });
    return c;
  } // ... bind(...)
}; // class GridProvider


} // namespace bindings
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // PYTHON_DUNE_XT_GRID_GRIDPROVIDER_HH
