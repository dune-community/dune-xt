// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017, 2019 - 2020)
//   René Fritze     (2018 - 2020)
//   Tobias Leibner  (2018, 2020)

#ifndef PYTHON_DUNE_XT_GRID_GRIDPROVIDER_HH
#define PYTHON_DUNE_XT_GRID_GRIDPROVIDER_HH

#include <algorithm>

#include <dune/geometry/type.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/pybindxi/functional.h>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/common/parallel/mpi_comm_wrapper.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/functions/generic/grid-function.hh>
#include <dune/xt/functions/visualization.hh>
#include <dune/xt/grid/element.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/exceptions.hh>
#include <dune/xt/grid/filters/intersection.hh>
#include <dune/xt/grid/gridprovider/dgf.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/mapper.hh>
#include <dune/xt/la/container/common/vector/dense.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune::XT::Grid::bindings {

template <class G>
class GridProvider
{
public:
  using type = Grid::GridProvider<G>;
  using bound_type = pybind11::class_<type>;

  using GV = typename type::LeafGridViewType;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "grid_provider",
                         const std::string& grid_id = grid_name<G>::value())
  {
    namespace py = pybind11;
    using namespace pybind11::literals;
    constexpr const int dim = type::GridType::dimension;

    const std::string class_name = class_id + "_" + grid_id;
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), (XT::Common::to_camel_case(class_id) + " (" + grid_id + " variant)").c_str());
    c.def_property_readonly("dimension", [](type&) { return dim; });
    c.def_property_readonly("max_level", &type::max_level);
    c.def(
        "size",
        [](type& self, const int codim) {
          DUNE_THROW_IF(
              codim < 0 || codim > dim, Exceptions::wrong_codimension, "dim = " << dim << "\n   codim = " << codim);
          DUNE_THROW_IF(codim != dim && codim != 0 && !G::LeafGridView::conforming,
                        XT::Common::Exceptions::requirements_not_met,
                        "This is not yet implemented for non-conforming grids and codim " << codim << "!");
          const LeafMultipleCodimMultipleGeomTypeMapper<G> mapper(self.grid(), [codim](GeometryType gt, int dimgrid) {
            return dimgrid - Common::numeric_cast<int>(gt.dim()) == codim;
          });
          return mapper.size();
        },
        "codim"_a);
    c.def(
        "centers",
        [](type& self, const int codim) {
          DUNE_THROW_IF(
              codim < 0 || codim > dim, Exceptions::wrong_codimension, "dim = " << dim << "\n   codim = " << codim);
          DUNE_THROW_IF(codim != dim && codim != 0 && !G::LeafGridView::conforming,
                        XT::Common::Exceptions::requirements_not_met,
                        "This is not yet implemented for non-conforming grids and codim " << codim << "!");
          auto grid_view = self.leaf_view();
          const LeafMultipleCodimMultipleGeomTypeMapper<G> mapper(self.grid(), [codim](GeometryType gt, int dimgrid) {
            return dimgrid - Common::numeric_cast<int>(gt.dim()) == codim;
          });
          auto centers =
              std::make_unique<XT::LA::CommonDenseMatrix<double>>(mapper.size(), Common::numeric_cast<size_t>(dim), 0.);
          for (auto&& element : elements(grid_view)) {
            for (auto&& ii : Common::value_range(element.subEntities(codim))) {
              auto index = sub_entity_index(mapper, element, codim, ii);
              auto center = sub_entity_center(element, codim, ii);
              for (size_t jj = 0; jj < Common::numeric_cast<size_t>(dim); ++jj)
                centers->set_entry(index, jj, center[jj]);
            }
          }
          return centers;
        },
        "codim"_a = 0,
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "inner_intersection_indices",
        [](type& self) {
          DUNE_THROW_IF(!G::LeafGridView::conforming,
                        XT::Common::Exceptions::requirements_not_met,
                        "This is not yet implemented for non-conforming grids!");
          auto grid_view = self.leaf_view();
          const LeafMultipleCodimMultipleGeomTypeMapper<G> mapper(self.grid(), mcmgLayout(Codim<1>()));
          std::set<size_t> global_indices;
          Grid::ApplyOn::InnerIntersections<GV> filter;
          for (auto&& element : elements(grid_view)) {
            for (auto&& intersection : intersections(grid_view, element)) {
              if (filter.contains(grid_view, intersection)) {
                const auto intersection_entity = element.template subEntity<1>(intersection.indexInInside());
                global_indices.insert(Common::numeric_cast<size_t>(mapper.index(intersection_entity)));
              }
            }
          }
          LA::CommonDenseVector<size_t> ret(global_indices.size());
          size_t ii = 0;
          for (const auto& index : global_indices) {
            ret[ii] = index;
            ++ii;
          }
          return ret;
        },
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "inside_element_indices",
        [](type& self) {
          DUNE_THROW_IF(!G::LeafGridView::conforming,
                        XT::Common::Exceptions::requirements_not_met,
                        "This is not yet implemented for non-conforming grids!");
          auto grid_view = self.leaf_view();
          const LeafMultipleCodimMultipleGeomTypeMapper<G> element_mapper(self.grid(), mcmgElementLayout());
          const LeafMultipleCodimMultipleGeomTypeMapper<G> intersection_mapper(self.grid(), mcmgLayout(Codim<1>()));
          LA::CommonDenseVector<size_t> element_indices(intersection_mapper.size(), std::numeric_limits<size_t>::max());
          Grid::ApplyOn::AllIntersectionsOnce<GV> filter;
          for (auto&& element : elements(grid_view)) {
            const auto element_index = Common::numeric_cast<size_t>(element_mapper.index(element));
            for (auto&& intersection : intersections(grid_view, element)) {
              if (filter.contains(grid_view, intersection)) {
                const auto intersection_entity = element.template subEntity<1>(intersection.indexInInside());
                const auto intersection_index =
                    Common::numeric_cast<size_t>(intersection_mapper.index(intersection_entity));
                element_indices[intersection_index] = element_index;
              }
            }
          }
          return element_indices;
        },
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "outside_element_indices",
        [](type& self) {
          DUNE_THROW_IF(!G::LeafGridView::conforming,
                        XT::Common::Exceptions::requirements_not_met,
                        "This is not yet implemented for non-conforming grids!");
          auto grid_view = self.leaf_view();
          const LeafMultipleCodimMultipleGeomTypeMapper<G> element_mapper(self.grid(), mcmgElementLayout());
          const LeafMultipleCodimMultipleGeomTypeMapper<G> intersection_mapper(self.grid(), mcmgLayout(Codim<1>()));
          LA::CommonDenseVector<size_t> element_indices(intersection_mapper.size(), std::numeric_limits<size_t>::max());
          Grid::ApplyOn::InnerIntersectionsOnce<GV> filter;
          for (auto&& element : elements(grid_view)) {
            for (auto&& intersection : intersections(grid_view, element)) {
              if (filter.contains(grid_view, intersection)) {
                const auto intersection_entity = element.template subEntity<1>(intersection.indexInInside());
                const auto intersection_index =
                    Common::numeric_cast<size_t>(intersection_mapper.index(intersection_entity));
                const auto outside_element_index =
                    Common::numeric_cast<size_t>(element_mapper.index(intersection.outside()));
                element_indices[intersection_index] = outside_element_index;
              }
            }
          }
          return element_indices;
        },
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "boundary_intersection_indices",
        [](type& self, const Grid::IntersectionFilter<GV>& filter) {
          DUNE_THROW_IF(!G::LeafGridView::conforming,
                        XT::Common::Exceptions::requirements_not_met,
                        "This is not yet implemented for non-conforming grids!");
          auto grid_view = self.leaf_view();
          const LeafMultipleCodimMultipleGeomTypeMapper<G> mapper(self.grid(), mcmgLayout(Codim<1>()));
          std::set<size_t> global_indices;
          for (auto&& element : elements(grid_view)) {
            for (auto&& intersection : intersections(grid_view, element)) {
              if (filter.contains(grid_view, intersection)) {
                const auto intersection_entity = element.template subEntity<1>(intersection.indexInInside());
                global_indices.insert(Common::numeric_cast<size_t>(mapper.index(intersection_entity)));
              }
            }
          }
          LA::CommonDenseVector<size_t> ret(global_indices.size());
          size_t ii = 0;
          for (const auto& index : global_indices) {
            ret[ii] = index;
            ++ii;
          }
          return ret;
        },
        "intersection_filter"_a = Grid::ApplyOn::BoundaryIntersections<GV>(),
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "visualize",
        [](type& self, const std::string& filename) {
          const LeafMultipleCodimMultipleGeomTypeMapper<G> mapper(self.grid(), [](GeometryType gt, int dimgrid) {
            return dimgrid - Common::numeric_cast<int>(gt.dim()) == 0;
          });
          double element_index = 0; // not thread safe!
          Functions::GenericGridFunction<extract_entity_t<typename G::LeafGridView>> element_index_function(
              /*order=*/[](const auto&) { return 0; },
              /*post_bind=*/
              [&mapper, &element_index](const auto& element) { element_index = mapper.index(element); },
              /*evaluate=*/
              [&element_index](const auto&, const auto&) { return element_index; },
              /*param_type=*/{},
              /*name=*/"Element index");
          Functions::visualize(element_index_function,
                               self.leaf_view(),
                               filename,
                               /*subsampling=*/false);
        },
        "filename"_a,
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "global_refine",
        [](type& self, const int count) { self.global_refine(count); },
        "count"_a = 1,
        py::call_guard<py::gil_scoped_release>());
    c.def("refine_steps_for_half", [](type& /*self*/) { return DGFGridInfo<G>::refineStepsForHalf(); });
    c.def("apply_on_each_element",
          [](type& self, std::function<void(const XT::Grid::extract_entity_t<GV>&)> generic_function) {
            for (auto&& element : elements(self.leaf_view()))
              generic_function(element);
          });
    c.def("apply_on_each_intersection",
          [](type& self, std::function<void(const XT::Grid::extract_intersection_t<GV>&)> generic_function) {
            auto gv = self.leaf_view();
            for (auto&& element : elements(gv))
              for (auto&& intersection : intersections(gv, element))
                generic_function(intersection);
          });

    return c;
  } // ... bind(...)
}; // class GridProvider

} // namespace Dune::XT::Grid::bindings

#endif // PYTHON_DUNE_XT_GRID_GRIDPROVIDER_HH
