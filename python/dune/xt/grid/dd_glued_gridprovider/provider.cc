// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:

#include "config.h"

#include <algorithm>

#include <dune/geometry/type.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/common/parallel/mpi_comm_wrapper.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/la/container/common/vector/dense.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/element.hh>
#include <dune/xt/grid/exceptions.hh>
#include <dune/xt/grid/gridprovider/dgf.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/grid/filters/intersection.hh>
#include <dune/xt/grid/mapper.hh>
#include <dune/xt/functions/generic/grid-function.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace bindings {


template <class G>
class GluedGridProvider
{
public:
  using type = Grid::DD::Glued<G, G, XT::Grid::Layers::leaf>;
  using bound_type = pybind11::class_<type>;

  using GV = typename type::LocalViewType;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "glued_grid_provider",
                         const std::string& grid_id = grid_name<G>::value())
  {
    namespace py = pybind11;
    using namespace pybind11::literals;
    constexpr const int dim = type::dimDomain;

    const std::string class_name = class_id + "_" + grid_id;
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), (XT::Common::to_camel_case(class_id) + " (" + grid_id + " variant)").c_str());
    c.def_property_readonly("dimension", [](type&) { return dim; });
    c.def(
          "local_grid",
          [](type& self, const size_t macro_entity_index) {
          auto local_grid = self.local_grid(macro_entity_index);
          return local_grid;
          },
          "macro_entity_index"_a);
    c.def_property_readonly("num_subdomains", [](type& self) { return self.num_subdomains(); });
    c.def_property_readonly("boundary_subdomains", [](type& self) {
        std::vector<size_t> boundary_subdomains;
        for (auto&& macro_element : elements(self.macro_grid_view()))
          if (self.boundary(macro_element))
            boundary_subdomains.push_back(self.subdomain(macro_element));
          return boundary_subdomains;
      });
    c.def("neighbors", [](type& self, const size_t ss) {
        DUNE_THROW_IF(ss >= self.num_subdomains(),
                        XT::Common::Exceptions::index_out_of_range,
                        "ss = " << ss << "\n   self.num_subdomains() = " << self.num_subdomains());
         std::vector<size_t> neighboring_subdomains;
         for (auto&& macro_element : elements(self.macro_grid_view())) {
           if (self.subdomain(macro_element) == ss) {
             for (auto&& macro_intersection : intersections(self.macro_grid_view(), macro_element))
               if (macro_intersection.neighbor())
                 neighboring_subdomains.push_back(self.subdomain(macro_intersection.outside()));
             break;
           }
         }
         return neighboring_subdomains;
    });

    c.def(
        "write_global_visualization",
        [](type& self,
           const std::string& filename_prefix,
           const XT::Functions::FunctionInterface<dim>& func,
           const std::list<int>& subdomains)
            { self.write_global_visualization(filename_prefix, func, subdomains); },
        "filename_prefix"_a,
        "function"_a,
        "subdomains"_a);
    return c;
  } // ... bind(...)
}; // class GridProvider


} // namespace bindings
} // namespace Grid
} // namespace XT
} // namespace Dune



#include <dune/xt/grid/grids.hh>
// #include <python/dune/xt/grid/filters/intersection.hh>


//template <template <class> class Filter, class GridTypes = Dune::XT::Grid::AvailableGridTypes>
//struct InitlessIntersectionFilter_for_all_grids
//{
//  static void bind(pybind11::module& m, const std::string& class_id)
//  {
//    Dune::XT::Grid::bindings::InitlessIntersectionFilter<Filter, typename GridTypes::head_type>::bind(m, class_id);
//    InitlessIntersectionFilter_for_all_grids<Filter, typename GridTypes::tail_type>::bind(m, class_id);
//  }
//};

//template <template <class> class Filter>
//struct InitlessIntersectionFilter_for_all_grids<Filter, boost::tuples::null_type>
//{
//  static void bind(pybind11::module& /*m*/, const std::string& /*class_id*/) {}
//};


//template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
//struct CustomBoundaryIntersectionFilter_for_all_grids
//{
//  static void bind(pybind11::module& m)
//  {
//    Dune::XT::Grid::bindings::CustomBoundaryIntersectionsFilter<typename GridTypes::head_type>::bind(m);
//    CustomBoundaryIntersectionFilter_for_all_grids<typename GridTypes::tail_type>::bind(m);
//  }
//};

//template <>
//struct CustomBoundaryIntersectionFilter_for_all_grids<boost::tuples::null_type>
//{
//  static void bind(pybind11::module& /*m*/) {}
//};


template <class GridTypes = Dune::XT::Grid::Available2dGridTypes>  // grid-glue only working 2d
struct GluedGridProvider_for_all_grids
{
  static void bind(pybind11::module& m)
  {
    Dune::XT::Grid::bindings::GluedGridProvider<Dune::XT::Common::tuple_head_t<GridTypes>>::bind(m);
    GluedGridProvider_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct GluedGridProvider_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_grid_dd_glued_gridprovider_provider, m)
{
  namespace py = pybind11;
  using namespace Dune::XT::Grid;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid._grid_boundaryinfo_interfaces");
  py::module::import("dune.xt.grid._grid_boundaryinfo_types");
  py::module::import("dune.xt.grid._grid_filters_base");

//#define BIND_(NAME) InitlessIntersectionFilter_for_all_grids<ApplyOn::NAME>::bind(m, std::string("ApplyOn") + #NAME)

//  BIND_(AllIntersections);
//  BIND_(AllIntersectionsOnce);
//  BIND_(NoIntersections);
//  BIND_(InnerIntersections);
//  BIND_(InnerIntersectionsOnce);
//  //  BIND_(PartitionSetInnerIntersectionsOnce); <- requires partition set as template argument
//  BIND_(BoundaryIntersections);
//  BIND_(NonPeriodicBoundaryIntersections);
//  BIND_(PeriodicBoundaryIntersections);
//  BIND_(PeriodicBoundaryIntersectionsOnce);
//  //  BIND_(GenericFilteredIntersections); <- requires lambda in init
//  //  BIND_(CustomBoundaryAndProcessIntersections); <- requires boundary type and info in init
//  BIND_(ProcessIntersections);

//#undef BIND_

//  CustomBoundaryIntersectionFilter_for_all_grids<>::bind(m);

  GluedGridProvider_for_all_grids<>::bind(m);
}
