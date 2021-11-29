// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tim Keil (2020 - 2021)

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
#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/element.hh>
#include <dune/xt/grid/exceptions.hh>
#include <dune/xt/grid/filters/intersection.hh>
#include <dune/xt/grid/gridprovider/dgf.hh>
#include <dune/xt/grid/gridprovider/coupling.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/mapper.hh>
#include <dune/xt/grid/view/coupling.hh>
#include <dune/xt/functions/generic/grid-function.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>


namespace Dune::XT::Grid {

/**
 * TODO: MOVE THIS ELSEWHERE !
 */
template <class MacroGV, class MicroGV>
class MacroGridBasedBoundaryInfo : public Dune::XT::Grid::BoundaryInfo<Dune::XT::Grid::extract_intersection_t<MicroGV>>
{
  using BaseType = Dune::XT::Grid::BoundaryInfo<Dune::XT::Grid::extract_intersection_t<MicroGV>>;

public:
  using typename BaseType::IntersectionType;
  using MacroBoundaryInfoType = Dune::XT::Grid::BoundaryInfo<Dune::XT::Grid::extract_intersection_t<MacroGV>>;
  using MacroElementType = Dune::XT::Grid::extract_entity_t<MacroGV>;

  MacroGridBasedBoundaryInfo(const MacroGV macro_grid_view,
                             const MacroElementType macro_element,
                             const MacroBoundaryInfoType& macro_boundary_info)
    : macro_grid_view_(macro_grid_view)
    , macro_element_(macro_element)
    , macro_boundary_info_(macro_boundary_info)
  {}

  const Dune::XT::Grid::BoundaryType& type(const IntersectionType& intersection) const override final
  {
    // find out if this micro intersection lies within the macro element or on a macro intersection
    for (auto&& macro_intersection : intersections(macro_grid_view_, macro_element_)) {
      const size_t num_corners = intersection.geometry().corners();
      size_t num_corners_inside = 0;
      size_t num_corners_outside = 0;
      for (size_t cc = 0; cc < num_corners; ++cc) {
        const auto micro_corner = intersection.geometry().corner(cc);
        if (XT::Grid::contains(macro_intersection, micro_corner))
          ++num_corners_inside;
        else
          ++num_corners_outside;
      }
      if (num_corners_inside == num_corners && num_corners_outside == 0) {
        // we found the macro intersection this micro intersection belongs to
        return macro_boundary_info_.type(macro_intersection);
      }
    }
    // we could not find a macro intersection this micro intersection belongs to
    return no_boundary_;
  } // ... type(...)

  const MacroGV macro_grid_view_;
  const MacroElementType macro_element_;
  const MacroBoundaryInfoType& macro_boundary_info_;
  const XT::Grid::NoBoundary no_boundary_;
}; // class MacroGridBasedBoundaryInfo


} // namespace Dune::XT::Grid

namespace Dune::XT::Grid::bindings {


template <class MG, class G>
class GluedGridProvider
{
public:
  using type = Grid::DD::Glued<MG, G, XT::Grid::Layers::leaf>;
  using bound_type = pybind11::class_<type>;

  using GV = typename type::LocalViewType;
  using MacroGV = typename type::MacroGridViewType;
  using ElementType = typename type::MacroEntityType;
  using CouplingGridViewType = Dune::XT::Grid::CouplingGridView<type>;
  using IntersectionType = typename MacroGV::Intersection;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "glued_grid_provider",
                         const std::string& macro_grid_id = grid_name<MG>::value(),
                         const std::string& micro_grid_id = grid_name<G>::value())
  {
    namespace py = pybind11;
    using namespace pybind11::literals;
    constexpr const int dim = type::dimDomain;

    const std::string class_name = class_id + "_" + macro_grid_id + "_" + micro_grid_id;
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(
        m,
        ClassName.c_str(),
        (XT::Common::to_camel_case(class_id) + " (" + macro_grid_id + "_" = micro_grid_id + " variant)").c_str());
    // dim cannot be binded directly because it is static const
    c.def_property_readonly("dimension", [](type&) { return dim; });
    c.def("local_grid", py::overload_cast<size_t>(&type::local_grid));
    c.def_property_readonly("num_subdomains", &type::num_subdomains);
    c.def_property_readonly("boundary_subdomains", [](type& self) {
      std::vector<size_t> boundary_subdomains;
      for (auto&& macro_element : elements(self.macro_grid_view()))
        if (self.boundary(macro_element))
          boundary_subdomains.push_back(self.subdomain(macro_element));
      return boundary_subdomains;
    });
    c.def(
        "neighbors",
        [](type& self, const size_t ss) {
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
        },
        "ss"_a);
    c.def(
        "coupling_grid",
        [](type& self, const size_t ss, const size_t nn) {
          for (auto&& inside_macro_element : elements(self.macro_grid_view())) {
            if (self.subdomain(inside_macro_element) == ss) {
              // this is the subdomain we are interested in
              bool found_correct_macro_intersection = false;
              for (auto&& macro_intersection : intersections(self.macro_grid_view(), inside_macro_element)) {
                if (macro_intersection.neighbor()) {
                  const auto outside_macro_element = macro_intersection.outside();
                  if (self.subdomain(outside_macro_element) == nn) {
                    found_correct_macro_intersection = true;
                    // these are the subdomains we are interested in
                    auto cgv = Dune::XT::Grid::make_coupling_grid_view<type, ElementType, IntersectionType>(
                        inside_macro_element, outside_macro_element, self, macro_intersection);
                    return std::make_unique<Dune::XT::Grid::CouplingGridProvider<CouplingGridViewType>>(cgv);
                  }
                }
              }
              DUNE_THROW_IF(!found_correct_macro_intersection,
                            XT::Common::Exceptions::index_out_of_range,
                            "ss = " << ss << "\n   nn = " << nn);
            }
          }
        },
        "ss"_a,
        "nn"_a);
    c.def(
        "macro_based_boundary_info",
        [](type& self, const size_t ss, const XT::Grid::BoundaryInfo<IntersectionType>& macro_boundary_info) {
          DUNE_THROW_IF(ss >= self.num_subdomains(),
                        XT::Common::Exceptions::index_out_of_range,
                        "ss = " << ss << "\n   self.num_subdomains() = " << self.num_subdomains());
          for (auto&& macro_element : elements(self.macro_grid_view())) {
            if (self.subdomain(macro_element) == ss) {
              // this is the subdomain we are interested in, create space
              return MacroGridBasedBoundaryInfo<MacroGV, GV>(
                  self.macro_grid_view(), macro_element, macro_boundary_info);
            }
          }
        },
        "ss"_a,
        "macro_boundary_info"_a);
    c.def("write_global_visualization", &type::write_global_visualization);
    return c;
  } // ... bind(...)
}; // class GluedGridProvider

template <class CGV>
class CouplingGridProvider
{
public:
  using type = Dune::XT::Grid::CouplingGridProvider<CGV>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id = grid_name<typename CGV::MacroGridType>::value(),
                         const std::string& class_id = "coupling_grid_provider")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const std::string class_name = class_id + "_" + grid_id;
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), (XT::Common::to_camel_case(class_id) + " (" + grid_id + " variant)").c_str());
    c.def_property_readonly("static_id", [](type& self) { return self.static_id(); });
    return c;
  } // ... bind(...)
}; // class CouplingGridProvider

template <class G>
struct MacroGridBasedBoundaryInfo
{
  using GV = typename G::LeafGridView;
  using I = Dune::XT::Grid::extract_intersection_t<GV>;
  using E = Dune::XT::Grid::extract_entity_t<GV>;
  using type = Dune::XT::Grid::MacroGridBasedBoundaryInfo<GV, GV>; // assumes to have the same GV in micro und macro
  using base_type = Dune::XT::Grid::BoundaryInfo<I>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id = grid_name<G>::value(),
                         const std::string& class_id = "macro_grid_based_boundary_info")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;
    const std::string class_name = class_id + "_" + grid_id;
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), (XT::Common::to_camel_case(class_id) + " (" + grid_id + " variant)").c_str());
    return c;
  } // ... bind(...)

  // no factory needed: will only be called from GluedGridProvider

}; // struct MacroGridBasedBoundaryInfo

} // namespace Dune::XT::Grid::bindings

/**
 * \note Available grid types for DD::Glued. So far, we only use Alugrid and Yasp
 */

using AvailableGridGlueGridTypes = std::tuple<YASP_2D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                              ,
                                              ALU_2D_SIMPLEX_CONFORMING,
                                              ALU_2D_CUBE
#endif
                                              >;

using GridGlue2dYaspYasp =
    Dune::XT::Grid::DD::Glued<YASP_2D_EQUIDISTANT_OFFSET, YASP_2D_EQUIDISTANT_OFFSET, Dune::XT::Grid::Layers::leaf>;
using CouplingGridView2dYaspYasp = Dune::XT::Grid::CouplingGridView<GridGlue2dYaspYasp>;
#if HAVE_DUNE_ALUGRID
using GridGlue2dAluSimplexConformingAluSimplexConforming =
    Dune::XT::Grid::DD::Glued<ALU_2D_SIMPLEX_CONFORMING, ALU_2D_SIMPLEX_CONFORMING, Dune::XT::Grid::Layers::leaf>;
using GridGlue2dAluConformingAluConforming =
    Dune::XT::Grid::DD::Glued<ALU_2D_CUBE, ALU_2D_CUBE, Dune::XT::Grid::Layers::leaf>;
using CouplingGridView2dAluSimplexConformingAluSimplexConforming =
    Dune::XT::Grid::CouplingGridView<GridGlue2dAluSimplexConformingAluSimplexConforming>;
using CouplingGridView2dAluConformingAluConforming =
    Dune::XT::Grid::CouplingGridView<GridGlue2dAluConformingAluConforming>;
#endif

using AvailableCouplingGridViewTypes = std::tuple<CouplingGridView2dYaspYasp
#if HAVE_DUNE_ALUGRID
                                                  ,
                                                  CouplingGridView2dAluSimplexConformingAluSimplexConforming,
                                                  CouplingGridView2dAluConformingAluConforming
#endif
                                                  >;


template <class GridTypes = AvailableGridGlueGridTypes>
struct GluedGridProvider_for_all_available_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;

  static void bind(pybind11::module& m)
  {
    Dune::XT::Grid::bindings::GluedGridProvider<G, G>::bind(m);
    GluedGridProvider_for_all_available_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct GluedGridProvider_for_all_available_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


template <class CouplingGridViewTypes = AvailableCouplingGridViewTypes>
struct CouplingGridProvider_for_all_available_grids
{
  using CGV = Dune::XT::Common::tuple_head_t<CouplingGridViewTypes>;

  static void bind(pybind11::module& m)
  {
    Dune::XT::Grid::bindings::CouplingGridProvider<CGV>::bind(m);
    CouplingGridProvider_for_all_available_grids<Dune::XT::Common::tuple_tail_t<CouplingGridViewTypes>>::bind(m);
  }
};

template <>
struct CouplingGridProvider_for_all_available_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct MacroGridBasedBoundaryInfo_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;

  static void bind(pybind11::module& m)
  {
    Dune::XT::Grid::bindings::MacroGridBasedBoundaryInfo<G>::bind(m);
    MacroGridBasedBoundaryInfo_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct MacroGridBasedBoundaryInfo_for_all_grids<Dune::XT::Common::tuple_null_type>
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

#if HAVE_DUNE_GRID_GLUE
  GluedGridProvider_for_all_available_grids<>::bind(m);
  CouplingGridProvider_for_all_available_grids<>::bind(m);
  MacroGridBasedBoundaryInfo_for_all_grids<>::bind(m);
#endif
}
