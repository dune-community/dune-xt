// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_XT_GRID_WALKER_APPLYON_BINDGINS_HH
#define DUNE_XT_GRID_WALKER_APPLYON_BINDGINS_HH

#include <sstream>
#include <type_traits>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/operators.h>

#include <dune/xt/grid/type_traits.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/layers.bindings.hh>
#include <dune/xt/grid/boundaryinfo/interfaces.hh>

#include <dune/xt/grid/filters.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace bindings {


template <class Imp, bool ctor_expects_boundary_info>
class IntersectionFilter
{
  typedef typename Imp::GridViewType GL;
  typedef typename extract_grid<GL>::type G;
  typedef XT::Grid::IntersectionFilter<GL> InterfaceType;

public:
  typedef Imp type;
  typedef pybind11::class_<type, InterfaceType> bound_type;

private:
  typedef XT::Grid::BoundaryInfo<extract_intersection_t<GL>> BoundaryInfoType;

  static std::string makename(std::string class_name, std::string layer_name)
  {
    return std::string("make_apply_on_") + class_name + "_" + layer_name + "_"
           + XT::Grid::bindings::grid_name<G>::value();
  }

  template <bool with_bi = ctor_expects_boundary_info, bool anything = true>
  struct addbind // with_bi = false
  {
    void operator()(pybind11::module& m, bound_type& c, const std::string& class_name, const std::string& layer_name)
    {
      c.def(pybind11::init<>());


      m.def(makename(class_name, layer_name).c_str(), []() { return type(); });
    }
  };

  template <bool anything>
  struct addbind<true, anything>
  {
    void operator()(pybind11::module& m, bound_type& c, const std::string& class_name, const std::string& layer_name)
    {
      using namespace pybind11::literals;

      c.def(pybind11::init<const BoundaryInfoType&, XT::Grid::BoundaryType*&&>());

      m.def(makename(class_name, layer_name).c_str(),
            [](const BoundaryInfoType& boundary_info, XT::Grid::BoundaryType*&& boundary_type) {
              return type(boundary_info, std::move(boundary_type));
            },
            "boundary_info"_a,
            "boundary_type"_a);
    }
  };

public:
  static bound_type bind(pybind11::module& m, const std::string& class_name, const std::string& layer_name)
  {
    namespace py = pybind11;

    const auto grid_name = XT::Grid::bindings::grid_name<G>::value();
    const auto InterfaceName = Common::to_camel_case("ApplyOnWhichIntersection_" + grid_name + "_" + layer_name);

    // bind interface, guard since we might not be the first to do so for this intersection
    try {
      py::class_<InterfaceType>(m, InterfaceName.c_str(), InterfaceName.c_str());
    } catch (std::runtime_error&) {
    }

    // bind class
    const auto ClassName = Common::to_camel_case("apply_on_" + class_name + "_" + grid_name + "_" + layer_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    addbind<>()(m, c, class_name, layer_name);

    return c;
  } // ... bind(...)
}; // class WhichIntersection


} // namespace bindings
} // namespace Grid
} // namespace XT
} // namespace Dune


// begin: this is what we need for the .so

#define _DUNE_XT_GRID_WALKER_APPLYON_BIND(_m, _W, _w, _G, _layer, _backend, _class_name)                               \
  Dune::XT::Grid::bindings::IntersectionFilter<                                                                        \
      _W<typename Dune::XT::Grid::Layer<_G,                                                                            \
                                        Dune::XT::Grid::Layers::_layer,                                                \
                                        Dune::XT::Grid::Backends::_backend,                                            \
                                        Dune::XT::Grid::DD::SubdomainGrid<_G>>::type>,                                 \
      _w>::bind(_m,                                                                                                    \
                _class_name,                                                                                           \
                Dune::XT::Grid::layer_names[Dune::XT::Grid::Layers::_layer] + "_"                                      \
                    + Dune::XT::Grid::bindings::backend_name<Dune::XT::Grid::Backends::_backend>::value())

/*#if HAVE_ALBERTA
#define _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALBERTA(_m, _W, _w, _layer, _backend, _class_name) \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND(_m, _W, _w, ALBERTA_2D, _layer, _backend, _class_name)
#else*/
#define _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALBERTA(_m, _W, _w, _layer, _backend, _class_name)
//#endif

#if HAVE_DUNE_ALUGRID
#  define _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALU(_m, _W, _w, _layer, _backend, _class_name)                             \
    _DUNE_XT_GRID_WALKER_APPLYON_BIND(_m, _W, _w, ALU_2D_SIMPLEX_CONFORMING, _layer, _backend, _class_name)
#else
#  define _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALU(_m, _W, _w, _layer, _backend, _class_name)
#endif

/*#if HAVE_DUNE_UGGRID
#define _DUNE_XT_GRID_WALKER_APPLYON_BIND_UG(_m, _W, _w, _layer, _backend, _class_name) \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND(_m, _W, _w, UG_2D, _layer, _backend, _class_name)
#else*/
#define _DUNE_XT_GRID_WALKER_APPLYON_BIND_UG(_m, _W, _w, _layer, _backend, _class_name)
//#endif

#define _DUNE_XT_GRID_WALKER_APPLYON_BIND_YASP(_m, _W, _w, _layer, _backend, _class_name)                              \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND(_m, _W, _w, YASP_1D_EQUIDISTANT_OFFSET, _layer, _backend, _class_name);            \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND(_m, _W, _w, YASP_2D_EQUIDISTANT_OFFSET, _layer, _backend, _class_name)

#define _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALL_GRIDS(_m, _W, _w, _layer, _backend, _class_name)                         \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALBERTA(_m, Dune::XT::Grid::ApplyOn::_W, _w, _layer, _backend, _class_name);       \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALU(_m, Dune::XT::Grid::ApplyOn::_W, _w, _layer, _backend, _class_name);           \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND_UG(_m, Dune::XT::Grid::ApplyOn::_W, _w, _layer, _backend, _class_name);            \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND_YASP(_m, Dune::XT::Grid::ApplyOn::_W, _w, _layer, _backend, _class_name)

#define _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALL(_m, _W, _w, _class_name)                                                 \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALL_GRIDS(_m, _W, _w, leaf, view, _class_name);                                    \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALL_GRIDS(_m, _W, _w, level, view, _class_name);                                   \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALL_GRIDS(_m, _W, _w, dd_subdomain, view, _class_name);                            \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALL_GRIDS(_m, _W, _w, dd_subdomain_boundary, view, _class_name);                   \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALL_GRIDS(_m, _W, _w, dd_subdomain_coupling, view, _class_name)

#define DUNE_XT_GRID_WALKER_APPLYON_BIND(_m)                                                                           \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALL(_m, AllIntersections, false, "all_intersections");                             \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALL(_m, InnerIntersections, false, "inner_intersections");                         \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALL(_m, InnerIntersectionsOnce, false, "inner_intersections_once");                \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALL(_m, BoundaryIntersections, false, "boundary_intersections");                   \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALL(                                                                               \
      _m, NonPeriodicBoundaryIntersections, false, "non_periodic_boundary_intersections");                             \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALL(_m, PeriodicBoundaryIntersections, false, "periodic_intersections");           \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALL(_m, PeriodicBoundaryIntersectionsOnce, false, "periodic_intersections_once");  \
  _DUNE_XT_GRID_WALKER_APPLYON_BIND_ALL(                                                                               \
      _m, CustomBoundaryAndProcessIntersections, true, "custom_boundary_and_process_intersections");

// end: this is what we need for the .so


#endif // DUNE_XT_GRID_WALKER_APPLYON_BINDGINS_HH
