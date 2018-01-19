// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2017 - 2018)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_BINDINGS_HH
#define DUNE_XT_GRID_BOUNDARYINFO_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <sstream>
#include <type_traits>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/configuration.pbh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/layers.bindings.hh>

#include "boundaryinfo.hh"

namespace Dune {
namespace XT {
namespace Grid {
namespace bindings {
namespace internal {


template <class I, class GP, Layers layer, bool anything = true>
class BoundaryInfoFactory
{
public:
  static void bind(pybind11::module& m)
  {
    using namespace pybind11::literals;

    try { // guard since we might not be the first to do so for this grid/intersection
      m.def(std::string("available_boundary_infos_on_" + layer_name<layer>::value() + "_layer").c_str(),
            [](const GP& /*grid_provider*/) { return XT::Grid::BoundaryInfoFactory<I>::available(); },
            "grid_provider"_a);
      m.def(std::string("default_boundary_info_config_on_" + layer_name<layer>::value() + "_layer").c_str(),
            [](const GP& /*grid_provider*/, const std::string& type) {
              return XT::Grid::BoundaryInfoFactory<I>::default_config(type);
            },
            "grid_provider"_a,
            "type"_a);
      m.def(std::string("make_boundary_info_on_" + layer_name<layer>::value() + "_layer").c_str(),
            [](const GP& /*grid_provider*/, const std::string& type, const Common::Configuration& cfg) {
              return XT::Grid::BoundaryInfoFactory<I>::create(type, cfg).release();
            },
            "grid_provider"_a,
            "type"_a,
            "cfg"_a = Common::Configuration());
      m.def(std::string("make_boundary_info_on_" + layer_name<layer>::value() + "_layer").c_str(),
            [](const GP& /*grid_provider*/, const Common::Configuration& cfg) {
              return XT::Grid::BoundaryInfoFactory<I>::create(cfg).release();
            },
            "grid_provider"_a,
            "cfg"_a);
    } catch (std::runtime_error&) {
    }
  } // ... bind(...)
}; // class BoundaryInfoFactory


template <class I, class G, bool anything>
class BoundaryInfoFactory<I, GridProvider<G>, Layers::dd_subdomain, anything>
{
public:
  static void bind(pybind11::module& /*m*/)
  {
  }
};

template <class I, class G, bool anything>
class BoundaryInfoFactory<I, GridProvider<G>, Layers::dd_subdomain_boundary, anything>
{
public:
  static void bind(pybind11::module& /*m*/)
  {
  }
};

template <class I, class G, bool anything>
class BoundaryInfoFactory<I, GridProvider<G>, Layers::dd_subdomain_coupling, anything>
{
public:
  static void bind(pybind11::module& /*m*/)
  {
  }
};


} // namespace internal


template <class Imp, class G, Layers layer>
class BoundaryInfo
{
  typedef typename Imp::IntersectionType I;
  typedef XT::Grid::BoundaryInfo<I> InterfaceType;

public:
  typedef Imp type;
  typedef pybind11::class_<type, InterfaceType> bound_type;

  static void bind(pybind11::module& m, const std::string& class_name, const std::string& layer_name)
  {
    namespace py = pybind11;

    const auto grid_name = bindings::grid_name<G>::value();
    const auto InterfaceName = Common::to_camel_case("BoundaryInfo_" + layer_name + "_" + grid_name);

    // bind interface, guard since we might not be the first to do so for this intersection
    try {
      py::class_<InterfaceType>(m, InterfaceName.c_str(), InterfaceName.c_str(), py::metaclass());
    } catch (std::runtime_error&) {
    }

    // bind factory
    internal::BoundaryInfoFactory<I, GridProvider<G>, layer>::bind(m);
    internal::BoundaryInfoFactory<I, GridProvider<G, DD::SubdomainGrid<G>>, layer>::bind(m);

    // bind class, guard since we might not be the first to do so for this intersection
    try {
      const auto ClassName = Common::to_camel_case(class_name + "_" + layer_name + "_" + grid_name);
      bound_type c(m, ClassName.c_str(), ClassName.c_str(), py::metaclass());
      c.def_property_readonly_static("static_id", [](const type& /*self*/) { return type::static_id(); });
      c.def("__repr__", [ClassName](const type& /*self*/) { return ClassName; });
    } catch (std::runtime_error&) {
    }
  } // ... bind(...)
}; // class BoundaryInfo


} // namespace bindings
} // namespace Grid
} // namespace XT
} // namespace Dune


// begin: this is what we need for the lib

#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, _G, _layer, _backend)                                          \
  prefix class Dune::XT::Grid::bindings::                                                                              \
      BoundaryInfo<_B<Dune::XT::Grid::extract_intersection_t<                                                          \
                       typename Dune::XT::Grid::Layer<_G,                                                              \
                                                      Dune::XT::Grid::Layers::_layer,                                  \
                                                      Dune::XT::Grid::Backends::_backend,                              \
                                                      Dune::XT::Grid::DD::SubdomainGrid<_G>>::type>>,                  \
                   _G,                                                                                                 \
                   Dune::XT::Grid::Layers::_layer>
#if HAVE_DUNE_FEM
#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_YASP(prefix, _B)                                                           \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, YASP_1D_EQUIDISTANT_OFFSET, leaf, view);                             \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, YASP_2D_EQUIDISTANT_OFFSET, leaf, view);                             \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, YASP_1D_EQUIDISTANT_OFFSET, dd_subdomain, part);                     \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, YASP_2D_EQUIDISTANT_OFFSET, dd_subdomain, part)
#else // HAVE_DUNE_FEM
#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_YASP(prefix, _B)                                                           \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, YASP_1D_EQUIDISTANT_OFFSET, leaf, view);                             \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, YASP_2D_EQUIDISTANT_OFFSET, leaf, view);
#endif // HAVE_DUNE_FEM

#if HAVE_DUNE_ALUGRID
#if HAVE_DUNE_FEM
#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_ALU(prefix, _B)                                                            \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, ALU_2D_SIMPLEX_CONFORMING, leaf, view);                              \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, ALU_2D_SIMPLEX_CONFORMING, level, view);                             \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, ALU_2D_SIMPLEX_CONFORMING, dd_subdomain, part);                      \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, ALU_2D_SIMPLEX_CONFORMING, dd_subdomain_boundary, part)
#else // HAVE_DUNE_FEM
#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_ALU(prefix, _B)                                                            \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, ALU_2D_SIMPLEX_CONFORMING, leaf, view);                              \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, ALU_2D_SIMPLEX_CONFORMING, level, view);
#endif // HAVE_DUNE_FEM
#else // HAVE_DUNE_ALUGRID
#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_ALU(prefix, _B)
#endif

//#if HAVE_DUNE_UGGRID || HAVE_UG
//#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_UG(prefix, _B)                                                           \
//  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, UG_2D, leaf, view);                                                \
//  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, UG_2D, level, view);                                               \
//  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, UG_2D, dd_subdomain, part)
//#else
#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_UG(prefix, _B)
//#endif

//#if HAVE_ALBERTA
//#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_ALBERTA(prefix, _B)                                                      \
//  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, ALBERTA_2D, leaf, view);                                           \
//  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix, _B, ALBERTA_2D, dd_subdomain, part)
//#else
#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_ALBERTA(prefix, _B)
//#endif

#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_ALL(prefix, _B)                                                            \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_YASP(prefix, _B);                                                                \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_UG(prefix, _B);                                                                  \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_ALBERTA(prefix, _B);                                                             \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_ALU(prefix, _B)

#define DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(prefix)                                                                     \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_ALL(prefix, Dune::XT::Grid::AllDirichletBoundaryInfo);                           \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_ALL(prefix, Dune::XT::Grid::AllNeumannBoundaryInfo);                             \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_ALL(prefix, Dune::XT::Grid::BoundarySegmentIndexBasedBoundaryInfo);              \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB_ALL(prefix, Dune::XT::Grid::NormalBasedBoundaryInfo)

DUNE_XT_GRID_BOUNDARYINFO_BIND_LIB(extern template);

// end: this is what we need for the lib


// begin: this is what we need for the .so

#define _DUNE_XT_GRID_BOUNDARYINFO_BIND(_m, _B, _G, _layer, _backend, _class_name, _layer_name)                        \
  Dune::XT::Grid::bindings::                                                                                           \
      BoundaryInfo<_B<Dune::XT::Grid::extract_intersection_t<                                                          \
                       typename Dune::XT::Grid::Layer<_G,                                                              \
                                                      Dune::XT::Grid::Layers::_layer,                                  \
                                                      Dune::XT::Grid::Backends::_backend,                              \
                                                      Dune::XT::Grid::DD::SubdomainGrid<_G>>::type>>,                  \
                   _G,                                                                                                 \
                   Dune::XT::Grid::Layers::_layer>::bind(_m, _class_name, _layer_name)

#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_YASP(_m, _B, _class_name)                                                      \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND(_m, _B, YASP_1D_EQUIDISTANT_OFFSET, leaf, view, _class_name, "");                    \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND(                                                                                     \
      _m, _B, YASP_1D_EQUIDISTANT_OFFSET, dd_subdomain, part, _class_name, "dd_subdomain");                            \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND(                                                                                     \
      _m, _B, YASP_1D_EQUIDISTANT_OFFSET, dd_subdomain_boundary, part, _class_name, "dd_subdomain_boundary");          \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND(                                                                                     \
      _m, _B, YASP_1D_EQUIDISTANT_OFFSET, dd_subdomain_coupling, part, _class_name, "dd_subdomain_coupling");          \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND(_m, _B, YASP_2D_EQUIDISTANT_OFFSET, leaf, view, _class_name, "");                    \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND(                                                                                     \
      _m, _B, YASP_2D_EQUIDISTANT_OFFSET, dd_subdomain, part, _class_name, "dd_subdomain");                            \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND(                                                                                     \
      _m, _B, YASP_2D_EQUIDISTANT_OFFSET, dd_subdomain_boundary, part, _class_name, "dd_subdomain_boundary");          \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND(                                                                                     \
      _m, _B, YASP_2D_EQUIDISTANT_OFFSET, dd_subdomain_coupling, part, _class_name, "dd_subdomain_coupling")

#if HAVE_DUNE_ALUGRID
#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_ALU(_m, _B, _class_name)                                                       \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND(_m, _B, ALU_2D_SIMPLEX_CONFORMING, leaf, view, _class_name, "leaf");                 \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND(_m, _B, ALU_2D_SIMPLEX_CONFORMING, level, view, _class_name, "level");               \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND(_m, _B, ALU_2D_SIMPLEX_CONFORMING, dd_subdomain, part, _class_name, "dd_subdomain"); \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND(                                                                                     \
      _m, _B, ALU_2D_SIMPLEX_CONFORMING, dd_subdomain_boundary, part, _class_name, "dd_subdomain_boundary");           \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND(                                                                                     \
      _m, _B, ALU_2D_SIMPLEX_CONFORMING, dd_subdomain_coupling, part, _class_name, "dd_subdomain_coupling")
#else
#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_ALU(_m, _B, _class_name)
#endif

//#if HAVE_DUNE_UGGRID
//#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_UG(_m, _B, _class_name)                                                      \
//  _DUNE_XT_GRID_BOUNDARYINFO_BIND(_m, _B, UG_2D, leaf, view, _class_name, "leaf");                                   \
//  _DUNE_XT_GRID_BOUNDARYINFO_BIND(_m, _B, UG_2D, level, view, _class_name, "level");                                 \
//  _DUNE_XT_GRID_BOUNDARYINFO_BIND(_m, _B, UG_2D, dd_subdomain, part, _class_name, "dd_subdomain")
//#else
#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_UG(_m, _B, _class_name)
//#endif

//#if HAVE_ALBERTA
//#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_ALBERTA(_m, _B, _class_name)                                                 \
//  _DUNE_XT_GRID_BOUNDARYINFO_BIND(_m, _B, ALBERTA_2D, leaf, view, _class_name, "");                                  \
//  _DUNE_XT_GRID_BOUNDARYINFO_BIND(_m, _B, ALBERTA_2D, dd_subdomain, part, _class_name, "dd_subdomain")
//#else
#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_ALBERTA(_m, _B, _class_name)
//#endif

#define _DUNE_XT_GRID_BOUNDARYINFO_BIND_ALL(_m, _B, _class_name)                                                       \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_YASP(_m, Dune::XT::Grid::_B, _class_name);                                           \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_UG(_m, Dune::XT::Grid::_B, _class_name);                                             \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_ALBERTA(_m, Dune::XT::Grid::_B, _class_name);                                        \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_ALU(_m, Dune::XT::Grid::_B, _class_name)

#define DUNE_XT_GRID_BOUNDARYINFO_BIND(_m)                                                                             \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_ALL(_m, AllDirichletBoundaryInfo, "all_dirichlet_boundary_info");                    \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_ALL(_m, AllNeumannBoundaryInfo, "all_neumann_boundary_info");                        \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_ALL(                                                                                 \
      _m, BoundarySegmentIndexBasedBoundaryInfo, "boundary_segment_index_based_boundary_info");                        \
  _DUNE_XT_GRID_BOUNDARYINFO_BIND_ALL(_m, NormalBasedBoundaryInfo, "normal_based_boundary_info")

// end: this is what we need for the .so


#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_XT_GRID_BOUNDARYINFO_BINDINGS_HH
