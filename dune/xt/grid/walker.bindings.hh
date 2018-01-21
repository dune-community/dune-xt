// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_GRID_WALKER_BINDINGS_HH
#define DUNE_XT_GRID_WALKER_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include "dd/subdomains/grid.hh"
#include "gridprovider/provider.hh"
#include "grids.bindings.hh"
#include "layers.bindings.hh"
#include "type_traits.hh"
#include "walker.hh"

namespace Dune {
namespace XT {
namespace Grid {
namespace bindings {


template <class G, Layers layer, Backends backend>
class Walker
{
  static_assert(is_grid<G>::value, "");
  typedef typename Layer<G, layer, backend>::type GL;

public:
  typedef XT::Grid::Walker<GL> type;
  typedef pybind11::class_<type> bound_type;

private:
  template <bool is_dd = (layer == Layers::dd_subdomain || layer == Layers::dd_subdomain_boundary
                          || layer == Layers::dd_subdomain_coupling
                          || layer == Layers::dd_subdomain_oversampled),
            bool anything = true>
  struct factory_method
  {
    static void addbind(pybind11::module& m)
    {
      using namespace pybind11::literals;

      m.def(std::string("make_walker_on_" + layer_name<layer>::value() + "_" + backend_name<backend>::value()).c_str(),
            [](GridProvider<G, DD::SubdomainGrid<G>>& grid_provider, const int level_or_subdomain) {
              return type(grid_provider.template layer<layer, backend>(level_or_subdomain));
            },
            "grid_provider"_a,
            "level_or_subdomain"_a = -1);
    }
  }; // struct factory_method<true, ...>

  template <bool anything>
  struct factory_method<false, anything>
  {
    static void addbind(pybind11::module& m)
    {
      using namespace pybind11::literals;

      m.def(std::string("make_walker_on_" + layer_name<layer>::value() + "_" + backend_name<backend>::value()).c_str(),
            [](GridProvider<G>& grid_provider, const int level) {
              return type(grid_provider.template layer<layer, backend>(level));
            },
            "grid_provider"_a,
            "level"_a = -1);

      factory_method<true>::addbind(m);
    }
  }; // struct factory_method<false, ...>

public:
  static bound_type bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    // we need to add the factory methods first, since adding the class below might fail (if someone added it before)
    factory_method<>::addbind(m);

    const auto ClassName = Common::to_camel_case(
        "walker_" + grid_name<G>::value() + "_" + layer_name<layer>::value() + "_" + backend_name<backend>::value());

    bound_type c(m, ClassName.c_str());
    c.def("clear", &type::clear);
    c.def("append", [](type& self, type& other) { self.append(other); }, "other"_a, py::keep_alive<1, 2>());
    c.def("prepare", &type::prepare);
    c.def("finalize", &type::finalize);
    c.def("walk", [](type& self, const bool use_tbb) { self.walk(use_tbb); }, "use_tbb"_a = false);

    return c;
  } // ... bind(...)
}; // class Walker

} // namespace bindings
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_XT_GRID_WALKER_BINDINGS_HH
