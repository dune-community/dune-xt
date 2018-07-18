// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_FUNCTIONS_INTERFACE_PBH
#define DUNE_XT_FUNCTIONS_INTERFACE_PBH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/operators.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/gridprovider/provider.hh>

#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/interfaces/function.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace internal {


template <class L, class R, internal::Combination comb>
struct get_combined
{
}; // struct get_combined

template <class L, class R>
struct get_combined<L, R, internal::Combination::difference>
{
  typedef DifferenceFunction<L, R> type;

  static std::string id()
  {
    return "DifferenceFunction";
  }

  static std::string doc()
  {
    return "difference";
  }

  static std::string op()
  {
    return "__sub__";
  }

  static auto call(const L& l, const R& r) -> decltype(l - r)
  {
    return l - r;
  }
}; // struct get_combined

template <class L, class R>
struct get_combined<L, R, internal::Combination::sum>
{
  typedef SumFunction<L, R> type;

  static std::string id()
  {
    return "SumFunction";
  }

  static std::string doc()
  {
    return "sum";
  }

  static std::string op()
  {
    return "__add__";
  }

  static auto call(const L& l, const R& r) -> decltype(l + r)
  {
    return l + r;
  }
}; // struct get_combined

template <class L, class R>
struct get_combined<L, R, internal::Combination::product>
{
  typedef ProductFunction<L, R> type;

  static std::string id()
  {
    return "ProductFunction";
  }

  static std::string doc()
  {
    return "product";
  }

  static std::string op()
  {
    return "__mul__";
  }

  static auto call(const L& l, const R& r) -> decltype(l * r)
  {
    return l * r;
  }
}; // struct get_combined


/**
 * grid_combined
 */

template <class L, class R, internal::Combinations comb>
struct get_grid_combined
{
}; // struct get_grid_combined

template <class L, class R>
struct get_grid_combined<L, R, internal::Combinations::difference>
{
  typedef DifferenceGridFunction<L, R> type;

  static std::string id()
  {
    return "DifferenceGridFunction";
  }

  static std::string doc()
  {
    return "difference";
  }

  static std::string op()
  {
    return "__sub__";
  }

  static auto call(const L& l, const R& r) -> decltype(l - r)
  {
    return l - r;
  }
}; // struct get_grid_combined

template <class L, class R>
struct get_grid_combined<L, R, internal::Combinations::sum>
{
  typedef SumGridFunction<L, R> type;

  static std::string id()
  {
    return "SumGridFunction";
  }

  static std::string doc()
  {
    return "sum";
  }

  static std::string op()
  {
    return "__add__";
  }

  static auto call(const L& l, const R& r) -> decltype(l + r)
  {
    return l + r;
  }
}; // struct get_grid_combined

template <class L, class R>
struct get_grid_combined<L, R, internal::Combinations::product>
{
  typedef ProductGridFunction<L, R> type;

  static std::string id()
  {
    return "ProductGridFunction";
  }

  static std::string doc()
  {
    return "product";
  }

  static std::string op()
  {
    return "__mul__";
  }

  static auto call(const L& l, const R& r) -> decltype(l * r)
  {
    return l * r;
  }
}; // struct get_grid_combined

// template <class G>
// struct Divergence
//{
//  template <size_t d, size_t r, size_t rC, bool dims_match = (d == r) && (rC == 1)>
//  struct helper
//  {
//    template <class E, class R>
//    static void addbind(pybind11::module& m, pybind11::class_<GridFunctionInterface<E, d, 1, R>>& c)
//    {
//      namespace py = pybind11;
//      using namespace pybind11::literals;
//      using Common::to_string;

//      try { // guard since we might not be the first to do bind this combination
//        py::class_<DivergenceFunction<GridFunctionInterface<E, d, 1, R>>,
//                   GridFunctionInterface<E, d, 1, R>>(
//            m,
//            Common::to_camel_case(
//                "divergence_of_function_from_" + Grid::bindings::grid_name<G>::value() + "_to_" + to_string(r) + "x"
//                + to_string(rC))
//                .c_str(),
//            "DivergenceFunction");
//      } catch (std::runtime_error&) {
//      }

//      c.def("divergence",
//            [](const GridFunctionInterface<E, d, 1, R>& self, const std::string& name) {
//              return new DivergenceFunction<GridFunctionInterface<E, d, 1, R>>(self, name);
//            },
//            "name"_a = "",
//            py::keep_alive<0, 1>());
//    } // ... addbind(...)
//  }; // struct helper<..., true>

//  template <size_t d, size_t r, size_t rC>
//  struct helper<d, r, rC, false>
//  {
//    template <class E, class R>
//    static void addbind(pybind11::module& /*m*/, pybind11::class_<GridFunctionInterface<E, r, rC, R>>& /*c*/)
//    {
//    }
//  }; // struct helper<..., false>

//  template <class E, size_t d, class R, size_t r, size_t rC>
//  static void addbind(pybind11::module& m, pybind11::class_<GridFunctionInterface<E, r, rC, R>>& c)
//  {
//    helper<d, r, rC>::addbind(m, c);
//  } // ... addbind(...)
//}; // struct Divergence


} // namespace internal


template <class G, size_t r, size_t rC>
pybind11::class_<GridFunctionInterface<typename G::template Codim<0>::Entity, r, rC, double>>
bind_GridFunctionInterface(pybind11::module& m, const std::string& grid_id)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  typedef GridFunctionInterface<typename G::template Codim<0>::Entity, r, rC, double> C;

  py::class_<C> c(
      m,
      std::string("GridFunctionInterface__" + grid_id + "_to_" + Common::to_string(r) + "x" + Common::to_string(rC))
          .c_str(),
      std::string("GridFunctionInterface__" + grid_id + "_to_" + Common::to_string(r) + "x" + Common::to_string(rC))
          .c_str());

  c.def_property_readonly("static_id", [](const C& /*self*/) { return C::static_id(); });
  c.def_property_readonly("type", [](const C& self) { return self.type(); });
  c.def_property_readonly("name", [](const C& self) { return self.name(); });

  c.def("visualize",
        [](const C& self,
           const Grid::GridProvider<G, Grid::none_t>& grid_provider,
           const std::string& layer,
           const ssize_t lvl,
           const std::string& path,
           const bool subsampling) {
          const auto level = XT::Common::numeric_cast<int>(lvl);
          if (layer == "leaf")
            self.visualize(grid_provider.leaf_view(), path, subsampling);
          else if (layer == "level")
            self.visualize(grid_provider.template layer<XT::Grid::Layers::level, XT::Grid::Backends::view>(level),
                           path,
                           subsampling);
          else
            DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
                       "Given layer has to be one of ('leaf', 'level'), is '" << layer << "'!");
        },
        "grid_provider"_a,
        "layer"_a = "leaf",
        "level"_a = -1,
        "path"_a,
        "subsampling"_a = true);
  c.def("visualize",
        [](const C& self,
           const Grid::GridProvider<G, Grid::DD::SubdomainGrid<G>>& dd_grid_provider,
           const std::string& layer,
           const ssize_t lvl_or_sbdmn,
           const std::string& path,
           const bool subsampling) {
          const auto level_or_subdomain = XT::Common::numeric_cast<int>(lvl_or_sbdmn);
          if (layer == "leaf")
            self.visualize(dd_grid_provider.leaf_view(), path, subsampling);
          else if (layer == "level")
            self.visualize(
                dd_grid_provider.template layer<XT::Grid::Layers::level, XT::Grid::Backends::view>(level_or_subdomain),
                path,
                subsampling);
          else if (layer == "dd_subdomain")
            self.visualize(dd_grid_provider.template layer<XT::Grid::Layers::dd_subdomain, XT::Grid::Backends::view>(
                               level_or_subdomain),
                           path,
                           subsampling);
          else if (layer == "dd_subdomain_oversampled")
            self.visualize(
                dd_grid_provider.template layer<XT::Grid::Layers::dd_subdomain_oversampled, XT::Grid::Backends::view>(
                    level_or_subdomain),
                path,
                subsampling);
          else
            DUNE_THROW(
                XT::Common::Exceptions::wrong_input_given,
                "Given layer has to be one of ('leaf', 'level', 'dd_subdomain', 'dd_subdomain_oversampled'), is '"
                    << layer
                    << "'!");
        },
        "dd_grid_provider"_a,
        "layer"_a = "leaf",
        "level_or_subdomain"_a = -1,
        "path"_a,
        "subsampling"_a = true);

  // internal::Divergence<G>::addbind(m, c);

  return c;
} // ... bind_GridFunctionInterface(...)


/**
 * \note We would like to drop the d template paremter and use either of
\code
static const           size_t d = G::dimension;
static const constexpr size_t d = G::dimension;
\endcode
 *       but this triggers a bug in gcc-4.9, see e.g.: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=59937
 */
template <class G, size_t d, internal::Combinations comb, size_t lr, size_t lrC, size_t rr, size_t rrC>
pybind11::class_<
    typename internal::get_grid_combined<GridFunctionInterface<typename G::template Codim<0>::Entity, lr, lrC, double>,
                                    GridFunctionInterface<typename G::template Codim<0>::Entity, rr, rrC, double>,
                                    comb>::type>
bind_combined_GridFunction(pybind11::module& m, const std::string& grid_id)
{
  namespace py = pybind11;

  typedef typename G::template Codim<0>::Entity E;
  typedef double R;
  typedef GridFunctionInterface<E, lr, lrC, R> Left;
  typedef GridFunctionInterface<E, rr, rrC, R> Right;
  typedef typename internal::get_grid_combined<Left, Right, comb>::type C;
  static const size_t r = C::range_dim;
  static const size_t rC = C::range_dim_cols;
  const std::string id = internal::get_grid_combined<Left, Right, comb>::id();
  const std::string op = internal::get_grid_combined<Left, Right, comb>::doc();
  const std::string class_name = id + "__" + grid_id + "_to_" + Common::to_string(r) + "x" + Common::to_string(rC);
  const std::string doc = class_name + " (as a " + op + " of functions of dimensions " + Common::to_string(lr) + "x"
                          + Common::to_string(lrC) + " and " + Common::to_string(rr) + "x" + Common::to_string(rrC)
                          + ")";

  py::class_<C, GridFunctionInterface<E, r, rC, R>> c(m, std::string(class_name).c_str(), doc.c_str());

  c.def_property_readonly("static_id", [](const C& /*self*/) { return C::static_id(); });

  return c;
} // ... bind_combined_GridFunction(...)


template <size_t d, internal::Combination comb, size_t lr, size_t lrC, size_t rr, size_t rrC>
pybind11::class_<
    typename internal::get_combined<FunctionInterface<d, lr, lrC, double>,
                                    FunctionInterface<d, rr, rrC, double>,
                                    comb>::type>
bind_combined_Function(pybind11::module& m)
{
  namespace py = pybind11;

  typedef double R;
  typedef FunctionInterface<d, lr, lrC, R> Left;
  typedef FunctionInterface<d, rr, rrC, R> Right;
  typedef typename internal::get_combined<Left, Right, comb>::type C;
  static const size_t r = C::range_dim;
  static const size_t rC = C::range_dim_cols;
  const std::string id = internal::get_combined<Left, Right, comb>::id();
  const std::string op = internal::get_combined<Left, Right, comb>::doc();
  const std::string class_name = id + "_to_" + Common::to_string(r) + "x" + Common::to_string(rC);
  const std::string doc = class_name + " (as a " + op + " of functions of dimensions " + Common::to_string(lr) + "x"
                          + Common::to_string(lrC) + " and " + Common::to_string(rr) + "x" + Common::to_string(rrC)
                          + ")";

  py::class_<C, FunctionInterface<d, r, rC, R>> c(m, std::string(class_name).c_str(), doc.c_str());

  c.def_property_readonly("static_id", [](const C& /*self*/) { return C::static_id(); });

  return c;
} // ... bind_combined_Function(...)


/**
 * \note We would like to drop the d template paremter and use either of
\code
static const           size_t d = G::dimension;
static const constexpr size_t d = G::dimension;
\endcode
 *       but this triggers a bug in gcc-4.9, see e.g.: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=59937
 */
template <class G, size_t d, internal::Combinations comb, size_t r, size_t rC, size_t oR, size_t orC, class C>
void addbind_GridFunctionInterface_combined_op(C& c)
{
  namespace py = pybind11;

  typedef typename G::template Codim<0>::Entity E;
  typedef GridFunctionInterface<E, r, rC, double> S;
  typedef GridFunctionInterface<E, oR, orC, double> O;

  c.def(internal::get_grid_combined<S, O, comb>::op().c_str(),
        [](const S& self, const O& other) { return internal::get_grid_combined<S, O, comb>::call(self, other); },
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>());
} // ... addbind_GridFunctionInterface_combined_op(...)


template <size_t d, internal::Combination comb, size_t r, size_t rC, size_t oR, size_t orC, class C>
void addbind_FunctionInterface_combined_op(C& c)
{
  namespace py = pybind11;

  typedef FunctionInterface<d, r, rC, double> S;
  typedef FunctionInterface<d, oR, orC, double> O;

  c.def(internal::get_combined<S, O, comb>::op().c_str(),
        [](const S& self, const O& other) { return internal::get_combined<S, O, comb>::call(self, other); },
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>());
} // ... addbind_FunctionInterface_combined_op(...)


template <size_t d, size_t r, size_t rC>
pybind11::class_<FunctionInterface<d, r, rC, double>> bind_FunctionInterface(pybind11::module& m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  typedef FunctionInterface<d, r, rC, double> C;

  py::class_<C> c(m,
                  std::string("FunctionInterface__" + Common::to_string(d) + "d_to_" + Common::to_string(r) + "x"
                              + Common::to_string(rC))
                      .c_str(),
                  std::string("FunctionInterface__" + Common::to_string(d) + "d_to_" + Common::to_string(r) + "x"
                              + Common::to_string(rC))
                      .c_str());

  c.def_property_readonly("static_id", [](const C& /*self*/) { return C::static_id(); });
  c.def_property_readonly("type", [](const C& self) { return self.type(); });
  c.def_property_readonly("name", [](const C& self) { return self.name(); });

  //  c.def("visualize",
  //        [](const C& self,
  //           const Grid::GridProvider<G>& grid_provider,
  //           const std::string& layer,
  //           const ssize_t lvl,
  //           const std::string& path,
  //           const bool subsampling) {
  //          const auto level = XT::Common::numeric_cast<int>(lvl);
  //          if (layer == "leaf")
  //            self.visualize(grid_provider.leaf_view(), path, subsampling);
  //          else if (layer == "level")
  //            self.visualize(grid_provider.template layer<XT::Grid::Layers::level, XT::Grid::Backends::view>(level),
  //                           path,
  //                           subsampling);
  //          else
  //            DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
  //                       "Given layer has to be one of ('leaf', 'level'), is '" << layer << "'!");
  //        },
  //        "grid_provider"_a,
  //        "layer"_a = "leaf",
  //        "level"_a = -1,
  //        "path"_a,
  //        "subsampling"_a = true);
  //  c.def("visualize",
  //        [](const C& self,
  //           const Grid::GridProvider<G, Grid::DD::SubdomainGrid<G>>& dd_grid_provider,
  //           const std::string& layer,
  //           const ssize_t lvl_or_sbdmn,
  //           const std::string& path,
  //           const bool subsampling) {
  //          const auto level_or_subdomain = XT::Common::numeric_cast<int>(lvl_or_sbdmn);
  //          if (layer == "leaf")
  //            self.visualize(dd_grid_provider.leaf_view(), path, subsampling);
  //          else if (layer == "level")
  //            self.visualize(
  //                dd_grid_provider.template layer<XT::Grid::Layers::level,
  //                XT::Grid::Backends::view>(level_or_subdomain),
  //                path,
  //                subsampling);
  //          else if (layer == "dd_subdomain")
  //            self.visualize(dd_grid_provider.template layer<XT::Grid::Layers::dd_subdomain,
  //            XT::Grid::Backends::view>(
  //                               level_or_subdomain),
  //                           path,
  //                           subsampling);
  //          else if (layer == "dd_subdomain_oversampled")
  //            self.visualize(
  //                dd_grid_provider.template layer<XT::Grid::Layers::dd_subdomain_oversampled,
  //                XT::Grid::Backends::view>(
  //                    level_or_subdomain),
  //                path,
  //                subsampling);
  //          else
  //            DUNE_THROW(
  //                XT::Common::Exceptions::wrong_input_given,
  //                "Given layer has to be one of ('leaf', 'level', 'dd_subdomain', 'dd_subdomain_oversampled'), is '"
  //                    << layer
  //                    << "'!");
  //        },
  //        "dd_grid_provider"_a,
  //        "layer"_a = "leaf",
  //        "level_or_subdomain"_a = -1,
  //        "path"_a,
  //        "subsampling"_a = true);

  // internal::Divergence<G>::addbind(m, c);

  return c;
} // ... bind_GridFunctionInterface(...)

} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_INTERFACE_PBH
