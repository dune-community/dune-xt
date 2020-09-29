// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2019)
//   Tobias Leibner  (2020)

#ifndef PYTHON_DUNE_XT_FUNCTIONS_FUNCTION_INTERFACE_HH
#define PYTHON_DUNE_XT_FUNCTIONS_FUNCTION_INTERFACE_HH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/operators.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/gridprovider/provider.hh>

#include <dune/xt/functions/interfaces/function.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace internal {


template <class L, class R, class comb>
struct get_combined
{}; // struct get_combined

template <class L, class R>
struct get_combined<L, R, CombinationType::difference>
{
  using type = DifferenceFunction<L, R>;

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
struct get_combined<L, R, CombinationType::sum>
{
  using type = SumFunction<L, R>;

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
struct get_combined<L, R, CombinationType::product>
{
  using type = ProductFunction<L, R>;

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


} // namespace internal


template <size_t d,
          class comb,
          size_t lr,
          size_t lrC,
          size_t rr,
          size_t rrC,
          class C = typename internal::
              get_combined<FunctionInterface<d, lr, lrC, double>, FunctionInterface<d, rr, rrC, double>, comb>::type>
pybind11::class_<C, FunctionInterface<d, C::range_dim, C::range_dim_cols, double>>
bind_combined_Function(pybind11::module& m)
{
  namespace py = pybind11;

  using R = double;
  using Left = FunctionInterface<d, lr, lrC, R>;
  using Right = FunctionInterface<d, rr, rrC, R>;
  static constexpr size_t r = C::range_dim;
  static constexpr size_t rC = C::range_dim_cols;
  const std::string id = internal::get_combined<Left, Right, comb>::id();
  const std::string op = internal::get_combined<Left, Right, comb>::doc();
  const std::string class_name =
      id + "__" + Common::to_string(d) + "d_to_" + Common::to_string(r) + "x" + Common::to_string(rC);
  const std::string doc = class_name + " (as a " + op + " of functions of dimensions " + Common::to_string(lr) + "x"
                          + Common::to_string(lrC) + " and " + Common::to_string(rr) + "x" + Common::to_string(rrC)
                          + ")";

  py::class_<C, FunctionInterface<d, r, rC, R>> c(m, std::string(class_name).c_str(), doc.c_str());

  c.def_property_readonly("static_id", [](const C& /*self*/) { return C::static_id(); });

  return c;
} // ... bind_combined_Function(...)


template <size_t d, class comb, size_t r, size_t rC, size_t oR, size_t orC, class C>
void addbind_FunctionInterface_combined_op(C& c)
{
  namespace py = pybind11;

  using S = FunctionInterface<d, r, rC, double>;
  using O = FunctionInterface<d, oR, orC, double>;

  c.def(
      internal::get_combined<S, O, comb>::op().c_str(),
      [](const S& self, const O& other) { return internal::get_combined<S, O, comb>::call(self, other); },
      py::is_operator());
} // ... addbind_FunctionInterface_combined_op(...)


template <size_t d, size_t r, size_t rC>
pybind11::class_<FunctionInterface<d, r, rC, double>> bind_FunctionInterface(pybind11::module& m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  using C = FunctionInterface<d, r, rC, double>;

  py::class_<C> c(m,
                  std::string("FunctionInterface__" + Common::to_string(d) + "d_to_" + Common::to_string(r) + "x"
                              + Common::to_string(rC))
                      .c_str(),
                  std::string("FunctionInterface__" + Common::to_string(d) + "d_to_" + Common::to_string(r) + "x"
                              + Common::to_string(rC))
                      .c_str());

  c.def_property_readonly("dim_domain", [](const C& /*self*/) { return size_t(d); });
  if (rC == 1)
    c.def_property_readonly("dim_range", [](const C& /*self*/) { return size_t(r); });
  else
    c.def_property_readonly("dim_range", [](const C& /*self*/) { return std::make_pair(size_t(r), size_t(rC)); });
  c.def_property_readonly("static_id", [](const C& /*self*/) { return C::static_id(); });
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

  // internal::Divergence<G>::addbind(m, c);

  return c;
} // ... bind_FunctionInterface(...)


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // PYTHON_DUNE_XT_FUNCTIONS_FUNCTION_INTERFACE_HH
