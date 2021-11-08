// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019 - 2020)
//   René Fritze     (2019 - 2020)
//   Tobias Leibner  (2020)

#ifndef PYTHON_DUNE_XT_FUNCTIONS_FUNCTION_INTERFACE_HH
#define PYTHON_DUNE_XT_FUNCTIONS_FUNCTION_INTERFACE_HH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/numpy.h>
#include <dune/pybindxi/operators.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/la/container/common.hh>
#include <dune/xt/grid/gridprovider/provider.hh>

#include <dune/xt/functions/interfaces/function.hh>

namespace Dune::XT::Functions {
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

  c.def(
      "evaluate",
      [](const C& self, const typename C::DomainType& x, const XT::Common::Parameter& mu) {
        return self.evaluate(x, mu);
      },
      "x"_a,
      "mu"_a = XT::Common::Parameter());
  if constexpr (rC == 1) {
    c.def(
        "evaluate",
        [](const C& self, py::array_t<double> list_of_points, const XT::Common::Parameter& mu) {
          DUNE_THROW_IF(list_of_points.ndim() != 2,
                        XT::Common::Exceptions::shapes_do_not_match,
                        "list_of_points.ndim() = " << list_of_points.ndim() << " (has to be 2)!");
          DUNE_THROW_IF(list_of_points.shape(1) != d,
                        XT::Common::Exceptions::shapes_do_not_match,
                        "list_of_points.shape(1) = " << list_of_points.shape(1) << " (has to be " << size_t(d) << ")!");
          const auto& access_to_list_of_points = list_of_points.unchecked<2>();
          const size_t num_points = access_to_list_of_points.shape(0);
          py::array_t<double> values(/*shape=*/{num_points, r});
          auto access_to_values = values.mutable_unchecked<2>();
          typename C::DomainType point;
          for (size_t ii = 0; ii < num_points; ++ii) {
            for (size_t dd = 0; dd < d; ++dd)
              point[dd] = access_to_list_of_points(ii, dd);
            const auto value = self.evaluate(point, mu);
            for (size_t rr = 0; rr < r; ++rr)
              access_to_values(ii, rr) = value[rr];
          }
          return values;
        },
        "x"_a,
        "mu"_a = XT::Common::Parameter());
  }

  // internal::Divergence<G>::addbind(m, c);

  return c;
} // ... bind_FunctionInterface(...)


} // namespace Dune::XT::Functions

#endif // PYTHON_DUNE_XT_FUNCTIONS_FUNCTION_INTERFACE_HH
