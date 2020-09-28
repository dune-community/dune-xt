// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License
// (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_XT_FUNCTIONS_INTERFACES_GRID_FUNCTION_HH
#define PYTHON_DUNE_XT_FUNCTIONS_INTERFACES_GRID_FUNCTION_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/visualization.hh>
#include <dune/xt/grid/gridprovider/provider.hh>

#include <python/dune/xt/common/parameter.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace bindings {


template <class G, class E, size_t r = 1, size_t rC = 1, class R = double>
class GridFunctionInterface
{
  using GP = XT::Grid::GridProvider<G>;
  static const constexpr size_t d = G::dimension;

  template <bool vector = (r != 1 && rC == 1), bool matrix = (rC != 1), bool anything = false>
  struct product_helper // <true, false, ...>
  {
    template <class T, typename... options>
    static void addbind(pybind11::class_<T, options...>& c)
    {
      namespace py = pybind11;

      c.def(
          "__mul__",
          [](const T& self, const Functions::GridFunctionInterface<E, r, 1, R>& other) {
            return std::make_unique<decltype(self * other)>(self * other);
          },
          py::is_operator());
    }
  };

  template <bool anything>
  struct product_helper<false, true, anything>
  {
    template <class T, typename... options>
    static void addbind(pybind11::class_<T, options...>& c)
    {
      namespace py = pybind11;

      c.def(
          "__mul__",
          [](const T& self, const Functions::GridFunctionInterface<E, rC, 1, R>& other) {
            return std::make_unique<decltype(self * other)>(self * other);
          },
          py::is_operator());
      c.def(
          "__mul__",
          [](const T& self, const Functions::GridFunctionInterface<E, rC, 2, R>& other) {
            return std::make_unique<decltype(self * other)>(self * other);
          },
          py::is_operator());
      c.def(
          "__mul__",
          [](const T& self, const Functions::GridFunctionInterface<E, rC, 3, R>& other) {
            return std::make_unique<decltype(self * other)>(self * other);
          },
          py::is_operator());
    }
  };

  template <bool scalar = (r == 1 && rC == 1), bool anything = false>
  struct fraction_helper // <true, ...>
  {
    template <class T, typename... options>
    static void addbind(pybind11::class_<T, options...>& c)
    {
      namespace py = pybind11;

      c.def(
          "__truediv__",
          [](const T& self, const type& other) { return std::make_unique<decltype(other / self)>(other / self); },
          py::is_operator());
    }
  };

  template <bool anything>
  struct fraction_helper<false, anything>
  {
    template <class T, typename... options>
    static void addbind(pybind11::class_<T, options...>& /*c*/)
    {}
  };

public:
  using type = Functions::GridFunctionInterface<E, r, rC, R>;
  using bound_type = pybind11::class_<type>;

  static std::string class_name(const std::string& grid_id, const std::string& layer_id, const std::string& class_id)
  {
    std::string ret = class_id;
    ret += "_" + grid_id;
    if (!layer_id.empty())
      ret += "_" + layer_id;
    ret += "_to_" + Common::to_string(r);
    if (rC > 1)
      ret += "x" + Common::to_string(rC);
    ret += "d";
    if (!std::is_same<R, double>::value)
      ret += "_" + Common::Typename<R>::value(/*fail_wo_typeid=*/true);
    return ret;
  } // ... class_name(...)

  template <class T, typename... options>
  static void addbind_methods(pybind11::class_<T, options...>& c)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    // our methods
    c.def(
        "visualize",
        [](const T& self, const GP& grid_provider, const std::string& filename, const bool subsampling) {
          Functions::visualize(self, grid_provider.leaf_view(), filename, subsampling);
        },
        "grid"_a,
        "filename"_a,
        "subsampling"_a = true);
    // our operators
    c.def(
        "__add__",
        [](const T& self, const type& other) { return std::make_unique<decltype(self + other)>(self + other); },
        py::is_operator());
    c.def(
        "__sub__",
        [](const T& self, const type& other) { return std::make_unique<decltype(self - other)>(self - other); },
        py::is_operator());
    // we can always multiply with a scalar from the right ...
    c.def(
        "__mul__",
        [](const T& self, const Functions::GridFunctionInterface<E, 1, 1, R>& other) {
          return std::make_unique<decltype(self * other)>(self * other);
        },
        py::is_operator());
    // .. and with lots of other dims
    product_helper<>::addbind(c);
    fraction_helper<>::addbind(c);

    if constexpr (r == 1 && rC == 1)
      c.def(
          "__pow__",
          [](const T& self) { return std::make_unique<decltype(self * self)>(self * self); },
          py::is_operator());

    // ParametricInterface methods
    c.def(
        "parse_parameter", [](const T& self, const Common::Parameter& mu) { return self.parse_parameter(mu); }, "mu"_a);
  } // ... addbind_methods(...)

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id,
                         const std::string& layer_id = "",
                         const std::string& class_id = "grid_function_interface")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = Common::to_camel_case(class_name(grid_id, layer_id, class_id));
    bound_type c(m, ClassName.c_str(), Common::to_camel_case(class_id).c_str());

    // our properties
    c.def_property_readonly("dim_domain", [](type&) { return size_t(d); });
    if (rC == 1)
      c.def_property_readonly("dim_range", [](type&) { return size_t(r); });
    else
      c.def_property_readonly("dim_range", [](type&) { return std::make_pair(size_t(r), size_t(rC)); });
    c.def_property_readonly("name", [](type& self) { return self.name(); });
    // ParametricInterface properties
    c.def_property_readonly("is_parametric", [](type& self) { return self.is_parametric(); });
    c.def_property_readonly("parameter_type", [](type& self) { return self.parameter_type(); });

    addbind_methods(c);

    return c;
  }
}; // class GridFunctionInterface


} // namespace bindings
} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // PYTHON_DUNE_XT_FUNCTIONS_INTERFACES_GRID_FUNCTION_HH
