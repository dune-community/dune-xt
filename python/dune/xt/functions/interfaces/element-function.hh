// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License
// (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_XT_FUNCTIONS_INTERFACES_ELEMENT_FUNCTION_HH
#define PYTHON_DUNE_XT_FUNCTIONS_INTERFACES_ELEMENT_FUNCTION_HH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/numpy.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/functions/interfaces/element-functions.hh>

#include <python/dune/xt/common/parameter.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace bindings {


template <class E, size_t r = 1, size_t rC = 1, class R = double>
class ElementFunctionSetInterface
{
  static const constexpr size_t d = E::dimension;

  // bind combined-element-functions.hh first to enable
#if 0
  template <bool vector = (r != 1 && rC == 1), bool matrix = (rC != 1), bool anything = false>
  struct product_helper // <true, false, ...>
  {
    template <class T, typename... options>
    static void addbind(pybind11::class_<T, options...>& c)
    {
      namespace py = pybind11;

      c.def(
          "__mul__",
          [](const T& self, const Functions::ElementFunctionSetInterface<E, r, 1, R>& other) {
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
          [](const T& self, const Functions::ElementFunctionSetInterface<E, rC, 1, R>& other) {
            return std::make_unique<decltype(self * other)>(self * other);
          },
          py::is_operator());
      c.def(
          "__mul__",
          [](const T& self, const Functions::ElementFunctionSetInterface<E, rC, 2, R>& other) {
            return std::make_unique<decltype(self * other)>(self * other);
          },
          py::is_operator());
      c.def(
          "__mul__",
          [](const T& self, const Functions::ElementFunctionSetInterface<E, rC, 3, R>& other) {
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
#endif // 0

public:
  using type = Functions::ElementFunctionSetInterface<E, r, rC, R>;
  using bound_type = pybind11::class_<type, Common::ParametricInterface>;

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
        "size", [](type& self, const Common::Parameter& mu) { return self.size(mu); }, "mu"_a = Common::Parameter());
    c.def(
        "max_size",
        [](type& self, const Common::Parameter& mu) { return self.max_size(mu); },
        "mu"_a = Common::Parameter());
    c.def(
        "order", [](type& self, const Common::Parameter& mu) { return self.order(mu); }, "mu"_a = Common::Parameter());
    c.def(
        "evaluate_set",
        [](type& self, py::array_t<double> points_in_reference_element, const Common::Parameter& mu) {
          py::array_t<R> result;
          if (points_in_reference_element.ndim() == 1) {
            // a single point
            DUNE_THROW_IF(points_in_reference_element.shape(0) != d,
                          XT::Common::Exceptions::shapes_do_not_match,
                          "points_in_reference_element.shape(0) = " << points_in_reference_element.shape(0)
                                                                    << " (has to be " << size_t(d) << ")!");
            const auto& access_to_points_in_reference_element = points_in_reference_element.unchecked<1>();
            XT::Common::FieldVector<double, d> point_in_reference_element_dune;
            for (size_t dd = 0; dd < d; ++dd)
              point_in_reference_element_dune[dd] = access_to_points_in_reference_element(dd);
            auto values = self.evaluate_set(point_in_reference_element_dune, mu);
            const size_t num_functions_in_set = values.size();
            if constexpr (rC == 1) {
              result = py::array_t<double>(/*shape=*/{num_functions_in_set, r});
              auto access_to_result = result.template mutable_unchecked<2>();
              for (size_t ff = 0; ff < num_functions_in_set; ++ff)
                for (size_t ii = 0; ii < r; ++ii)
                  access_to_result(ff, ii) = values[ff][ii];
            } else {
              result = py::array_t<double>(/*shape=*/{num_functions_in_set, r, rC});
              auto access_to_result = result.template mutable_unchecked<3>();
              for (size_t ff = 0; ff < num_functions_in_set; ++ff)
                for (size_t ii = 0; ii < r; ++ii)
                  for (size_t jj = 0; jj < rC; ++jj)
                    access_to_result(ff, ii, jj) = values[ff][ii][jj];
            }
          } else if (points_in_reference_element.ndim() == 2) {
            // a list of points
            DUNE_THROW_IF(points_in_reference_element.shape(1) != d,
                          XT::Common::Exceptions::shapes_do_not_match,
                          "points_in_reference_element.shape(1) = " << points_in_reference_element.shape(1)
                                                                    << " (has to be " << size_t(d) << ")!");
            const auto& access_to_points_in_reference_element = points_in_reference_element.unchecked<2>();
            const size_t num_points = access_to_points_in_reference_element.shape(0);
            const size_t num_functions_in_set = self.size(mu);
            XT::Common::FieldVector<double, d> point_in_reference_element_dune;
            std::vector<typename type::RangeType> values(num_functions_in_set);
            if constexpr (rC == 1) {
              result = py::array_t<double>(/*shape=*/{num_points, num_functions_in_set, r});
              auto access_to_result = result.template mutable_unchecked<3>();
              for (size_t pp = 0; pp < num_points; ++pp) {
                for (size_t ii = 0; ii < d; ++ii)
                  point_in_reference_element_dune[ii] = access_to_points_in_reference_element(pp, ii);
                self.evaluate(point_in_reference_element_dune, values);
                for (size_t ff = 0; ff < num_functions_in_set; ++ff)
                  for (size_t ii = 0; ii < r; ++ii)
                    access_to_result(pp, ff, ii) = values[ff][ii];
              }
            } else {
              result = py::array_t<double>(/*shape=*/{num_points, num_functions_in_set, r, rC});
              auto access_to_result = result.template mutable_unchecked<4>();
              for (size_t pp = 0; pp < num_points; ++pp) {
                for (size_t ii = 0; ii < d; ++ii)
                  point_in_reference_element_dune[ii] = access_to_points_in_reference_element(pp, ii);
                self.evaluate(point_in_reference_element_dune, values);
                for (size_t ff = 0; ff < num_functions_in_set; ++ff)
                  for (size_t ii = 0; ii < r; ++ii)
                    for (size_t jj = 0; jj < rC; ++jj)
                      access_to_result(pp, ff, ii, jj) = values[ff][ii][jj];
              }
            }
          } else
            DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                       "points_in_reference_element.ndim() = "
                           << points_in_reference_element.ndim()
                           << " (has to be 1 for a single point, 2 for a list of points)!");
          return result;
        },
        "points_in_reference_element"_a,
        "mu"_a = Common::Parameter());
    c.def(
        "jacobians_of_set",
        [](type& self, py::array_t<double> points_in_reference_element, const Common::Parameter& mu) {
          py::array_t<R> result;
          if (points_in_reference_element.ndim() == 1) {
            // a single point
            DUNE_THROW_IF(points_in_reference_element.shape(0) != d,
                          XT::Common::Exceptions::shapes_do_not_match,
                          "points_in_reference_element.shape(0) = " << points_in_reference_element.shape(0)
                                                                    << " (has to be " << size_t(d) << ")!");
            const auto& access_to_points_in_reference_element = points_in_reference_element.unchecked<1>();
            XT::Common::FieldVector<double, d> point_in_reference_element_dune;
            for (size_t dd = 0; dd < d; ++dd)
              point_in_reference_element_dune[dd] = access_to_points_in_reference_element(dd);
            auto jacobians = self.jacobians_of_set(point_in_reference_element_dune, mu);
            const size_t num_functions_in_set = jacobians.size();
            if constexpr (rC == 1) {
              result = py::array_t<double>(/*shape=*/{num_functions_in_set, r, d});
              auto access_to_result = result.template mutable_unchecked<3>();
              for (size_t ff = 0; ff < num_functions_in_set; ++ff)
                for (size_t ii = 0; ii < r; ++ii)
                  for (size_t jj = 0; jj < d; ++jj)
                    access_to_result(ff, ii, jj) = jacobians[ff][ii][jj];
            } else {
              result = py::array_t<double>(/*shape=*/{num_functions_in_set, r, rC, d});
              auto access_to_result = result.template mutable_unchecked<4>();
              for (size_t ff = 0; ff < num_functions_in_set; ++ff)
                for (size_t ii = 0; ii < r; ++ii)
                  for (size_t jj = 0; jj < rC; ++jj)
                    for (size_t kk = 0; kk < d; ++kk)
                      access_to_result(ff, ii, jj, kk) = jacobians[ff][ii][jj][kk];
            }
          } else if (points_in_reference_element.ndim() == 2) {
            // a list of points
            DUNE_THROW_IF(points_in_reference_element.shape(1) != d,
                          XT::Common::Exceptions::shapes_do_not_match,
                          "points_in_reference_element.shape(1) = " << points_in_reference_element.shape(1)
                                                                    << " (has to be " << size_t(d) << ")!");
            const auto& access_to_points_in_reference_element = points_in_reference_element.unchecked<2>();
            const size_t num_points = access_to_points_in_reference_element.shape(0);
            const size_t num_functions_in_set = self.size(mu);
            XT::Common::FieldVector<double, d> point_in_reference_element_dune;
            std::vector<typename type::DerivativeRangeType> jacobians(num_functions_in_set);
            if constexpr (rC == 1) {
              result = py::array_t<double>(/*shape=*/{num_points, num_functions_in_set, r, d});
              auto access_to_result = result.template mutable_unchecked<4>();
              for (size_t pp = 0; pp < num_points; ++pp) {
                for (size_t ii = 0; ii < d; ++ii)
                  point_in_reference_element_dune[ii] = access_to_points_in_reference_element(pp, ii);
                self.jacobians(point_in_reference_element_dune, jacobians);
                for (size_t ff = 0; ff < num_functions_in_set; ++ff)
                  for (size_t ii = 0; ii < r; ++ii)
                    for (size_t jj = 0; jj < d; ++jj)
                      access_to_result(pp, ff, ii, jj) = jacobians[ff][ii][jj];
              }
            } else {
              result = py::array_t<double>(/*shape=*/{num_points, num_functions_in_set, r, rC, d});
              auto access_to_result = result.template mutable_unchecked<5>();
              for (size_t pp = 0; pp < num_points; ++pp) {
                for (size_t ii = 0; ii < d; ++ii)
                  point_in_reference_element_dune[ii] = access_to_points_in_reference_element(pp, ii);
                self.jacobians(point_in_reference_element_dune, jacobians);
                for (size_t ff = 0; ff < num_functions_in_set; ++ff)
                  for (size_t ii = 0; ii < r; ++ii)
                    for (size_t jj = 0; jj < rC; ++jj)
                      for (size_t kk = 0; kk < d; ++kk)
                        access_to_result(pp, ff, ii, jj, kk) = jacobians[ff][ii][jj][kk];
              }
            }
          } else
            DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                       "points_in_reference_element.ndim() = "
                           << points_in_reference_element.ndim()
                           << " (has to be 1 for a single point, 2 for a list of points)!");
          return result;
        },
        "points_in_reference_element"_a,
        "mu"_a = Common::Parameter());
    // our operators
    // bind combined-element-functions.hh first to enable
#if 0
    c.def(
        "__add__",
        [](const T& self, const type& other) { return std::make_unique<decltype(self + other)>(self + other); },
        py::is_operator());
    c.def(
        "__sub__",
        [](const T& self, const type& other) { return std::make_unique<decltype(self - other)>(self - other); },
        py::is_operator());
    // we can always multiply with a scalar from the right
    // ...
    c.def(
        "__mul__",
        [](const T& self, const Functions::ElementFunctionSetInterface<E, 1, 1, R>& other) {
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
#endif // 0
  } // ... addbind_methods(...)

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id,
                         const std::string& layer_id = "",
                         const std::string& class_id = "element_function_set_interface")
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

    addbind_methods(c);

    return c;
  }
}; // class ElementFunctionSetInterface


} // namespace bindings
} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // PYTHON_DUNE_XT_FUNCTIONS_INTERFACES_ELEMENT_FUNCTION_HH
