// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017, 2020)
//   René Fritze     (2018 - 2020)
//   Tobias Leibner  (2018, 2020 - 2021)

#ifndef DUNE_XT_LA_CONTAINER_VECTOR_INTERFACE_PBH
#define DUNE_XT_LA_CONTAINER_VECTOR_INTERFACE_PBH

#include <sstream>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/operators.h>

#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/common/string.hh>

#include <dune/xt/la/type_traits.hh>
#include <python/dune/xt/la/container.bindings.hh>

#include <dune/xt/la/container/io.hh>
#include <dune/xt/la/container/vector-interface.hh>
#include <python/dune/xt/la/container/container-interface.hh>


namespace Dune::XT::LA {


template <class C>
auto bind_Vector(pybind11::module& m)
{
  static_assert(is_vector<C>::value);
  namespace py = pybind11;
  using namespace pybind11::literals;

  using S = typename C::ScalarType;
  using R = typename C::RealType;

  const auto ClassName = Common::to_camel_case(bindings::container_name<C>::value());

  py::class_<C> c = bind_ProvidesDataAccess<C>(m, ClassName, ClassName);

  c.def(py::init([](const ssize_t size, const S& value) { return new C(Common::numeric_cast<size_t>(size), value); }),
        "size"_a = 0,
        "value"_a = 0.0);
  c.def(py::init([](const py::iterable&& it) {
          std::vector<S> tmp;
          for (py::handle h : it)
            tmp.push_back(h.cast<S>());
          return new C(tmp);
        }),
        "Assigns the elements of the iterable to the vector.");

  c.def(
      "__repr__",
      [ClassName](const C& self) {
        std::stringstream ss;
        ss << ClassName << "([";
        if (self.size() > 0)
          ss << self[0];
        for (size_t ii = 1; ii < std::min(size_t(3), self.size()); ++ii)
          ss << " " << self[ii];
        if (self.size() > max_python_print_rows) {
          ss << " ...";
        } else {
          for (ssize_t ii = std::min(size_t(3), self.size()); ii < ssize_t(self.size()) - 3; ++ii)
            ss << " " << self[ii];
        }
        for (size_t ii = std::max(ssize_t(3), ssize_t(self.size()) - 3); ii < self.size(); ++ii)
          ss << " " << self[ii];
        ss << "])";
        return ss.str();
      },
      "A compact representation of the vector (only the first and last three elements).");
  c.def(
      "__str__",
      [ClassName](const C& self) {
        std::stringstream ss;
        ss << ClassName << "([";
        if (self.size() > 0)
          ss << self[0];
        for (size_t ii = 1; ii < self.size(); ++ii)
          ss << " " << self[ii];
        ss << "])";
        return ss.str();
      },
      "A full representation of the vector.");
  c.def("__len__", [](const C& self) { return self.size(); });
  c.def("__getitem__", [](const C& vec, size_t ii) -> S {
    if (ii >= vec.size())
      throw pybind11::index_error();
    return vec[ii];
  });
  c.def("__setitem__", [](C& vec, size_t ii, const S& value) {
    if (ii >= vec.size())
      throw pybind11::index_error();
    vec[ii] = value;
  });
  c.def(
      "__iter__",
      [](C& vec) {
        return py::
            make_iterator<py::return_value_policy::reference_internal, typename C::iterator, typename C::iterator, S>(
                vec.begin(), vec.end());
      },
      pybind11::keep_alive<0, 1>() /*Essential: keep object alive while iterator exists!*/);

  c.def(
      "__add__", [](const C& self, const C& other) { return std::make_unique<C>(self + other); }, py::is_operator());
  c.def("__iadd__", // function ptr signature required for the right return type
        static_cast<C& (C::*)(const C&)>(&C::operator+=),
        py::is_operator());
  c.def(
      "__sub__", [](const C& self, const C& other) { return std::make_unique<C>(self - other); }, py::is_operator());
  c.def("__isub__", // function ptr signature required for the right return type
        static_cast<C& (C::*)(const C&)>(&C::operator-=),
        py::is_operator());
  c.def(
      "__mul__",
      [](const C& self, const R& alpha) {
        auto ret = std::make_unique<C>(self.copy());
        (*ret) *= alpha;
        return ret;
      },
      py::is_operator());
  c.def(
      "__rmul__",
      [](const C& self, const R& alpha) {
        auto ret = std::make_unique<C>(self.copy());
        (*ret) *= alpha;
        return ret;
      },
      py::is_operator());
  // c.def("__imul__", // function ptr signature required for the right return type
  //       py::overload_cast<const C&>(&C::operator*=),
  //       py::is_operator());
  c.def(
      "__truediv__",
      [](const C& self, const R& alpha) {
        auto ret = std::make_unique<C>(self.copy());
        (*ret) /= alpha;
        return ret;
      },
      py::is_operator());
  // c.def("__itruediv__", // function ptr signature required for the right return type
  //       py::overload_cast<const C&>(&C::operator/=),
  //       py::is_operator());
  c.def(
      "neg",
      [](const C& self) {
        auto ret = std::make_unique<C>(self.copy());
        (*ret) *= -1;
        return ret;
      },
      py::is_operator());

  c.def_property_readonly("size", [](const C& self) { return self.size(); });
  c.def(
      "add_to_entry",
      [](C& self, const ssize_t ii, const S& value) { self.add_to_entry(Common::numeric_cast<size_t>(ii), value); },
      "ii"_a,
      "value"_a);
  c.def(
      "set_entry",
      [](C& self, const ssize_t ii, const S& value) { self.set_entry(Common::numeric_cast<size_t>(ii), value); },
      "ii"_a,
      "value"_a);
  c.def(
      "get_entry",
      [](const C& self, const ssize_t ii) { return self.get_entry(Common::numeric_cast<size_t>(ii)); },
      "jj"_a);
  c.def(
      "set_all", [](C& self, const S& value) { self.set_all(value); }, "value"_a);
  c.def("valid", [](const C& self) { return self.valid(); });
  c.def("dim", [](const C& self) { return self.size(); });
  c.def("mean", [](const C& self) { return self.mean(); });
  c.def("amax", [](const C& self) { return self.amax(); });
  c.def(
      "almost_equal",
      [](const C& self, const C& other, const S& epsilon) { return self.almost_equal(other, epsilon); },
      "other"_a,
      "epsilon"_a = Common::FloatCmp::DefaultEpsilon<S>::value());
  c.def("dot", [](const C& self, const C& other) { return self.dot(other); });
  c.def("l1_norm", [](const C& self) { return self.l1_norm(); });
  c.def("l2_norm", [](const C& self) { return self.l2_norm(); });
  c.def("sup_norm", [](const C& self) { return self.sup_norm(); });
  c.def("standard_deviation", [](const C& self) { return self.standard_deviation(); });
  c.def(
      "to_file",
      [](const C& self, const std::string& filename, const std::string& mode) { to_file(self, filename, mode); },
      "filename"_a,
      "mode"_a = "ascii");
  c.def_static(
      "from_file",
      [](const std::string& filename, const ssize_t min_size, const std::string& mode) {
        return from_file<C>(filename, min_size, mode);
      },
      "filename"_a,
      "min_size"_a = -1,
      "mode"_a = "ascii");

  c.def(py::pickle([](const C& self) { return py::make_tuple(std::vector<S>(self)); },
                   [](const py::tuple&& t) {
                     if (t.size() != 1)
                       throw std::runtime_error("Invalid state!");
                     const auto data = t[0].cast<std::vector<S>>();
                     C* ret = new C(data.size());
                     /* Assign any additional state */
                     for (size_t ii = 0; ii < ret->size(); ++ii)
                       (*ret)[ii] = data[ii];
                     return ret;
                   }));

  addbind_ContainerInterface(c);

  return c;
} // ... bind_Vector(...)


} // namespace Dune::XT::LA

#endif // DUNE_XT_LA_CONTAINER_VECTOR_INTERFACE_PBH
