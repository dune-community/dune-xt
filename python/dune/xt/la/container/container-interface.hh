// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017, 2020)
//   René Fritze     (2018 - 2020)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_LA_CONTAINER_INTERFACE_PBH
#define DUNE_XT_LA_CONTAINER_INTERFACE_PBH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/operators.h>

#include <dune/xt/la/type_traits.hh>

#include <dune/xt/la/container/container-interface.hh>

namespace Dune {
namespace XT {
namespace LA {

constexpr size_t max_python_print_rows{8};

pybind11::enum_<Backends> bind_Backends(pybind11::module& m)
{
  namespace py = pybind11;

  py::enum_<Backends> c(m, "Backends");
  c.value("common_dense", Backends::common_dense);
  c.value("common_sparse", Backends::common_sparse);
  c.value("istl_sparse", Backends::istl_sparse);
  c.value("eigen_dense", Backends::eigen_dense);
  c.value("eigen_sparse", Backends::eigen_sparse);
  c.value("none", Backends::none);

  m.attr("default_backend") = py::cast(default_backend);
  m.attr("default_sparse_backend") = py::cast(default_sparse_backend);
  m.attr("default_dense_backend") = py::cast(default_dense_backend);

  return c;
} // ... bind_Backends(...)


template <class C>
void addbind_ContainerInterface(pybind11::class_<C>& c)
{
  static_assert(is_container<C>::value);
  namespace py = pybind11;
  using namespace pybind11::literals;

  using S = typename C::ScalarType;
#ifndef __clang_analyzer__ // tidy throws a false positive on the lambda
  c.def(
      "copy",
      [](C& self, const bool deep) {
        if (deep)
          return self.copy();
        return C(self);
      },
      "deep"_a = false);
#endif
  c.def(
      "scal", [](C& self, const S& alpha) { self.scal(alpha); }, "alpha"_a);
  c.def("axpy", [](C& self, const S& alpha, const C& xx) { self.axpy(alpha, xx); });
  c.def(
      "has_equal_shape", [](const C& self, const C& other) { return self.has_equal_shape(other); }, "other"_a);
  c.def(py::self *= S());
  c.def(py::self * S());
  c.def("__sub__", [](const C& self, const C& other) {
    auto ret = self.copy();
    ret.axpy(-1, other);
    return ret;
  });

} // ... addbind_ContainerInterface(...)

template <class C>
auto bind_ProvidesDataAccess(pybind11::module& m, const std::string& class_id, const std::string& help_id)
{
  namespace py = pybind11;
  if constexpr (!provides_data_access<C>::value) {
    return py::class_<C>(m, class_id.c_str(), help_id.c_str());
  } else {
    namespace py = pybind11;
    using D = typename C::DataType;
    py::class_<C> c(m, class_id.c_str(), help_id.c_str(), py::buffer_protocol());
    if constexpr (is_vector<C>::value) {
      c.def_buffer([](C& vec) -> py::buffer_info {
        return py::buffer_info(
            vec.data(), sizeof(D), py::format_descriptor<D>::format(), 1, {vec.data_size()}, {sizeof(D)});
      });
    } else if constexpr (is_matrix<C>::value) {
      c.def_buffer([](C& mat) -> py::buffer_info {
        return py::buffer_info(mat.data(), /* Pointer to buffer */
                               sizeof(D), /* Size of one scalar */
                               py::format_descriptor<D>::format(), /* Python struct-style format descriptor */
                               2, /* Number of dimensions */
                               {mat.rows(), mat.cols()}, /* Buffer dimensions */
                               {sizeof(D) * mat.cols(), /* Strides (in bytes) for each index */
                                sizeof(D)});
      });
    }
    return c;
  }
}

} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_INTERFACE_PBH
