// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017, 2020)
//   René Fritze     (2018 - 2020)
//   Tobias Leibner  (2020 - 2021)

#ifndef DUNE_XT_LA_CONTAINER_MATRIX_INTERFACE_PBH
#define DUNE_XT_LA_CONTAINER_MATRIX_INTERFACE_PBH

#include <sstream>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/operators.h>

#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/common/string.hh>

#include <dune/xt/la/container/io.hh>
#include <dune/xt/la/type_traits.hh>

#include <dune/xt/la/container/matrix-interface.hh>

#include <python/dune/xt/la/container/container-interface.hh>
#include <python/dune/xt/la/container.bindings.hh>

namespace Dune::XT::LA {
namespace internal {


template <class T, bool sparse = false>
struct addbind_Matrix
{
  template <class S, class C>
  static void ctor(C& c)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    c.def(py::init([](const ssize_t rows, const ssize_t cols, const S& value) {
            return new T(Common::numeric_cast<size_t>(rows), Common::numeric_cast<size_t>(cols), value);
          }),
          "rows"_a = 0,
          "cols"_a = 0,
          "value"_a = 0.0);
  }
};

template <class T>
struct addbind_Matrix<T, true>
{
  template <class S, class C>
  static void ctor(C& /*c*/)
  {}
};


template <class M>
void print_row_sparsely(const M& self, const size_t row, std::stringstream& ss)
{
  ss << "[";
  if (self.cols() > 0)
    ss << self.get_entry(row, 0);
  for (size_t ii = 1; ii < std::min(size_t(3), self.cols()); ++ii)
    ss << " " << self.get_entry(row, ii);
  if (self.cols() > max_python_print_rows) {
    ss << " ...";
  } else {
    for (ssize_t ii = std::min(size_t(3), self.cols()); ii < ssize_t(self.cols()) - 3; ++ii)
      ss << " " << self.get_entry(row, ii);
  }
  for (size_t ii = std::max(ssize_t(3), ssize_t(self.cols()) - 3); ii < self.cols(); ++ii)
    ss << " " << self.get_entry(row, ii);
  ss << "]";
} // ... print_row_sparsely(...)


} // namespace internal


template <class C, bool sparse>
auto bind_Matrix(pybind11::module& m)
{
  static_assert(is_matrix<C>::value);
  namespace py = pybind11;
  using namespace pybind11::literals;

  using S = typename C::ScalarType;

  const auto ClassName = Common::to_camel_case(bindings::container_name<C>::value());

  py::class_<C> c = bind_ProvidesDataAccess<C>(m, ClassName, ClassName);

  internal::addbind_Matrix<C, sparse>::template ctor<S>(c);
  c.def(py::init([](const ssize_t rows, const ssize_t cols, const SparsityPatternDefault& pattern) {
          return new C(Common::numeric_cast<size_t>(rows), Common::numeric_cast<size_t>(cols), pattern);
        }),
        "rows"_a,
        "cols"_a,
        "pattern"_a);

  c.def(
      "__repr__",
      [ClassName](const C& self) {
        std::stringstream ss;
        ss << ClassName << "(";
        if (self.rows() > 0) {
          ss << "\n  ";
          internal::print_row_sparsely(self, 0, ss);
          ss << "\n";
        }
        for (size_t ii = 1; ii < std::min(size_t(3), self.rows()); ++ii) {
          ss << "  ";
          internal::print_row_sparsely(self, ii, ss);
          ss << "\n";
        }
        if (self.rows() > max_python_print_rows) {
          ss << "   ...\n";
        } else {
          for (ssize_t ii = std::min(size_t(3), self.rows()); ii < ssize_t(self.rows()) - 3; ++ii) {
            ss << "  ";
            internal::print_row_sparsely(self, ii, ss);
            ss << "\n";
          }
        }
        for (size_t ii = std::max(ssize_t(3), ssize_t(self.rows()) - 3); ii < self.rows(); ++ii) {
          ss << "  ";
          internal::print_row_sparsely(self, ii, ss);
          ss << "\n";
        }
        ss << ")";
        return ss.str();
      },
      "A compact representation of the matrix (only the outer three elements).");
  c.def(
      "__str__",
      [ClassName](const C& self) {
        std::stringstream ss;
        ss << ClassName << "(\n" << self << "\n)";
        return ss.str();
      },
      "A full representation of the matrix.");

  c.def_property_readonly("rows", [](const C& self) { return self.rows(); });
  c.def_property_readonly("cols", [](const C& self) { return self.cols(); });
  c.def(
      "add_to_entry",
      [](C& self, const ssize_t ii, const ssize_t jj, const S& value) {
        self.add_to_entry(Common::numeric_cast<size_t>(ii), Common::numeric_cast<size_t>(jj), value);
      },
      "ii"_a,
      "jj"_a,
      "value"_a);
  c.def(
      "set_entry",
      [](C& self, const ssize_t ii, const ssize_t jj, const S& value) {
        self.set_entry(Common::numeric_cast<size_t>(ii), Common::numeric_cast<size_t>(jj), value);
      },
      "ii"_a,
      "jj"_a,
      "value"_a);
  c.def(
      "get_entry",
      [](const C& self, const ssize_t ii, const ssize_t jj) {
        return self.get_entry(Common::numeric_cast<size_t>(ii), Common::numeric_cast<size_t>(jj));
      },
      "ii"_a,
      "jj"_a);
  c.def(
      "clear_row", [](C& self, const ssize_t ii) { self.clear_row(Common::numeric_cast<size_t>(ii)); }, "ii"_a);
  c.def(
      "clear_col", [](C& self, const ssize_t jj) { self.clear_row(Common::numeric_cast<size_t>(jj)); }, "jj"_a);
  c.def(
      "unit_row", [](C& self, const ssize_t ii) { self.unit_row(Common::numeric_cast<size_t>(ii)); }, "ii"_a);
  c.def(
      "unit_col", [](C& self, const ssize_t jj) { self.unit_col(Common::numeric_cast<size_t>(jj)); }, "jj"_a);
  c.def("valid", [](const C& self) { return self.valid(); });
  c.def("sup_norm", [](const C& self) { return self.sup_norm(); });
  c.def_property_readonly("non_zeros", [](const C& self) { return self.non_zeros(); });
  c.def(
      "pattern",
      [](const C& self, const bool prune, const S& eps) { return self.pattern(prune, eps); },
      "prune"_a = false,
      "eps"_a = Common::FloatCmp::DefaultEpsilon<S>::value());
  c.def(
      "pruned",
      [](const C& self, const S& eps) { return self.pruned(eps); },
      "eps"_a = Common::FloatCmp::DefaultEpsilon<S>::value());
  c.def("almost_equal", [](const C& self, const C& other) { return self.almost_equal(other); });
  c.def(
      "to_file",
      [](const C& self, const std::string& filename, const std::string& mode) { to_file(self, filename, mode); },
      "filename"_a,
      "mode"_a = "ascii");
  c.def_static(
      "from_file",
      [](const std::string& filename, const ssize_t min_rows, const ssize_t min_cols, const std::string& mode) {
        return from_file<C>(filename, min_rows, min_cols, mode);
      },
      "filename"_a,
      "min_rows"_a = -1,
      "min_cols"_a = -1,
      "mode"_a = "ascii");
  c.def(
      "to_csr",
      [](const C& self, const bool prune, const S& eps) {
        const auto pattern = self.pattern(prune, eps);
        size_t nnz = 0;
        for (const auto& row : pattern)
          nnz += row.size();
        std::vector<S> data;
        std::vector<size_t> row_indices;
        std::vector<size_t> col_indices;
        data.reserve(nnz);
        row_indices.reserve(nnz);
        col_indices.reserve(nnz);
        for (size_t row_index = 0; row_index < self.rows(); ++row_index) {
          for (const auto& col_index : pattern.inner(row_index)) {
            data.push_back(self.get_entry(row_index, col_index));
            row_indices.push_back(row_index);
            col_indices.push_back(col_index);
          }
        }
        return std::make_tuple(std::move(data), std::move(row_indices), std::move(col_indices));
      },
      py::call_guard<py::gil_scoped_release>(),
      "prune"_a = false,
      "eps"_a = 1e-15);
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
      [](const C& self, const S& alpha) {
        auto ret = std::make_unique<C>(self.copy());
        (*ret) *= alpha;
        return ret;
      },
      py::is_operator());
  c.def(
      "__rmul__",
      [](const C& self, const S& alpha) {
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
      [](const C& self, const S& alpha) {
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

  addbind_ContainerInterface(c);

  return c;
} // ... bind_Matrix(...)


template <class M, class V>
void addbind_Matrix_Vector_interaction(pybind11::class_<M>& mat, pybind11::class_<V>& vec)
{
  mat.def(
      "mv", [](const M& self, const V& xx, V& yy) { self.mv(xx, yy); }, "source", "range");
  mat.def(
      "dot",
      [](const M& self, const V& xx) {
        V yy(self.rows());
        self.mv(xx, yy);
        return yy;
      },
      "source");
  mat.def(
      "__matmul__",
      [](const M& self, const V& xx) {
        V yy(self.rows());
        self.mv(xx, yy);
        return yy;
      },
      "source");
  mat.def(
      "mtv", [](const M& self, const V& xx, V& yy) { self.mtv(xx, yy); }, "source", "range");

  mat.def("vector_type", [vec](M& /*self*/) { return vec; });
  vec.def("matrix_type", [mat](V& /*self*/) { return mat; });
} // ... addbind_Matrix_Vector_interaction(...)


} // namespace Dune::XT::LA

#endif // DUNE_XT_LA_CONTAINER_MATRIX_INTERFACE_PBH
