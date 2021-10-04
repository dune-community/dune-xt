// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017, 2019 - 2020)
//   René Fritze     (2018 - 2020)
//   Tim Keil        (2018)
//   Tobias Leibner  (2018, 2020)

#ifndef PYTHON_DUNE_XT_FUNCTIONS_SPE10_HH
#define PYTHON_DUNE_XT_FUNCTIONS_SPE10_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/gridprovider/provider.hh>

#include <dune/xt/functions/spe10/model1.hh>
#include <dune/xt/functions/spe10/model2.hh>

#include <python/dune/xt/common/fvector.hh>

namespace Dune::XT::Functions {


/**
 * \note We would like to drop the d template parameter and use either of
\code
static constexpr           size_t d = G::dimension;
static constexpr size_t d = G::dimension;
\endcode
 *       but this triggers a bug in gcc-4.9 and we thus need to use G::dimension
 *       everywhere: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=59937
 */
template <class G, size_t d, size_t r = 1, size_t rC = 1>
auto bind_Spe10Model1Function_2D(pybind11::module& m, const std::string& grid_id)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  using E = typename G::template Codim<0>::Entity;
  using D = typename G::ctype;
  using R = double;

  using I = GridFunctionInterface<E, r, rC, R>;
  using C = Spe10::Model1Function<E, r, rC, R>;

  py::class_<C, I> c(
      m,
      std::string("Spe10Model1Function__" + grid_id + "_to_" + Common::to_string(r) + "x" + Common::to_string(rC))
          .c_str(),
      std::string("Spe10Model1Function__" + grid_id + "_to_" + Common::to_string(r) + "x" + Common::to_string(rC))
          .c_str());

  c.def_property_readonly_static("num_elements", [=](py::object /*self*/) {
    return std::array<size_t, 2>{Spe10::internal::model1_x_elements /*, Spe10::internal::model1_y_elements*/,
                                 Spe10::internal::model1_z_elements};
  });
  c.def_property_readonly_static("domain", [=](py::object /*self*/) {
    return std::make_tuple(std::array<double, 2>{0., 0.},
                           std::array<double, 2>{Spe10::internal::model_2_length_x,
                                                 /*Spe10::internal::model_2_length_y,*/
                                                 Spe10::internal::model_2_length_z});
  });

  c.def(py::init<const std::string&,
                 const Common::FieldVector<D, d>&,
                 const Common::FieldVector<D, d>&,
                 const R&,
                 const R&,
                 const std::string>(),
        "filename"_a,
        "lower_left"_a,
        "upper_right"_a,
        "min"_a = Spe10::internal::model1_min_value,
        "max"_a = Spe10::internal::model1_max_value,
        "name"_a = C::static_id());
#if HAVE_SPE10_DATA
  c.def(py::init([](const Common::FieldVector<D, d>& ll,
                    const Common::FieldVector<D, d>& up,
                    const R& min,
                    const R& max,
                    const std::string nm) { return std::make_unique<C>(SPE10_MODEL1_FILENAME, ll, up, min, max, nm); }),
        "lower_left"_a,
        "upper_right"_a,
        "min"_a = Spe10::internal::model1_min_value,
        "max"_a = Spe10::internal::model1_max_value,
        "name"_a = C::static_id());
#endif // HAVE_SPE10_DATA

  c.def_property_readonly("static_id", [](const C& /*self*/) { return C::static_id(); });

  m.def(
      "Spe10Model1Function",
      [](const Grid::GridProvider<G>& /*grid*/,
         const std::string& filename,
         const Common::FieldVector<D, d>& lower_left,
         const Common::FieldVector<D, d>& upper_right,
         const R& min,
         const R& max,
         const std::string& name) { return C(filename, lower_left, upper_right, min, max, name); },
      "grid_provider"_a,
      "filename"_a,
      "lower_left"_a,
      "upper_right"_a,
      "min"_a = Spe10::internal::model1_min_value,
      "max"_a = Spe10::internal::model1_max_value,
      "name"_a = C::static_id());
#if HAVE_SPE10_DATA
  m.def(
      "Spe10Model1Function",
      [](const Grid::GridProvider<G>& /*grid*/,
         const std::string& filename,
         const Common::FieldVector<D, d>& lower_left,
         const Common::FieldVector<D, d>& upper_right,
         const R& min,
         const R& max,
         const std::string& name) { return C(filename, lower_left, upper_right, min, max, name); },
      "grid_provider"_a,
      "filename"_a,
      "lower_left"_a,
      "upper_right"_a,
      "min"_a = Spe10::internal::model1_min_value,
      "max"_a = Spe10::internal::model1_max_value,
      "name"_a = C::static_id());
#endif // HAVE_SPE10_DATA

  return c;
} // ... bind_Spe10Model1Function(...)

template <class G, size_t d, size_t r = 1, size_t rC = 1>
void bind_Spe10Model1Function(pybind11::module& m, const std::string& grid_id)
{
  if constexpr (Grid::is_grid<G>::value && d == 2) {
    bind_Spe10Model1Function_2D<G, d, r, rC>(m, grid_id);
  }
}


template <class G>
auto bind_Spe10Model2Function(pybind11::module& m, const std::string& grid_id)
{
  static_assert(G::dimension == 3, "");
  namespace py = pybind11;
  using namespace pybind11::literals;

  using E = typename G::template Codim<0>::Entity;
  using D = typename G::ctype;
  static constexpr size_t d = 3;
  static constexpr size_t r = 3;
  static constexpr size_t rC = 3;
  using R = double;

  using I = GridFunctionInterface<E, 3, 3, R>;
  using C = Spe10::Model2Function<E, 3, 3, R>;

  py::class_<C, I> c(
      m,
      std::string("Spe10Model2Function__" + grid_id + "_to_" + Common::to_string(r) + "x" + Common::to_string(rC))
          .c_str(),
      std::string("Spe10Model2Function__" + grid_id + "_to_" + Common::to_string(r) + "x" + Common::to_string(rC))
          .c_str());

  c.def_property_readonly_static("num_elements", [=](py::object /*self*/) {
    return std::array<size_t, 3>{
        Spe10::internal::model2_x_elements, Spe10::internal::model2_y_elements, Spe10::internal::model2_z_elements};
  });
  c.def_property_readonly_static("domain", [=](py::object /*self*/) {
    return std::make_tuple(std::array<double, 3>{0., 0., 0.},
                           std::array<double, 3>{Spe10::internal::model_2_length_x,
                                                 Spe10::internal::model_2_length_y,
                                                 Spe10::internal::model_2_length_z});
  });

  c.def(py::init<const std::string&, const Common::FieldVector<D, d>&, const Common::FieldVector<D, d>&>(),
        "filename"_a,
        "lower_left"_a,
        "upper_right"_a);
#if HAVE_SPE10_DATA
  c.def(py::init([](const Common::FieldVector<D, d>& ll, const Common::FieldVector<D, d>& up) {
          return std::make_unique<C>(SPE10_MODEL2_FILENAME, ll, up);
        }),
        "lower_left"_a,
        "upper_right"_a);
#endif // HAVE_SPE10_DATA

  c.def_property_readonly("static_id", [](const C& /*self*/) { return C::static_id(); });

  m.def("Spe10Model2FunctionType", [=](const Grid::GridProvider<G>& /*grid*/) { return c; });

  m.def(
      "Spe10Model2Function",
      [](const Grid::GridProvider<G>& /*grid*/,
         const std::string& filename,
         const Common::FieldVector<D, d>& lower_left,
         const Common::FieldVector<D, d>& upper_right,
         const std::string& name) {
        return C(filename,
                 lower_left,
                 upper_right,
                 {Spe10::internal::model2_x_elements,
                  Spe10::internal::model2_y_elements,
                  Spe10::internal::model2_z_elements},
                 name);
      },
      "grid_provider"_a,
      "filename"_a,
      "lower_left"_a,
      "upper_right"_a,
      "name"_a = C::static_id());
#if HAVE_SPE10_DATA
  m.def(
      "Spe10Model2Function",
      [](const Grid::GridProvider<G>& /*grid*/,
         const Common::FieldVector<D, d>& lower_left,
         const Common::FieldVector<D, d>& upper_right,
         const std::string& name) {
        return C(SPE10_MODEL2_FILENAME,
                 lower_left,
                 upper_right,
                 {Spe10::internal::model2_x_elements,
                  Spe10::internal::model2_y_elements,
                  Spe10::internal::model2_z_elements},
                 name);
      },
      "grid_provider"_a,
      "lower_left"_a,
      "upper_right"_a,
      "name"_a = C::static_id());
#endif // HAVE_SPE10_DATA

  return c;


} // namespace Dune::XT::Functions

#endif // PYTHON_DUNE_XT_FUNCTIONS_SPE10_HH
