// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017, 2019)
//   René Fritze     (2018 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2018, 2020)

#ifndef PYTHON_DUNE_XT_FUNCTIONS_SPE10_HH
#define PYTHON_DUNE_XT_FUNCTIONS_SPE10_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/gridprovider/provider.hh>

#include <dune/xt/functions/spe10/model1.hh>

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
template <class G, size_t d, size_t r, size_t rC>
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
#if HAVE_DUNE_XT_DATA
  c.def(py::init([](const Common::FieldVector<D, d>& ll,
                    const Common::FieldVector<D, d>& up,
                    const R& min,
                    const R& max,
                    const std::string nm) {
          return std::make_unique<C>(XT::Data::spe10_model1_filename(), ll, up, min, max, nm);
        }),
        "lower_left"_a,
        "upper_right"_a,
        "min"_a = Spe10::internal::model1_min_value,
        "max"_a = Spe10::internal::model1_max_value,
        "name"_a = C::static_id());
#endif // HAVE_DUNE_XT_DATA

  c.def_property_readonly("static_id", [](const C& /*self*/) { return C::static_id(); });

  const std::string make_name = "make_spe10_model1_function_" + Common::to_string(r) + "x" + Common::to_string(rC);
  m.def(
      std::string(make_name).c_str(),
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

  return c;
} // ... bind_Spe10Model1Function(...)

template <class G, size_t d, size_t r, size_t rC>
void bind_Spe10Model1Function(pybind11::module& m, const std::string& grid_id)
{
  if constexpr (Grid::is_grid<G>::value && d == 2) {
    bind_Spe10Model1Function_2D<G, d, r, rC>(m, grid_id);
  }
}

} // namespace Dune::XT::Functions

#endif // PYTHON_DUNE_XT_FUNCTIONS_SPE10_HH
