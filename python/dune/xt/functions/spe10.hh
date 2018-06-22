// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_XT_FUNCTIONS_SPE10_PBH
#define DUNE_XT_FUNCTIONS_SPE10_PBH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>

#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/gridprovider/provider.hh>

#include <dune/xt/functions/spe10/model1.hh>

#include <python/dune/xt/common/fvector.hh>

namespace Dune {
namespace XT {
namespace Functions {


template <class G, size_t d, size_t r, size_t rC>
typename std::enable_if<Grid::is_grid<G>::value && d != 2, void>::type
bind_Spe10Model1Function(pybind11::module& /*m*/, const std::string& /*grid_id*/)
{
}

/**
 * \note We would like to drop the d template parameter and use either of
\code
static const           size_t d = G::dimension;
static const constexpr size_t d = G::dimension;
\endcode
 *       but this triggers a bug in gcc-4.9 and we thus need to use G::dimension
 *       everywhere: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=59937
 */
template <class G, size_t d, size_t r, size_t rC>
typename std::enable_if<Grid::is_grid<G>::value && d == 2,
                        pybind11::class_<Spe10::Model1Function<typename G::template Codim<0>::Entity, r, rC, double>>>::
    type
    bind_Spe10Model1Function(pybind11::module& m, const std::string& grid_id)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  typedef double R;

  typedef GridFunctionInterface<E, r, rC, R> I;
  typedef Spe10::Model1Function<E, r, rC, R> C;

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

  c.def_property_readonly("static_id", [](const C& /*self*/) { return C::static_id(); });

  const std::string make_name = "make_spe10_model1_function_" + Common::to_string(r) + "x" + Common::to_string(rC);
  m.def(std::string(make_name).c_str(),
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
  m.def(std::string(make_name).c_str(),
        [](const Grid::GridProvider<G, Grid::DD::SubdomainGrid<G>>& /*grid*/,
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


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_SPE10_PBH
