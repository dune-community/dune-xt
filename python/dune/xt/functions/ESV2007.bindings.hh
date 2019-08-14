// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018)
//   Tim Keil        (2018)

#ifndef DUNE_XT_FUNCTIONS_ESV2007_BINDINGS_HH
#define DUNE_XT_FUNCTIONS_ESV2007_BINDINGS_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/xt/functions/ESV2007.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace bindings {
namespace ESV2007 {


template <class G>
class CutoffFunction
{
  static_assert(XT::Grid::is_grid<G>::value, "");
  typedef typename G::ctype D;
  static const size_t d = G::dimension;

  template <bool is_correct_dim = (d == 2), bool anything = false>
  struct helper
  {
    static void bind(pybind11::module& m)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      typedef typename G::template Codim<0>::Entity E;
      typedef double R;
      typedef GridFunctionInterface<E, 1, 1, R> ScalarFunction;
      typedef GridFunctionInterface<E, d, d, R> TensorFunction;
      typedef Functions::ESV2007::CutoffFunction<ScalarFunction, void> type_single_diffusion;
      typedef Functions::ESV2007::CutoffFunction<ScalarFunction, TensorFunction> type;

      py::class_<type_single_diffusion, ScalarFunction> c_single_diffusion(
          m,
          Common::to_camel_case("esv2007_cutoff_function_single_diffusion_from_"
                                + XT::Grid::bindings::grid_name<G>::value() + "_to_1x1")
              .c_str(),
          "ESV2007::CutoffFunction");

      c_single_diffusion.def(py::init<const ScalarFunction&, const R, const std::string>(),
                             "diffusion"_a,
                             "poincare_constant"_a = 1.0 / (M_PIl * M_PIl),
                             "nm"_a = type_single_diffusion::static_id());

      c_single_diffusion.def_property_readonly(
          "static_id", [](const type_single_diffusion& /*self*/) { return type_single_diffusion::static_id(); });

      py::class_<type, ScalarFunction> c(
          m,
          Common::to_camel_case("esv2007_cutoff_function_diffusion_factor_and_tensor_from_"
                                + XT::Grid::bindings::grid_name<G>::value() + "_to_1x1")
              .c_str(),
          "ESV2007::CutoffFunction");

      c.def(py::init<const ScalarFunction&, const TensorFunction&, const R, const std::string>(),
            "diffusion_factor"_a,
            "diffusion_tensor"_a,
            "poincare_constant"_a = 1.0 / (M_PIl * M_PIl),
            "nm"_a = type::static_id());
      c.def_property_readonly("static_id", [](const type& /*self*/) { return type::static_id(); });

      const std::string make_name = "make_esv2007_cutoff_function";
      m.def(std::string(make_name + "_single_diffusion_to_1x1").c_str(),
            [](const Grid::GridProvider<G>& /*grid*/,
               const ScalarFunction& diffusion,
               const R& poincare_constant,
               const std::string& name) { return type_single_diffusion(diffusion, poincare_constant, name); },
            "grid_provider"_a,
            "diffusion"_a,
            "poincare_constant"_a = 1.0 / (M_PIl * M_PIl),
            "name"_a = type_single_diffusion::static_id(),
            py::keep_alive<0, 2>());
      m.def(std::string(make_name + "_diffusion_factor_and_tensor_to_1x1").c_str(),
            [](const Grid::GridProvider<G>& /*grid*/,
               const ScalarFunction& diffusion_factor,
               const TensorFunction& diffusion_tensor,
               const R& poincare_constant,
               const std::string& name) { return type(diffusion_factor, diffusion_tensor, poincare_constant, name); },
            "grid_provider"_a,
            "diffusion_factor"_a,
            "diffusion_tesor"_a,
            "poincare_constant"_a = 1.0 / (M_PIl * M_PIl),
            "name"_a = type::static_id(),
            py::keep_alive<0, 2>(),
            py::keep_alive<0, 3>());
    }
  }; // struct helper<true, ...>

  template <bool anything>
  struct helper<false, anything>
  {
    static void bind(pybind11::module& /*m*/) {}
  };

public:
  static void bind(pybind11::module& m)
  {
    helper<>::bind(m);
  }
}; // class CutoffFunction


} // namespace ESV2007
} // namespace bindings
} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_ESV2007_BINDINGS_HH
