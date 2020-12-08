// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/expression/parametric.hh>

#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/fmatrix.hh>
#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/grid/traits.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>

namespace Dune::XT::Functions::bindings {


template <size_t d, size_t r = 1, class R = double>
class ParametricExpressionFunction
{
  static constexpr size_t rC = 1;

  using type = Functions::ParametricExpressionFunction<d, r, rC, R>;
  using base_type = Functions::FunctionInterface<d, r, rC, R>;
  using bound_type = pybind11::class_<type, base_type>;

public:
  static bound_type bind(pybind11::module& m, const std::string& class_id = "parametric_expression_function")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id + "_" + Common::to_string(d) + "d_to_" + Common::to_string(r);
    if (rC > 1)
      class_name += "x" + Common::to_string(rC);
    class_name += "d";
    const auto ClassName = Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), XT::Common::to_camel_case(class_id).c_str());
    c.def(py::init<const std::string&,
                   const Common::ParameterType&,
                   const Common::FieldVector<std::string, r>&,
                   const size_t,
                   const std::string>(),
          "variable"_a,
          "param_type"_a,
          "expressions"_a,
          "order"_a,
          "name"_a = type::static_id());

    m.def(
        Common::to_camel_case(class_id).c_str(),
        [](Grid::bindings::Dimension<d> /*dim_domain*/,
           const std::string& variable,
           const Common::ParameterType& param_type,
           const Common::FieldVector<std::string, r>& expressions,
           const size_t order,
           const std::string& name) { return new type(variable, param_type, expressions, order, name); },
        "dim_domain"_a,
        "variable"_a,
        "param_type"_a,
        "expressions"_a,
        "order"_a,
        "name"_a = type::static_id());

    return c;
  }
}; // class ParametricExpressionFunction


} // namespace Dune::XT::Functions::bindings


PYBIND11_MODULE(_functions_parametric_expression, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.functions._functions_function_interface_1d");
  py::module::import("dune.xt.functions._functions_function_interface_2d");
  py::module::import("dune.xt.functions._functions_function_interface_3d");

  Dune::XT::Functions::bindings::ParametricExpressionFunction<1, 1>::bind(m);
  Dune::XT::Functions::bindings::ParametricExpressionFunction<1, 2>::bind(m);
  Dune::XT::Functions::bindings::ParametricExpressionFunction<1, 3>::bind(m);

  Dune::XT::Functions::bindings::ParametricExpressionFunction<2, 1>::bind(m);
  Dune::XT::Functions::bindings::ParametricExpressionFunction<2, 2>::bind(m);
  Dune::XT::Functions::bindings::ParametricExpressionFunction<2, 3>::bind(m);

  Dune::XT::Functions::bindings::ParametricExpressionFunction<3, 1>::bind(m);
  Dune::XT::Functions::bindings::ParametricExpressionFunction<3, 2>::bind(m);
  Dune::XT::Functions::bindings::ParametricExpressionFunction<3, 3>::bind(m);
} // PYBIND11_MODULE(...)
