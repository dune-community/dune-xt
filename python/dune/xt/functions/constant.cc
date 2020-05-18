// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#include "config.h"

#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/constant.hh>

#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/fmatrix.hh>
#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/grid/traits.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace bindings {


template <size_t d, size_t r = 1, size_t rC = 1, class R = double>
class ConstantFunction
{
  using type = Functions::ConstantFunction<d, r, rC, R>;
  using base_type = Functions::FunctionInterface<d, r, rC, R>;
  using bound_type = pybind11::class_<type, base_type>;

public:
  static bound_type bind(pybind11::module& m, const std::string& class_id = "constant_function")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id + "_" + Common::to_string(d) + "_to_" + Common::to_string(r);
    if (rC > 1)
      class_name += "x" + Common::to_string(rC);
    class_name += "d";
    const auto ClassName = Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), XT::Common::to_camel_case(class_id).c_str());
    c.def(
        py::init<const typename type::RangeReturnType&, const std::string&>(), "value"_a, "name"_a = type::static_id());

    if (rC == 1)
      m.def(XT::Common::to_camel_case(class_id).c_str(),
            [](Grid::bindings::Dimension<d> /*dim_domain*/,
               Grid::bindings::Dimension<r> /*dim_range*/,
               const typename type::RangeReturnType& value,
               const std::string& name) { return type(value, name); },
            "dim_domain"_a,
            "dim_range"_a,
            "value"_a,
            "name"_a = type::static_id());
    else
      m.def(XT::Common::to_camel_case(class_id).c_str(),
            [](Grid::bindings::Dimension<d> /*dim_domain*/,
               std::pair<Grid::bindings::Dimension<r>, Grid::bindings::Dimension<rC>> /*dim_range*/,
               const typename type::RangeReturnType& value,
               const std::string& name) { return type(value, name); },
            "dim_domain"_a,
            "dim_range"_a,
            "value"_a,
            "name"_a = type::static_id());

    return c;
  }
}; // class ConstantFunction


} // namespace bindings
} // namespace Functions
} // namespace XT
} // namespace Dune


PYBIND11_MODULE(_functions_constant, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.functions._functions_function_interface_1d");
  py::module::import("dune.xt.functions._functions_function_interface_2d");
  py::module::import("dune.xt.functions._functions_function_interface_3d");

  Dune::XT::Functions::bindings::ConstantFunction<1, 1, 1>::bind(m);
  Dune::XT::Functions::bindings::ConstantFunction<1, 2, 1>::bind(m);
  Dune::XT::Functions::bindings::ConstantFunction<1, 2, 2>::bind(m);
  Dune::XT::Functions::bindings::ConstantFunction<1, 3, 1>::bind(m);
  Dune::XT::Functions::bindings::ConstantFunction<1, 3, 3>::bind(m);

  Dune::XT::Functions::bindings::ConstantFunction<2, 1, 1>::bind(m);
  Dune::XT::Functions::bindings::ConstantFunction<2, 2, 1>::bind(m);
  Dune::XT::Functions::bindings::ConstantFunction<2, 2, 2>::bind(m);
  Dune::XT::Functions::bindings::ConstantFunction<2, 3, 1>::bind(m);
  Dune::XT::Functions::bindings::ConstantFunction<2, 3, 3>::bind(m);

  Dune::XT::Functions::bindings::ConstantFunction<3, 1, 1>::bind(m);
  Dune::XT::Functions::bindings::ConstantFunction<3, 2, 1>::bind(m);
  Dune::XT::Functions::bindings::ConstantFunction<3, 2, 2>::bind(m);
  Dune::XT::Functions::bindings::ConstantFunction<3, 3, 1>::bind(m);
  Dune::XT::Functions::bindings::ConstantFunction<3, 3, 3>::bind(m);
} // PYBIND11_MODULE(...)
