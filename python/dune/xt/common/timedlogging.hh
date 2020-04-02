// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_XT_COMMON_TIMEDLOGGING_HH
#define PYTHON_DUNE_XT_COMMON_TIMEDLOGGING_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/timedlogging.hh>

namespace Dune {
namespace XT {
namespace Common {
namespace bindings {


class DefaultLogger
{
public:
  using type = Common::DefaultLogger;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m, const std::string& class_id = "default_logger")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto ClassName = Common::to_camel_case(class_id);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def_readwrite("info_enabled", &type::info_enabled);
    c.def_readwrite("debug_enabled", &type::debug_enabled);
    c.def_readwrite("warn_enabled", &type::warn_enabled);
    c.def("enable", [](type& self, const std::string& prefix){self.enable(prefix);}, ""_a = "");
    c.def("disable", [](type& self){self.disable();});
    c.def("info", [](type& self, const std::string& to_print){self.info() << to_print << std::endl;});
    c.def("debug", [](type& self, const std::string& to_print){self.debug() << to_print << std::endl;});
    c.def("warn", [](type& self, const std::string& to_print){self.warn() << to_print << std::endl;});

    return c;
  }
}; // class DefaultLogger


} // namespace bindings
} // namespace Common
} // namespace XT
} // namespace Dune

#endif // PYTHON_DUNE_XT_COMMON_TIMEDLOGGING_HH
