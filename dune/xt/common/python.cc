// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018)
//   Tobias Leibner  (2017)

#include <config.h>

#include "python.hh"


#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/xt/common/math.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/common/string.hh>

namespace Dune {
namespace XT {
namespace Common {
namespace bindings {


void add_initialization(pybind11::module& m, std::string logger_name, std::string so_name)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  m.def("_init_mpi",
        [](const std::vector<std::string>& args) {
          int argc = Dune::XT::Common::numeric_cast<int>(args.size());
          char** argv = Dune::XT::Common::vector_to_main_args(args);
          Dune::MPIHelper::instance(argc, argv);
        },
        "args"_a = std::vector<std::string>());

  m.def(
      "_init_logger",
      [](const ssize_t max_info_level,
         const ssize_t max_debug_level,
         const bool enable_warnings,
         const bool enable_colors,
         const std::string& info_color,
         const std::string& debug_color,
         const std::string& warning_color) {
        try {
          Dune::XT::Common::TimedLogger().create(
              max_info_level, max_debug_level, enable_warnings, enable_colors, info_color, debug_color, warning_color);
        } catch (Dune::XT::Common::Exceptions::logger_error&) {
        }
      },
      "max_info_level"_a = std::numeric_limits<ssize_t>::max(),
      "max_debug_level"_a = std::numeric_limits<ssize_t>::max(),
      "enable_warnings"_a = true,
      "enable_colors"_a = true,
      "info_color"_a = "blue",
      "debug_color"_a = "darkgray",
      "warning_color"_a = "red");

  m.def("_test_logger",
        [=](const bool info, const bool debug, const bool warning) {
          auto logger = Dune::XT::Common::TimedLogger().get(logger_name);
          auto prefix = so_name.empty() ? "" : "(from " + so_name + ".so) ";
          if (info)
            logger.info() << prefix << "info logging works!" << std::endl;
          if (debug)
            logger.debug() << prefix << "debug logging works!" << std::endl;
          if (warning)
            logger.warn() << prefix << "warning logging works!" << std::endl;
        },
        "info"_a = true,
        "debug"_a = true,
        "warning"_a = true);
} // ... add_initialization(...)


void try_register(pybind11::module& m, const std::function<void(pybind11::module&)>& registrar)
{
  try {
    registrar(m);
  } catch (std::runtime_error& err) {
    const std::string what{err.what()};
    /*  pybind11 error msg format
     * ("generic_type: type \"" + std::string(rec->name) +
                          "\" is already registered!")
     */
    const auto reg_pos = what.find("already registered");
    const auto def_pos = what.find("already defined");
    size_t npos{std::string::npos}, left{0}, right{std::string::npos};
    if ((def_pos == npos && reg_pos == npos) || what.size() < 2)
      throw err;

    std::string type{};
    left = what.find("\"");
    right = what.rfind("\"");
    type = what.substr(left + 1, right - left - 1);
    DXTC_LOG_DEBUG << "try_register: added type " << type << std::endl;
  }
} // ... try_register(...)


} // namespace bindings
} // namespace Common
} // namespace XT
} // namespace Dune
