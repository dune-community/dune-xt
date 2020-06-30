// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017, 2019)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2017, 2020)

#include <config.h>

#include <dune/xt/common/logging.hh>
#include <dune/xt/common/string.hh>

#include "python.hh"

namespace Dune {
namespace XT {
namespace Common {
namespace bindings {


void guarded_bind(const std::function<void()>& registrar)
{
  try {
    registrar();
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
} // ... guarded_bind(...)


void add_initialization(pybind11::module& /*m*/, std::string /*logger_name*/, std::string /*so_name*/) {}


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
