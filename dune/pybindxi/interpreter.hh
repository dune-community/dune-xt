// This file is part of the dune-pybindxi project:
//   https://github.com/dune-community/dune-pybindxi
// Copyright 2009-2017 dune-pybindxi developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_PYBINDXI_INTERPRETER_HH
#define DUNE_PYBINDXI_INTERPRETER_HH

#include <map>
#include <string>

#include <dune/common/visibility.hh>

#include "embed.h"

namespace Dune {
namespace PybindXI {


/**
 * \note Most likely, you do not want to use this class directly, but
 * GlobalInterpreter instead!
 */
class DUNE_EXPORT ScopedInterpreter
{
public:
  pybind11::module import_module(const std::string& module_name);

private:
  pybind11::scoped_interpreter interpreter_;
  std::map<std::string, pybind11::module> modules_;
}; // class ScopedInterpreter


DUNE_EXPORT ScopedInterpreter& GlobalInterpreter();


} // namespace PybindXI
} // namespace Dune

#endif // DUNE_PYBINDXI_INTERPRETER_HH
