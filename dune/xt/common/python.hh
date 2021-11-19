// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2018 - 2020)
//   Tobias Leibner  (2020)
//
// Created by r_milk01 on 4/25/18.

#ifndef DUNE_XT_COMMON_PYTHON_HH
#define DUNE_XT_COMMON_PYTHON_HH

#include <functional>

#include <dune/pybindxi/pybind11.h>

namespace Dune::XT::Common::bindings {


void guarded_bind(const std::function<void()>& registrar);


[[deprecated("This is not required any more (08.08.2019)!")]] void
add_initialization(pybind11::module& /*m*/, const std::string& /*logger_name*/, const std::string& /*so_name*/ = "");


[[deprecated("use guarded_bind() instead (08.08.2019)!")]] void
try_register(pybind11::module& m, const std::function<void(pybind11::module&)>& registrar);


/** Use the python C-API to raise a warning
 * Be careful not to call this from bounded methods that release the GIL!
 * for the call sig see https://docs.python.org/3/c-api/exceptions.html#issuing-warnings
 * for possible categories see https://docs.python.org/3/c-api/exceptions.html#standard-warning-categories
 * example: Common::bindings::warning(PyExc_RuntimeWarning, "some Message");
 */
void warning(PyObject* category, const char* warning, Py_ssize_t stack = 1)
{
  if (PyErr_WarnEx(category, warning, stack) == -1)
    throw new pybind11::error_already_set;
}

} // namespace Dune::XT::Common::bindings

#endif // DUNE_XT_COMMON_PYTHON_HH
