// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019 - 2020)
//   René Fritze     (2019)
//   Tobias Leibner  (2020)

#include "config.h"

#include <exception>

#include <boost/exception/all.hpp>

#include <dune/common/exceptions.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>


PYBIND11_MODULE(_common_exceptions, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  static py::exception<boost::exception> boost_exc(m, "BoostError");
  static py::exception<Dune::Exception> dune_exc(m, "DuneError");

  py::register_exception_translator([](std::exception_ptr p) {
    try {
      if (p)
        std::rethrow_exception(p);
    } catch (const boost::exception& e) {
      boost_exc(boost::diagnostic_information_what(e));
    } catch (const Dune::Exception& e) {
      dune_exc(e.what());
    }
  });
}
