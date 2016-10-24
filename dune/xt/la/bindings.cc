// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <string>
#include <vector>

#include <dune/pybindxi/pybind11.h>
//#include <dune/pybindxi/stl_bind.h> // <- see dune/xt/common/bindings.cc

#include <dune/xt/la/container/vector-interface.pbh>

#include <dune/xt/la/container.hh>


namespace py = pybind11;


PYBIND11_PLUGIN(la)
{
  py::module m("la", "dune-xt-la");

  py::module::import("common");

  Dune::XT::LA::bind_Backends(m);
  Dune::XT::LA::bind_Vector<Dune::XT::LA::CommonDenseVector<double>>(m, "CommonDenseVector_double");

  return m.ptr();
} // PYBIND11_PLUGIN(la)

#endif // HAVE_DUNE_PYBINDXI
