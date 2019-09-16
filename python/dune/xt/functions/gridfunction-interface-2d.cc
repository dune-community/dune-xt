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

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

#include "gridfunction-interface.hh"


template <class Tuple = Dune::XT::Grid::Available2dGridTypes>
void bind_all_2d_grids(pybind11::module& m)
{
  Dune::XT::Common::bindings::try_register(m, [](auto& m_) { //  different grids but same entity
    Dune::XT::Functions::bindings::addbind_GridFunctionInterface_all_dims<typename Tuple::head_type>(m_);
  });
  bind_all_2d_grids<typename Tuple::tail_type>(m);
}

template <>
void bind_all_2d_grids<boost::tuples::null_type>(pybind11::module&)
{}


PYBIND11_MODULE(_gridfunction_interface_2d, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions._function_interface_2d");

  bind_all_2d_grids(m);
} // PYBIND11_MODULE(...)
