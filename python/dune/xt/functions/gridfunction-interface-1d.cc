// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2018)
//   Tim Keil        (2018)

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


template <class Tuple = Dune::XT::Grid::Available1dGridTypes>
void bind_all_1d_grids(pybind11::module& m)
{
  Dune::XT::Common::bindings::try_register(m, [](auto& m_) { //  different grids but same entity
    Dune::XT::Functions::bindings::addbind_GridFunctionInterface_all_dims<typename Tuple::head_type>(m_);
  });
  bind_all_1d_grids<typename Tuple::tail_type>(m);
}

template <>
void bind_all_1d_grids<boost::tuples::null_type>(pybind11::module&)
{}


PYBIND11_MODULE(_gridfunction_interface_1d, m)
{
  Dune::XT::Common::bindings::addbind_exceptions(m);
  Dune::XT::Common::bindings::add_initialization(m, "dune.xt.functions");

  bind_all_1d_grids(m);
} // PYBIND11_MODULE(...)
