// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2019 - 2020)

#include "config.h"

#include <string>
#include <vector>

#include <boost/tuple/tuple.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/grids.hh>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/grid/grids.bindings.hh>


template <class Tuple = Dune::XT::Grid::AvailableGridTypes>
void addbind_for_Grid(pybind11::module& m, std::vector<std::string>& available_types)
{
  using G = typename Tuple::head_type;
  const auto grid_id = Dune::XT::Grid::bindings::grid_name<G>::value();
  available_types.push_back(grid_id);
  addbind_for_Grid<typename Tuple::tail_type>(m, available_types);
} // ... addbind_for_Grid(...)


template <>
void addbind_for_Grid<boost::tuples::null_type>(pybind11::module&, std::vector<std::string>&)
{}

PYBIND11_MODULE(_types, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");

  std::vector<std::string> available_types;
  addbind_for_Grid(m, available_types);
  m.attr("available_types") = available_types;
}
