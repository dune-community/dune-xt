// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2019)
//   Tobias Leibner  (2019 - 2020)

#include "config.h"

#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

#include "spe10.hh"


template <class G>
void bind_model1_for_grid(pybind11::module& m)
{
  using namespace Dune::XT::Functions;
  const auto grid_id = Dune::XT::Grid::bindings::grid_name<G>::value();
  const auto g_dim = G::dimension;

  bind_Spe10Model1Function<G, g_dim>(m, grid_id);
} // ... bind_model1_for_grid(...)


template <class Tuple = Dune::XT::Grid::bindings::Available2dGridTypes>
void bind_model1_for_all_grids(pybind11::module& m)
{
  Dune::XT::Common::bindings::guarded_bind([&]() { //  different grids but same entity
    bind_model1_for_grid<typename Dune::XT::Common::list_content<Tuple>::head>(m);
  });
  bind_model1_for_all_grids<typename Dune::XT::Common::list_content<Tuple>::tail>(m);
}

template <>
void bind_model1_for_all_grids<Dune::XT::Common::tuple_null_type>(pybind11::module&)
{}


template <class G>
void bind_model2_for_grid(pybind11::module& m)
{
  Dune::XT::Functions::bind_Spe10Model2Function<G>(m, Dune::XT::Grid::bindings::grid_name<G>::value());
}


template <class Tuple = Dune::XT::Grid::bindings::Available3dGridTypes>
void bind_model2_for_all_grids(pybind11::module& m)
{
  Dune::XT::Common::bindings::guarded_bind([&]() { //  different grids but same entity
    bind_model2_for_grid<typename Dune::XT::Common::list_content<Tuple>::head>(m);
  });
  bind_model2_for_all_grids<typename Dune::XT::Common::list_content<Tuple>::tail>(m);
}

template <>
void bind_model2_for_all_grids<Dune::XT::Common::tuple_null_type>(pybind11::module&)
{}


PYBIND11_MODULE(_functions_spe10, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_1d");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_2d");
  py::module::import("dune.xt.functions._functions_interfaces_grid_function_3d");

  bind_model1_for_all_grids(m);
  bind_model2_for_all_grids(m);
}
