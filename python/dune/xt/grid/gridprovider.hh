// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_XT_GRID_GRIDPROVIDER_PBH
#define DUNE_XT_GRID_GRIDPROVIDER_PBH

#include <algorithm>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <dune/xt/common/parallel/mpi_comm_wrapper.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/grid/entity.hh>

#include <dune/xt/grid/gridprovider.hh>

namespace Dune {
namespace XT {
namespace Grid {


template <class G>
void bind_make_cube_grid(pybind11::module& m, const std::string& grid_id)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  m.def(std::string("make_cube_grid__" + grid_id).c_str(),
        [](const Common::Configuration& cfg, Common::MPI_Comm_Wrapper mpi_comm) {
          return make_cube_grid<G>(cfg, mpi_comm.get());
        },
        "cfg"_a = cube_gridprovider_default_config(),
        "mpi_comm"_a = Common::MPI_Comm_Wrapper());

  m.def(std::string("make_cube_grid__" + grid_id).c_str(),
        [](const FieldVector<typename G::ctype, G::dimension>& lower_left,
           const FieldVector<typename G::ctype, G::dimension>& upper_right,
           const std::array<unsigned int, G::dimension>& num_elements,
           const unsigned int num_refinements,
           const std::array<unsigned int, G::dimension>& overlap_size,
           Common::MPI_Comm_Wrapper mpi_comm) {
          return make_cube_grid<G>(
              lower_left, upper_right, num_elements, num_refinements, overlap_size, mpi_comm.get());
        },
        "lower_left"_a,
        "upper_right"_a,
        "num_elements"_a = XT::Common::make_array<unsigned int, G::dimension>(
            cube_gridprovider_default_config().template get<std::vector<unsigned int>>("num_elements")),
        "num_refinements"_a = cube_gridprovider_default_config().template get<unsigned int>("num_refinements"),
        "overlap_size"_a = XT::Common::make_array<unsigned int, G::dimension>(
            cube_gridprovider_default_config().template get<std::vector<unsigned int>>("overlap_size")),
        "mpi_comm"_a = Common::MPI_Comm_Wrapper());
} // ... bind_make_cube_grid(...)


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_GRIDPROVIDER_PBH
