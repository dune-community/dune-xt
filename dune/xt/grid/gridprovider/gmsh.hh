// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2012 - 2017, 2019)
//   Kirsten Weber   (2012)
//   René Fritze     (2012 - 2013, 2015 - 2016, 2018 - 2019)
//   Tobias Leibner  (2014 - 2016, 2020)

#ifndef DUNE_XT_GRID_GRIDPROVIDER_GMSH_HH
#define DUNE_XT_GRID_GRIDPROVIDER_GMSH_HH

#include <memory>
#include <type_traits>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>

#include <dune/xt/grid/grids.hh>

#include "provider.hh"

namespace Dune {
namespace XT {
namespace Grid {


static inline std::string gmsh_gridprovider_id()
{
  return "xt.grid.gridprovider.gmsh";
}


static inline Common::Configuration gmsh_gridprovider_default_config()
{
  Common::Configuration config;
  config["type"] = gmsh_gridprovider_id();
  config["filename"] = "g.msh";
  return config;
}


template <class GridType>
class GmshGridProviderFactory;

//! Gmsh grid creation relies on unstructured grid factory -> disable yasp
template <int dim, class Coordinates>
class GmshGridProviderFactory<Dune::YaspGrid<dim, Coordinates>>
{
public:
  static constexpr bool available = false;

  static std::string static_id()
  {
    return "dune.xt.grid.gmsh";
  }

  static Common::Configuration default_config()
  {
    DUNE_THROW(NotImplemented, "No default config for GmshGridProviderFactory Yaspgrid");
  }

  template <class... Args>
  static GridProvider<Dune::YaspGrid<dim, Coordinates>> create(Args&&...)
  {
    DUNE_THROW(NotImplemented, "No create implemeted for GmshGridProviderFactory Yaspgrid");
  }
}; // class GmshGridProviderFactory<Dune::YaspGrid<...>>


template <class GridType>
class GmshGridProviderFactory
{
public:
  static constexpr bool available = false;

  static std::string static_id()
  {
    return gmsh_gridprovider_id();
  }

  static Common::Configuration default_config()
  {
    auto cfg = gmsh_gridprovider_default_config();
#if HAVE_DUNE_ALUGRID
    if (std::is_same<ALUGrid<2, 2, simplex, conforming>, GridType>::value
        || std::is_same<ALUGrid<2, 2, simplex, nonconforming>, GridType>::value) {
      cfg["filename"] = "gmsh_2d_simplices.msh";
    }
#endif // HAVE_DUNE_ALUGRID
    return cfg;
  }

  static GridProvider<GridType> create(const std::string& filename, MPIHelper::MPICommunicator mpi_comm)
  {
    DUNE_THROW_IF(CollectiveCommunication<MPIHelper::MPICommunicator>(mpi_comm).size() > 1,
                  InvalidStateException,
                  "Gmsh reading not implemented in parallel");
    return GridProvider<GridType>(GmshReader<GridType>::read(filename));
  }

  static GridProvider<GridType> create(const Common::Configuration& cfg, MPIHelper::MPICommunicator mpi_comm)
  {
    return create(cfg.get("filename", default_config().template get<std::string>("filename")), mpi_comm);
  }
}; // class GmshGridProviderFactory


template <class GridType>
auto make_gmsh_grid(const std::string& filename, MPIHelper::MPICommunicator mpi_comm = MPIHelper::getCommunicator())
{
  static_assert(is_grid<GridType>::value);
  return GmshGridProviderFactory<GridType>::create(filename, mpi_comm);
}


template <class GridType>
auto make_gmsh_grid(const Common::Configuration& cfg = GmshGridProviderFactory<GridType>::default_config(),
                    MPIHelper::MPICommunicator mpi_comm = MPIHelper::getCommunicator())
{
  static_assert(is_grid<GridType>::value);
  return GmshGridProviderFactory<GridType>::create(cfg, mpi_comm);
}


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_GRIDPROVIDER_GMSH_HH
