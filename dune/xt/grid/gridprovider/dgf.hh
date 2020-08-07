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
//   Tobias Leibner  (2014, 2016, 2020)

#ifndef DUNE_XT_GRID_GRIDPROVIDER_DGF_HH
#define DUNE_XT_GRID_GRIDPROVIDER_DGF_HH

#include <memory>

#include <dune/grid/io/file/dgfparser/dgfgridfactory.hh> // How convenient that GridPtr requires DGFGridFactory but ...
#include <dune/grid/io/file/dgfparser.hh>
#if HAVE_DUNE_ALUGRID
#  include <dune/alugrid/dgf.hh>
#endif
#include <dune/grid/io/file/dgfparser/gridptr.hh> // ... does not include it!

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/string.hh>

#include "provider.hh"

namespace Dune::XT::Grid {


static inline std::string dgf_gridprovider_id()
{
  return "xt.grid.gridprovider.dgf";
}


static inline Common::Configuration dgf_gridprovider_default_config()
{
  Common::Configuration config;
  config["type"] = dgf_gridprovider_id();
  config["filename"] = "dgf_1d_interval.dgf";
  return config;
}


template <class GridType>
class DgfGridProviderFactory
{
  static_assert(is_grid<GridType>::value);

public:
  static const bool available = true;

  static std::string static_id()
  {
    return dgf_gridprovider_id();
  }

  static Common::Configuration default_config()
  {
    auto cfg = dgf_gridprovider_default_config();
    cfg["filename"] = std::string("dgf_") + Common::to_string(size_t(GridType::dimension)) + "d_interval.dgf";
    return cfg;
  }

  static GridProvider<GridType> create(const std::string& filename, MPIHelper::MPICommunicator mpi_comm)
  {
    return GridProvider<GridType>(GridPtr<GridType>(filename, mpi_comm).release());
  }

  static GridProvider<GridType> create(const Common::Configuration& cfg, MPIHelper::MPICommunicator mpi_comm)
  {
    return create(cfg.get("filename", default_config().template get<std::string>("filename")), mpi_comm);
  }
}; // class DgfGridProviderFactory


template <class GridType>
auto make_dgf_grid(const std::string& filename, MPIHelper::MPICommunicator mpi_comm = MPIHelper::getCommunicator())
{
  static_assert(is_grid<GridType>::value);
  return DgfGridProviderFactory<GridType>(filename, mpi_comm);
}


template <class GridType>
auto make_dgf_grid(const Common::Configuration& cfg = DgfGridProviderFactory<GridType>::default_config(),
                   MPIHelper::MPICommunicator mpi_comm = MPIHelper::getCommunicator())
{
  static_assert(is_grid<GridType>::value);
  return DgfGridProviderFactory<GridType>::create(cfg, mpi_comm);
}


} // namespace Dune::XT::Grid

#endif // DUNE_XT_GRID_GRIDPROVIDER_DGF_HH
