// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2017, 2019)
//   Kirsten Weber   (2013)
//   René Fritze     (2013 - 2020)
//   Tobias Leibner  (2014, 2017 - 2020)

#ifndef DUNE_XT_GRID_PROVIDER_HH
#define DUNE_XT_GRID_PROVIDER_HH

#include <memory>
#include <type_traits>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/print.hh>
#include <dune/xt/common/type_traits.hh>

#include "cube.hh"
#include "dgf.hh"
#include "gmsh.hh"

namespace Dune::XT::Grid {


template <class GridType>
class GridProviderFactory
{
  template <class G>
  static GridProvider<GridType> call_create(const Common::Configuration& cfg, MPIHelper::MPICommunicator mpi_comm)
  {
    if (cfg.empty())
      return G::create(G::default_config(), mpi_comm);
    return G::create(cfg, mpi_comm);
  }

  static std::string available_as_str()
  {
    std::stringstream ret;
    std::copy(available().begin(), available().end(), Common::PrefixOutputIterator<std::string>(ret, "\n   "));
    return ret.str();
  } // ... available_as_str(...)

  using CubeType = CubeGridProviderFactory<GridType>;
  using DgfType = DgfGridProviderFactory<GridType>;
  using GmshType = GmshGridProviderFactory<GridType>;

public:
  static std::vector<std::string> available()
  {
    std::vector<std::string> ret{CubeType::static_id(), DgfType::static_id()};
    if (GmshType::available)
      ret.push_back(GmshType::static_id());
    return ret;
  }

  static Common::Configuration default_config(const std::string type)
  {
    if (CubeType::static_id() == type)
      return CubeType::default_config();
    if (DgfType::static_id() == type)
      return DgfType::default_config();
    if (GmshType::static_id() == type)
      return GmshType::default_config();
    if (available().empty())
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "There is no grid provider available for " << Common::Typename<GridType>::value() << "!");
    else
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "Requested type '" << type << "' is not one of those avaible for "
                                    << Common::Typename<GridType>::value() << ":\n"
                                    << available_as_str());
  } // ... default_config(...)

  static GridProvider<GridType> create(const Common::Configuration& config,
                                       MPIHelper::MPICommunicator mpi_comm = MPIHelper::getCommunicator())
  {
    return create(config.get<std::string>("type"), config, mpi_comm);
  }

  static GridProvider<GridType> create(const std::string& type = available()[0],
                                       const Common::Configuration config = Common::Configuration(),
                                       MPIHelper::MPICommunicator mpi_comm = MPIHelper::getCommunicator())
  {
    if (is_alugrid<GridType>::value && mpi_comm != MPIHelper::getCommunicator()) {
      DUNE_THROW(InvalidStateException,
                 "Alugrid either ignores, or outright fails with non-world comes when used here");
    }
    if (CubeType::static_id() == type)
      return call_create<CubeType>(config, mpi_comm);
    if (DgfType::static_id() == type)
      return call_create<DgfType>(config, mpi_comm);
    if (GmshType::static_id() == type)
      return call_create<GmshType>(config, mpi_comm);
    if (available().empty())
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "There is no grid provider available for " << Common::Typename<GridType>::value() << "!");
    else
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "Requested type '" << type << "' is not one of those avaible for "
                                    << Common::Typename<GridType>::value() << ":\n"
                                    << available_as_str());
  } // ... create(...)
}; // class GridProviderFactory


} // namespace Dune::XT::Grid

#endif // DUNE_XT_GRID_PROVIDER_HH
