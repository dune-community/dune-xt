// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2017)
//   Kirsten Weber   (2013)
//   René Fritze     (2013 - 2018)
//   Tobias Leibner  (2014, 2017 - 2018)

#ifndef DUNE_XT_GRID_PROVIDER_HH
#define DUNE_XT_GRID_PROVIDER_HH

#include <memory>
#include <type_traits>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/parallel/communicator.hh>
#include <dune/xt/common/print.hh>
#include <dune/xt/common/type_traits.hh>

#include "cube.hh"
#include "dgf.hh"
#include "gmsh.hh"

namespace Dune {
namespace XT {
namespace Grid {


template <class GridType>
class GridProviderFactory
{
  template <class G>
  static GridProvider<GridType, none_t> call_create(const Common::Configuration& cfg,
                                                    MPIHelper::MPICommunicator mpi_comm)
  {
    if (cfg.empty())
      return G::create(G::default_config(), mpi_comm);
    else
      return G::create(cfg, mpi_comm);
  }

  static std::string available_as_str()
  {
    std::stringstream ret;
    std::copy(available().begin(), available().end(), Common::PrefixOutputIterator<std::string>(ret, "\n   "));
    return ret.str();
  } // ... available_as_str(...)

  typedef CubeGridProviderFactory<GridType> CubeType;
  typedef DgfGridProviderFactory<GridType> DgfType;
  typedef GmshGridProviderFactory<GridType> GmshType;

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
    else if (DgfType::static_id() == type)
      return DgfType::default_config();
    else if (GmshType::static_id() == type)
      return GmshType::default_config();
    else if (available().empty())
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "There is no grid provider available for " << Common::Typename<GridType>::value() << "!");
    else
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "Requested type '" << type << "' is not one of those avaible for "
                                    << Common::Typename<GridType>::value() << ":\n"
                                    << available_as_str());
  } // ... default_config(...)

  static GridProvider<GridType, none_t> create(const Common::Configuration& config,
                                               MPIHelper::MPICommunicator mpi_comm = MPIHelper::getCommunicator())
  {
    return create(config.get<std::string>("type"), config, mpi_comm);
  }

  static GridProvider<GridType, none_t> create(const std::string& type = available()[0],
                                               const Common::Configuration config = Common::Configuration(),
                                               MPIHelper::MPICommunicator mpi_comm = MPIHelper::getCommunicator())
  {
    if (is_alugrid<GridType>::value && mpi_comm != MPIHelper::getCommunicator()) {
      DUNE_THROW(InvalidStateException,
                 "Alugrid either ignores, or outright fails with non-world comes when used here");
    }
    if (CubeType::static_id() == type)
      return call_create<CubeType>(config, mpi_comm);
    else if (DgfType::static_id() == type)
      return call_create<DgfType>(config, mpi_comm);
    else if (GmshType::static_id() == type)
      return call_create<GmshType>(config, mpi_comm);
    else if (available().empty())
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "There is no grid provider available for " << Common::Typename<GridType>::value() << "!");
    else
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "Requested type '" << type << "' is not one of those avaible for "
                                    << Common::Typename<GridType>::value() << ":\n"
                                    << available_as_str());
  } // ... create(...)
}; // class GridProviderFactory


template <class GridType>
class DdSubdomainGridProviderFactory
{
  template <class G>
  static GridProvider<GridType, DD::SubdomainGrid<GridType>> call_create(const Common::Configuration& cfg,
                                                                         MPIHelper::MPICommunicator mpi_comm)
  {
    if (cfg.empty())
      return G::create(G::default_config(), mpi_comm);
    else
      return G::create(cfg, mpi_comm);
  }

  static std::string available_as_str()
  {
    std::stringstream ret;
    std::copy(available().begin(), available().end(), Common::PrefixOutputIterator<std::string>(ret, "\n   "));
    return ret.str();
  } // ... available_as_str(...)

  typedef CubeDdSubdomainsGridProviderFactory<GridType> CubeType;

public:
  static std::vector<std::string> available()
  {
    return std::vector<std::string>{CubeType::static_id()};
  }

  static Common::Configuration default_config(const std::string type)
  {
    if (CubeType::static_id() == type)
      return CubeType::default_config();
    else if (available().empty())
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "There is no grid provider available for " << Common::Typename<GridType>::value() << "!");
    else
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "Requested type '" << type << "' is not one of those avaible for "
                                    << Common::Typename<GridType>::value() << ":\n"
                                    << available_as_str());
  } // ... default_config(...)

  static GridProvider<GridType, DD::SubdomainGrid<GridType>>
  create(const Common::Configuration& config, MPIHelper::MPICommunicator mpi_comm = MPIHelper::getCommunicator())
  {
    return create(config.get<std::string>("type"), config, mpi_comm);
  }

  static GridProvider<GridType, DD::SubdomainGrid<GridType>>
  create(const std::string& type = available()[0],
         const Common::Configuration config = Common::Configuration(),
         MPIHelper::MPICommunicator mpi_comm = MPIHelper::getCommunicator())
  {
    if (is_alugrid<GridType>::value && mpi_comm != MPIHelper::getCommunicator()) {
      DUNE_THROW(InvalidStateException,
                 "Alugrid either ignores, or outright fails with non-world comes when used here");
    }

    if (CubeType::static_id() == type)
      return call_create<CubeType>(config, mpi_comm);
    else if (available().empty())
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "There is no grid provider available for " << Common::Typename<GridType>::value() << "!");
    else
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "Requested type '" << type << "' is not one of those avaible for "
                                    << Common::Typename<GridType>::value() << ":\n"
                                    << available_as_str());
  } // ... create(...)
}; // class DdSubdomainGridProviderFactory


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_PROVIDER_HH
