// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk       (2018)

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/grid/gridprovider/eoc.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/grid/uggrid.hh>

{% for name, type in config.all_grids %}

#if HAVE_MPI
GTEST_TEST(GridProvider_{{name}}, layers)
{
  using ProviderFactory = Dune::XT::Grid::GridProviderFactory<{{type}}>;
  const auto type = ProviderFactory::available()[0];
  const auto config = ProviderFactory::default_config(type);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm split_comm;
  MPI_Comm_split(MPI_COMM_WORLD, world_rank, world_rank, &split_comm);

  {% if "cube" in type or "UGGrid" in type %}
    // uggrid + alu cube grids cannot handle non-default comms
    EXPECT_THROW(ProviderFactory::create(config, split_comm), Dune::InvalidStateException);
    return;
  {% endif %}
  auto local_provider = ProviderFactory::create(config, split_comm);
  auto global_provider = ProviderFactory::create(config, MPI_COMM_WORLD);
  auto& local_grid = local_provider.grid();
  auto& global_grid = global_provider.grid();
  auto& local_colcom = local_grid.comm();
  auto& global_colcom = global_grid.comm();

  const bool actually_parallel = global_grid.comm().size() > 1;
  if(actually_parallel){
    EXPECT_NE(local_colcom.size(), global_colcom.size());
    EXPECT_GE(local_grid.size(0), global_grid.size(0));
  }
  else {
    EXPECT_EQ(local_colcom.size(), global_colcom.size());
    EXPECT_EQ(global_grid.size(0), local_grid.size(0));
  }
}
#else
GTEST_TEST(GridProvider_{{name}}, layers)
{
  using ProviderFactory = Dune::XT::Grid::GridProviderFactory<{{type}}>;
  const auto type = ProviderFactory::available()[0];
  const auto config = ProviderFactory::default_config(type);
  auto local_provider = ProviderFactory::create(config, Dune::MPIHelper::getLocalCommunicator());
  auto global_provider = ProviderFactory::create(config, Dune::MPIHelper::getCommunicator());
  auto& local_grid = local_provider.grid();
  auto& global_grid = global_provider.grid();
  auto& local_colcom = local_grid.comm();
  auto& global_colcom = global_grid.comm();

  EXPECT_EQ(local_colcom.size(), global_colcom.size());
  EXPECT_EQ(global_grid.size(0), local_grid.size(0));
}
#endif
{% endfor %}
