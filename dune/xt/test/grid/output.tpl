// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2019 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk       (2018)

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/grid/gridprovider/eoc.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/output/entity_visualization.hh>
#include <dune/grid/uggrid.hh>

{% for name, type in config.all_grids %}


GTEST_TEST(GridProvider_{{name}}, layers)
{
  using ProviderFactory = Dune::XT::Grid::GridProviderFactory<{{type}}>;
  const auto type = ProviderFactory::available()[0];
  const auto config = ProviderFactory::default_config(type);
  auto global_provider = ProviderFactory::create(config, Dune::MPIHelper::getCommunicator());
  auto& global_grid = global_provider.grid();

  Dune::XT::Grid::ElementVisualization::all(global_grid,  boost::filesystem::temp_directory_path().string());
}


{% endfor %}
