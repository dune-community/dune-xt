// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2016)
//   Rene Milk       (2012 - 2016)
//   Tobias Leibner  (2014 - 2016)

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/grid/gridprovider.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/search.hh>
#include <dune/xt/grid/information.hh>
#include <dune/xt/common/logging.hh>

static inline Dune::XT::Common::Configuration cube_gridprovider_default_config()
{
  auto config = Dune::XT::Grid::cube_gridprovider_default_config();
  config["lower_left"] = "[0 0 0 0]";
  config["upper_right"] = "[1 1 1 1]";
  config["num_elements"] = "[3 3 3 3]";
  return config;
}

struct InLevelSearch : public testing::Test
{
  typedef TESTGRIDTYPE GridType;
  typedef Dune::XT::Grid::GridProvider<GridType> GridProviderType;

  InLevelSearch()
    : grid_provider_(Dune::XT::Grid::make_cube_grid<TESTGRIDTYPE>())
  {
  }

  void check()
  {
    const auto view = grid_provider_.leaf_view();
    const auto dimensions = Dune::XT::Grid::dimensions(grid_provider_.leaf_view());
    Dune::XT::Grid::EntityInlevelSearch<decltype(view)> search(view);
    const auto center = dimensions.view_center();
    const auto result = search(std::vector<std::remove_const<decltype(center)>::type>{center});
    EXPECT_GE(result.size(), 1);
  }

  GridProviderType grid_provider_;
}; // class ConstGridProviderBase


TEST_F(InLevelSearch, check)
{
  this->check();
}
