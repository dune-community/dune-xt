// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2017)
//   Rene Milk       (2012 - 2016, 2018)
//   Tobias Leibner  (2014 - 2016)

#include <dune/xt/common/test/main.hxx>

#include <dune/common/shared_ptr.hh>

#include <dune/xt/common/logstreams.hh>
#include <dune/xt/common/math.hh>
#include <dune/xt/common/ranges.hh>

#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

using namespace std;
using namespace Dune::XT;
using namespace Dune::XT::Common;
using namespace Dune::XT::Grid;

struct GridInfoTest : public ::testing::Test
{
  static const size_t griddim = TESTGRIDDIM;
  static const unsigned int level = 1;
  typedef Dune::YaspGrid<griddim, Dune::EquidistantOffsetCoordinates<double, griddim>> GridType;
  typedef Dimensions<typename GridType::LeafGridView> DimensionsType;

  const GridProvider<GridType, none_t> grid_prv;
  GridInfoTest()
    : grid_prv(make_cube_grid<GridType>(0.f, 1.f, level))
  {
  }

  void check_dimensions(const DimensionsType& dim, const size_t entities)
  {
    EXPECT_DOUBLE_EQ(1.0 / double(entities), dim.entity_volume.min());
    EXPECT_DOUBLE_EQ(dim.entity_volume.min(), dim.entity_volume.max());
    EXPECT_DOUBLE_EQ(dim.entity_volume.min(), dim.entity_volume.average());
    EXPECT_DOUBLE_EQ(1.0, dim.volume_relation());
    const auto& dl = dim.coord_limits;
    for (auto i : value_range(griddim)) {
      EXPECT_DOUBLE_EQ(dl[i].max(), 1.0);
      EXPECT_DOUBLE_EQ(dl[i].min(), 0.0);
      EXPECT_DOUBLE_EQ(dl[i].average(), 0.5);
    }
  }

  void check()
  {
    const auto gv = grid_prv.grid().leafGridView();
    const auto entities = gv.size(0);
    check_dimensions(DimensionsType(grid_prv.grid().leafGridView()), entities);
    const auto first_entity = *(grid_prv.grid().leafGridView().template begin<0>());
    check_dimensions(DimensionsType(first_entity), 1u);
    const Statistics st(gv);
    const auto line = std::pow(2, level);
    EXPECT_EQ(line * (griddim), st.numberOfBoundaryIntersections);
    EXPECT_EQ(entities * (2 * griddim), st.numberOfIntersections);
    EXPECT_EQ(st.numberOfIntersections - st.numberOfBoundaryIntersections, st.numberOfInnerIntersections);
    EXPECT_EQ(griddim * 2, max_number_of_neighbors(gv));
  }

  void print(std::ostream& out)
  {
    const auto& gv = grid_prv.grid().leafGridView();
    print_info(gv, out);
  }
};

TEST_F(GridInfoTest, Misc)
{
  this->check();
  this->print(dev_null);
}
