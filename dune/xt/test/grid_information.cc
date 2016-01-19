// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2012 - 2015)
//   Rene Milk       (2012 - 2015)
//   Tobias Leibner  (2014)

#include <dune/xt/test/main.hxx>

#if HAVE_DUNE_GRID

#include <dune/common/shared_ptr.hh>

#include <dune/xt/common/logstreams.hh>
#include <dune/xt/common/math.hh>
#include <dune/xt/common/ranges.hh>

#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/provider/cube.hh>

using namespace std;
using namespace Dune::XT;
using namespace Dune::XT::Common;
using namespace Dune::XT::Grid;

struct GridInfoTest : public ::testing::Test
{
  static const size_t griddim     = TESTGRIDDIM;
  static const unsigned int level = 1;
  typedef Dune::YaspGrid<griddim, Dune::EquidistantOffsetCoordinates<double, griddim>> GridType;
  typedef Dimensions<typename GridType::LeafGridView> DimensionsType;

  const Providers::Cube<GridType> grid_prv;
  GridInfoTest()
    : grid_prv(0.f, 1.f, level)
  {
  }

  void check_dimensions(const DimensionsType& dim, const size_t entities)
  {
    EXPECT_DOUBLE_EQ(1.0 / double(entities), dim.entity_volume.min());
    EXPECT_DOUBLE_EQ(dim.entity_volume.min(), dim.entity_volume.max());
    EXPECT_DOUBLE_EQ(dim.entity_volume.min(), dim.entity_volume.average());
    EXPECT_DOUBLE_EQ(1.0, dim.volumeRelation());
    const auto& dl = dim.coord_limits;
    for (auto i : value_range(griddim)) {
      EXPECT_DOUBLE_EQ(dl[i].max(), 1.0);
      EXPECT_DOUBLE_EQ(dl[i].min(), 0.0);
      EXPECT_DOUBLE_EQ(dl[i].average(), 0.5);
    }
  }

  void check()
  {
    const auto gv       = grid_prv.grid().leafGridView();
    const auto entities = gv.size(0);
    check_dimensions(DimensionsType(grid_prv.grid().leafGridView()), entities);
    const auto& first_entity = *(grid_prv.grid().leafGridView().template begin<0>());
    check_dimensions(DimensionsType(first_entity), 1u);
    const Statistics st(gv);
    const auto line = std::pow(2, level);
    EXPECT_EQ(line * (griddim), st.numberOfBoundaryIntersections);
    EXPECT_EQ(entities * (2 * griddim), st.numberOfIntersections);
    EXPECT_EQ(st.numberOfIntersections - st.numberOfBoundaryIntersections, st.numberOfInnerIntersections);
    EXPECT_EQ(griddim * 2, maxNumberOfNeighbors(gv));
  }

  void print(std::ostream& out)
  {
    const auto& gv = grid_prv.grid().leafGridView();
    printInfo(gv, out);
  }
};

TEST_F(GridInfoTest, Misc)
{
  this->check();
  this->print(dev_null);
}

#endif // #if HAVE_DUNE_GRID
