// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "main.hxx"

#if HAVE_DUNE_GRID

#include <dune/stuff/grid/information.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/math.hh>
#include <dune/stuff/common/logstreams.hh>
#include <dune/common/shared_ptr.hh>

using namespace Dune::Stuff;
using namespace Dune::Stuff::Common;
using namespace Dune::Stuff::Grid;
using namespace std;

typedef testing::Types<Int<1>, Int<2>, Int<3>> GridDims;

template <class T>
struct GridInfoTest : public ::testing::Test
{
  static const size_t griddim     = T::value;
  static const unsigned int level = 1;
  typedef Dune::YaspGrid<griddim, Dune::EquidistantOffsetCoordinates<double, griddim>> GridType;
  typedef Dimensions<typename GridType::LeafGridView> DimensionsType;

  const DSG::Providers::Cube<GridType> grid_prv;
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
    for (auto i : valueRange(griddim)) {
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

TYPED_TEST_CASE(GridInfoTest, GridDims);
TYPED_TEST(GridInfoTest, Misc)
{
  this->check();
  this->print(dev_null);
}

#endif // #if HAVE_DUNE_GRID
