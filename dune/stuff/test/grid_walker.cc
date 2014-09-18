// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "main.hxx"

#if HAVE_DUNE_GRID

#include <dune/stuff/grid/walker.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/common/parallel/partitioner.hh>
#include <dune/stuff/common/logstreams.hh>

#if DUNE_VERSION_NEWER(DUNE_COMMON, 3, 9) && HAVE_TBB // EXADUNE
#include <dune/grid/utility/partitioning/seedlist.hh>
#endif

using namespace Dune::Stuff;
using namespace Dune::Stuff::Common;
using namespace Dune::Stuff::Grid;
using namespace std;

typedef testing::Types<Int<1>, Int<2>, Int<3>> GridDims;

template <class T>
struct GridWalkerTest : public ::testing::Test
{
  static const int griddim        = T::value;
  static const unsigned int level = 1;
  typedef Dune::SGrid<griddim, griddim> GridType;
  typedef typename GridType::LeafGridView GridViewType;
  typedef typename GridType::template Codim<0>::Entity EntityType;
  const DSG::Providers::Cube<GridType> grid_prv;
  GridWalkerTest()
    : grid_prv(0.f, 1.f, level)
  {
  }

  void check()
  {
    const auto gv = grid_prv.grid().leafGridView();
    Walker<GridViewType> walker(gv);
    size_t count = 0;
    auto counter = [&](const EntityType&) { count++; };
    walker.add(counter);
    walker.walk(false);
    EXPECT_EQ(count, gv.size(0));
#if DUNE_VERSION_NEWER(DUNE_COMMON, 3, 9) && HAVE_TBB // EXADUNE
    count                 = 0;
    const auto& index_set = gv.grid().leafIndexSet();
    IndexSetPartitioner<GridViewType> partitioner(index_set);
    Dune::SeedListPartitioning<GridType, 0> partitioning(gv, partitioner);
    walker.tbb_walk(partitioning);
    EXPECT_EQ(count, gv.size(0));
#endif // DUNE_VERSION_NEWER(DUNE_COMMON,3,9) && HAVE_TBB
  }
};

TYPED_TEST_CASE(GridWalkerTest, GridDims);
TYPED_TEST(GridWalkerTest, Misc)
{
  this->check();
}

#endif // #if HAVE_DUNE_GRID
