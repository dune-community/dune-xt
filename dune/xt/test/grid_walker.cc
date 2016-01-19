// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015)
//   Rene Milk       (2014 - 2015)
//   Tobias Leibner  (2015)

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
  static const size_t griddim = T::value;
  static const size_t level   = 4;
  typedef Dune::YaspGrid<griddim, Dune::EquidistantOffsetCoordinates<double, griddim>> GridType;
  typedef typename GridType::LeafGridView GridViewType;
  typedef typename DSG::Entity<GridViewType>::Type EntityType;
  typedef typename DSG::Intersection<GridViewType>::Type IntersectionType;
  const DSG::Providers::Cube<GridType> grid_prv;
  GridWalkerTest()
    : grid_prv(0.f, 1.f, level)
  {
  }

  void check_count()
  {
    const auto gv = grid_prv.grid().leafGridView();
    Walker<GridViewType> walker(gv);
    const auto correct_size = gv.size(0);
    atomic<size_t> count(0);
    auto counter = [&](const EntityType&) { count++; };
    auto test1 = [&] {
      walker.add(counter);
      walker.walk(false);
    };
    auto test2 = [&] {
      walker.add(counter);
      walker.walk(true);
    };
    auto test3 = [&] { walker.add(counter).walk(true); };
    list<function<void()>> tests({test1, test2, test3});
#if DUNE_VERSION_NEWER(DUNE_COMMON, 3, 9) && HAVE_TBB // EXADUNE
    auto test0        = [&] {
      const auto& set = gv.grid().leafIndexSet();
      IndexSetPartitioner<GridViewType> partitioner(set);
      EXPECT_EQ(set.size(0), partitioner.partitions());
      Dune::SeedListPartitioning<GridType, 0> partitioning(gv, partitioner);
      walker.add(counter);
      walker.walk(partitioning);
    };
    tests.push_back(test0);
#endif // DUNE_VERSION_NEWER(DUNE_COMMON, 3, 9) && HAVE_TBB

    for (const auto& test : tests) {
      count = 0;
      test();
      EXPECT_EQ(count, correct_size);
    }
  }

  void check_apply_on()
  {
    const auto gv = grid_prv.grid().leafGridView();
    Walker<GridViewType> walker(gv);

    size_t filter_count = 0, all_count = 0;
    auto boundaries = [=](const GridViewType&, const IntersectionType& inter) { return inter.boundary(); };
    auto filter_counter = [&](const IntersectionType&, const EntityType&, const EntityType&) { filter_count++; };
    auto all_counter    = [&](const IntersectionType&, const EntityType&, const EntityType&) { all_count++; };

    auto on_filter_boundaries = new DSG::ApplyOn::FilteredIntersections<GridViewType>(boundaries);
    auto on_all_boundaries = new DSG::ApplyOn::BoundaryIntersections<GridViewType>();
    walker.add(filter_counter, on_filter_boundaries);
    walker.add(all_counter, on_all_boundaries);
    walker.walk();
    EXPECT_EQ(filter_count, all_count);
  }
};

TYPED_TEST_CASE(GridWalkerTest, GridDims);
TYPED_TEST(GridWalkerTest, Misc)
{
  this->check_count();
  this->check_apply_on();
}

#else // HAVE_DUNE_GRID

TEST(DISABLED_GridWalkerTest, Misc){};

#endif // HAVE_DUNE_GRID
