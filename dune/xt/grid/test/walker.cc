// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2016)
//   Rene Milk       (2014 - 2016)
//   Tobias Leibner  (2015 - 2016)

#include <dune/xt/common/test/main.hxx>

#if DUNE_VERSION_NEWER(DUNE_COMMON, 3, 9) && HAVE_TBB // EXADUNE
#include <dune/grid/utility/partitioning/seedlist.hh>
#endif

#include <dune/xt/common/logstreams.hh>
#include <dune/xt/common/parallel/partitioner.hh>

#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/walker.hh>


using namespace Dune::XT::Common;
using namespace Dune::XT::Grid;
using namespace std;

typedef testing::Types<Int<1>, Int<2>, Int<3>> GridDims;

template <class T>
struct GridWalkerTest : public ::testing::Test
{
  static const size_t griddim = T::value;
  static const size_t level = 4;
  typedef Dune::YaspGrid<griddim, Dune::EquidistantOffsetCoordinates<double, griddim>> GridType;
  typedef typename GridType::LeafGridView GridViewType;
  typedef typename Entity<GridViewType>::Type EntityType;
  typedef typename Intersection<GridViewType>::Type IntersectionType;
  const GridProvider<GridType> grid_prv;
  GridWalkerTest()
    : grid_prv(make_cube_grid<GridType>(0.f, 1.f, level))
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
    auto test0 = [&] {
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
    auto all_counter = [&](const IntersectionType&, const EntityType&, const EntityType&) { all_count++; };

    auto on_filter_boundaries = new ApplyOn::FilteredIntersections<GridViewType>(boundaries);
    auto on_all_boundaries = new ApplyOn::BoundaryIntersections<GridViewType>();
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
