// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2014 - 2018)
//   Tobias Leibner  (2015 - 2016, 2018)

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
  typedef typename GridType::LeafGridView GridLayerType;
  using EntityType = extract_entity_t<GridLayerType>;
  using IntersectionType = extract_intersection_t<GridLayerType>;
  const GridProvider<GridType, none_t> grid_prv;
  GridWalkerTest()
    : grid_prv(make_cube_grid<GridType>(0.f, 1.f, level))
  {
  }

  void check_count()
  {
    const auto gv = grid_prv.grid().leafGridView();
    Walker<GridLayerType> walker(gv);
    const auto correct_size = gv.size(0);
    atomic<size_t> count(0);
    atomic<size_t> intersection_count(0);
    auto counter = [&](const EntityType&) { count++; };
    auto intersection_counter = [&](const IntersectionType&, const EntityType&, const EntityType&) {
      intersection_count++;
    };
    auto test1 = [&] {
      walker.append(counter);
      walker.walk(false);
    };
    auto test2 = [&] {
      walker.append(counter);
      walker.walk(true);
    };
    auto test3 = [&] { walker.append(counter).walk(true); };
    auto test4 = [&] { walker.append(intersection_counter).walk(false); };
    auto test5 = [&] { walker.append(intersection_counter).walk(true); };
    list<function<void()>> tests({test1, test2, test3});
#if DUNE_VERSION_NEWER(DUNE_COMMON, 3, 9) && HAVE_TBB // EXADUNE
    auto test0 = [&] {
      const auto& set = gv.grid().leafIndexSet();
      IndexSetPartitioner<GridLayerType> partitioner(set);
      EXPECT_EQ(set.size(0), partitioner.partitions());
      Dune::SeedListPartitioning<GridType, 0> partitioning(gv, partitioner);
      walker.append(counter);
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
    Walker<GridLayerType> walker(gv);

    size_t filter_count = 0, all_count = 0;
    auto boundaries = [=](const GridLayerType&, const IntersectionType& inter) { return inter.boundary(); };
    auto filter_counter = [&](const IntersectionType&, const EntityType&, const EntityType&) { filter_count++; };
    auto all_counter = [&](const IntersectionType&, const EntityType&, const EntityType&) { all_count++; };

    const auto* on_filter_boundaries = new ApplyOn::FilteredIntersections<GridLayerType>(boundaries);
    const auto* on_all_boundaries = new ApplyOn::BoundaryIntersections<GridLayerType>();
    walker.append(filter_counter, on_filter_boundaries);
    walker.append(all_counter, on_all_boundaries);
    walker.walk();
    EXPECT_EQ(filter_count, all_count);
  }

  void check_partitionsets()
  {
    const auto gv = grid_prv.grid().leafGridView();
    Walker<GridLayerType> walker(gv);

    size_t filter_count = 0, all_count = 0, inner_count = 0, inner_set_count = 0;
    auto filter_counter = [&](const IntersectionType&, const EntityType&, const EntityType&) { filter_count++; };
    auto inner_filter_counter = [&](const IntersectionType&, const EntityType&, const EntityType&) {
      inner_set_count++;
    };
    auto all_counter = [&](const IntersectionType&, const EntityType&, const EntityType&) { all_count++; };
    auto inner_counter = [&](const IntersectionType&, const EntityType& e, const EntityType&) {
      inner_count += e.partitionType() == Dune::PartitionType::InteriorEntity;
    };

    const auto* on_interior_partitionset =
        new ApplyOn::PartitionSetIntersections<GridLayerType, Dune::Partitions::Interior>();
    const auto* on_all_partitionset = new ApplyOn::PartitionSetIntersections<GridLayerType, Dune::Partitions::All>();
    const auto* on_all = new ApplyOn::AllIntersections<GridLayerType>();
    const auto* on_all2 =
        new ApplyOn::AllIntersections<GridLayerType>(); // second one needed because append owns pointer
    walker.append(filter_counter, on_all_partitionset);
    walker.append(inner_filter_counter, on_interior_partitionset);
    walker.append(all_counter, on_all);
    walker.append(inner_counter, on_all2);
    walker.walk();
    EXPECT_EQ(filter_count, all_count);
  }

  void check_partitioning()
  {
    const auto gv = grid_prv.grid().leafGridView();
    Walker<GridLayerType> walker(gv);

    size_t all_count = 0, inner_count = 0;
    auto all_set_counter = [&](const EntityType&) { all_count++; };
    auto inner_set_counter = [&](const EntityType&) { inner_count++; };
    const auto* on_interior_partitionset =
        new ApplyOn::PartitionSetEntities<GridLayerType, Dune::Partitions::Interior>();
    const auto* on_all_partitionset = new ApplyOn::PartitionSetEntities<GridLayerType, Dune::Partitions::All>();
    walker.append(inner_set_counter, on_interior_partitionset);
    walker.append(all_set_counter, on_all_partitionset);
    walker.walk();

    Dune::XT::Grid::RangedPartitioning<GridLayerType, 0, Dune::Interior_Partition> interior_part(gv, 1);
    Dune::XT::Grid::RangedPartitioning<GridLayerType, 0, Dune::All_Partition> all_part(gv, 1);
  }
};

TYPED_TEST_CASE(GridWalkerTest, GridDims);
TYPED_TEST(GridWalkerTest, Misc)
{
  this->check_count();
  this->check_apply_on();
  this->check_partitionsets();
  this->check_partitioning();
}
