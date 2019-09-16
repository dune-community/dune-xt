// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2019)
//   René Fritze     (2014 - 2019)
//   Tobias Leibner  (2015 - 2016, 2018)

#include <dune/xt/test/main.hxx>

#if DUNE_VERSION_GTE(DUNE_COMMON, 3, 9) && HAVE_TBB // EXADUNE
#  include <dune/grid/utility/partitioning/seedlist.hh>
#endif

#include <dune/xt/common/logstreams.hh>
#include <dune/xt/common/parallel/partitioner.hh>

#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/functors/boundary-detector.hh>
#include <dune/xt/grid/parallel/partitioning/ranged.hh>
#include <dune/xt/grid/walker.hh>


using namespace Dune::XT::Common;
using namespace Dune::XT::Grid;
using namespace std;

typedef testing::Types<Int<1>, Int<2>, Int<3>> GridDims;

template <class T>
struct GridWalkerTest : public ::testing::Test
{
  static const size_t griddim = T::value;
  static const size_t level = 128;
  typedef Dune::YaspGrid<griddim, Dune::EquidistantOffsetCoordinates<double, griddim>> GridType;
  typedef typename GridType::LeafGridView GridLayerType;
  using EntityType = extract_entity_t<GridLayerType>;
  using IntersectionType = extract_intersection_t<GridLayerType>;
  const GridProvider<GridType> grid_prv;
  GridWalkerTest()
    : grid_prv(make_cube_grid<GridType>(0.f, 1.f, level))
  {}

  void check_count()
  {
    const auto gv = grid_prv.grid().leafGridView();
    Walker<GridLayerType> walker(gv);
    const auto num_elements = gv.size(0);
    atomic<size_t> count(0);
    atomic<size_t> intersection_count(0);
    auto counter = GenericElementFunctor<GridLayerType>([] {}, [&count](const EntityType&) { count++; }, [] {});
    auto intersection_counter = GenericIntersectionFunctor<GridLayerType>(
        [] {}, [&](const IntersectionType&, const EntityType&, const EntityType&) { intersection_count++; }, [] {});
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

    list<function<void()>> element_tests({test1, test2, test3});
    list<function<void()>> intersection_tests({test4, test5});
#if DUNE_VERSION_GTE(DUNE_COMMON, 3, 9) && HAVE_TBB // EXADUNE
    // exadune guard for SeedListPartitioning
    auto test0 = [&] {
      const auto& set = gv.grid().leafIndexSet();
      IndexSetPartitioner<GridLayerType> partitioner(set);
      EXPECT_EQ(set.size(0), partitioner.partitions());
      Dune::SeedListPartitioning<GridType, 0> partitioning(gv, partitioner);
      walker.append(counter);
      walker.walk(partitioning);
    };
    tests.push_back(test0);
#endif // DUNE_VERSION_GTE(DUNE_COMMON, 3, 9) && HAVE_TBB

    for (const auto& test : element_tests) {
      count = 0;
      test();
      EXPECT_EQ(num_elements, count);
    }
    const auto faces_per_element = 2 * griddim; // only for cube grids
    for (const auto& test : intersection_tests) {
      intersection_count = 0;
      test();
      EXPECT_EQ(num_elements * faces_per_element, intersection_count);
    }
  }

  void check_apply_on()
  {
    const auto gv = grid_prv.grid().leafGridView();
    Walker<GridLayerType> walker(gv);

    size_t filter_count = 0, all_count = 0;
    auto boundaries = [=](const GridLayerType&, const IntersectionType& inter) { return inter.boundary(); };
    auto filter_counter = GenericIntersectionFunctor<GridLayerType>(
        [] {}, [&](const IntersectionType&, const EntityType&, const EntityType&) { filter_count++; }, [] {});
    auto all_counter = GenericIntersectionFunctor<GridLayerType>(
        [] {}, [&](const IntersectionType&, const EntityType&, const EntityType&) { all_count++; }, [] {});

    ApplyOn::GenericFilteredIntersections<GridLayerType> on_filter_boundaries{boundaries};
    ApplyOn::BoundaryIntersections<GridLayerType> on_all_boundaries{};
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
    auto filter_counter = GenericElementFunctor<GridLayerType>([] {}, [&](const auto&) { filter_count++; }, [] {});
    auto inner_filter_counter =
        GenericElementFunctor<GridLayerType>([] {}, [&](const auto&) { inner_set_count++; }, [] {});
    auto all_counter = GenericElementFunctor<GridLayerType>([] {}, [&](const auto&) { all_count++; }, [] {});
    auto inner_counter = GenericElementFunctor<GridLayerType>(
        [] {}, [&](const auto& e) { inner_count += e.partitionType() == Dune::PartitionType::InteriorEntity; }, [] {});

    ApplyOn::PartitionSetElements<GridLayerType, Dune::Partitions::Interior> on_interior_partitionset{};
    ApplyOn::PartitionSetElements<GridLayerType, Dune::Partitions::All> on_all_partitionset{};
    ApplyOn::AllElements<GridLayerType> on_all{};
    walker.append(filter_counter, on_all_partitionset);
    walker.append(inner_filter_counter, on_interior_partitionset);
    walker.append(all_counter, on_all);
    walker.append(inner_counter, on_all);
    walker.walk();
    EXPECT_EQ(filter_count, all_count);
    EXPECT_EQ(inner_set_count, inner_count);
  }

  void check_partitioning()
  {
    const auto gv = grid_prv.grid().leafGridView();
    Walker<GridLayerType> walker(gv);

    size_t all_count = 0, inner_count = 0;
    auto all_set_counter = GenericElementFunctor<GridLayerType>([] {}, [&](const auto&) { all_count++; }, [] {});
    auto inner_set_counter = GenericElementFunctor<GridLayerType>([] {}, [&](const auto&) { inner_count++; }, [] {});
    ApplyOn::PartitionSetElements<GridLayerType, Dune::Partitions::Interior> on_interior_partitionset{};
    ApplyOn::PartitionSetElements<GridLayerType, Dune::Partitions::All> on_all_partitionset{};
    walker.append(inner_set_counter, on_interior_partitionset);
    walker.append(all_set_counter, on_all_partitionset);
    walker.walk();

    Dune::XT::Grid::RangedPartitioning<GridLayerType, 0, Dune::Interior_Partition> interior_part(gv, 1);
    Dune::XT::Grid::RangedPartitioning<GridLayerType, 0, Dune::All_Partition> all_part(gv, 1);

    auto filter_inner = inner_count;

    ApplyOn::AllElements<GridLayerType> on_all_elements{};
    walker.append(all_set_counter, on_all_elements);
    walker.append(inner_set_counter, on_all_elements).walk(interior_part);
    EXPECT_EQ(inner_count, 2 * filter_inner);
    if (gv.grid().comm().size() > 1)
      EXPECT_LT(inner_count, all_count);
    else
      EXPECT_EQ(inner_count, all_count);
    walker.append(all_set_counter, on_all_elements);
    walker.append(inner_set_counter, on_all_elements).walk(all_part);
    EXPECT_EQ(inner_count, 3 * filter_inner);
    if (gv.grid().comm().size() > 1)
      EXPECT_LT(inner_count, all_count);
    else
      EXPECT_EQ(inner_count, all_count);
  }

  void check_boundaries()
  {
    const auto gv = grid_prv.grid().leafGridView();
    Walker<GridLayerType> walker(gv);

    std::atomic<size_t> filter_count{0};
    auto filter_counter = GenericIntersectionFunctor<GridLayerType>(
        [] {}, [&](const IntersectionType&, const EntityType&, const EntityType&) { filter_count++; }, [] {});
    const auto info = make_alldirichlet_boundaryinfo(gv);
    BoundaryDetectorFunctor<GridLayerType> detector(*info, new DirichletBoundary());

    ApplyOn::BoundaryIntersections<GridLayerType> on_all_boundaries;
    walker.append(filter_counter, ApplyOn::BoundaryIntersections<GridLayerType>());
    walker.walk(false);
    walker.clear();
    walker.append(detector, on_all_boundaries);
    Walker<GridLayerType> walker_copy(walker);
    walker.walk(true);
    EXPECT_EQ(filter_count, detector.result());
    walker_copy.walk(true);
    EXPECT_EQ(2 * filter_count, detector.result());
  }

  void check_walker_to_walker()
  {
    const auto gv = grid_prv.grid().leafGridView();

    Walker<GridLayerType> inner_walker(gv);
    // This functor is restricted to some elements of the grid view.
    size_t num_elements_applied = 0;
    auto functor =
        GenericElementFunctor<GridLayerType>([] {}, [&](const auto& /*element*/) { ++num_elements_applied; }, [] {});
    inner_walker.append(functor,
                        ApplyOn::GenericFilteredElements<GridLayerType>([](const auto& grid_view, const auto& element) {
                          return grid_view.indexSet().index(element) < grid_view.indexSet().size(0) / 2;
                        }));
    inner_walker.walk(false, /*clear_functors=*/false); // We want to keep the functor.
    ASSERT_LT(num_elements_applied, gv.indexSet().size(0));
    auto num_half_elements = num_elements_applied;
    num_elements_applied = 0;

    // When we add this walker to another walker without a restricting filter ...
    Walker<GridLayerType> walker(gv);
    walker.append(inner_walker);
    walker.walk();
    // ... we also expect the functor of inner_walker to be applied only to some elements of the grid view, according to
    // its original filter.
    EXPECT_EQ(num_elements_applied, num_half_elements);
  }
};

TYPED_TEST_CASE(GridWalkerTest, GridDims);
TYPED_TEST(GridWalkerTest, count)
{
  this->check_count();
}
TYPED_TEST(GridWalkerTest, apply_on)
{
  this->check_apply_on();
  this->check_partitionsets();
}
TYPED_TEST(GridWalkerTest, boundaries)
{
  this->check_boundaries();
  this->check_partitioning();
}
TYPED_TEST(GridWalkerTest, walker_to_walker)
{
  this->check_walker_to_walker();
}
