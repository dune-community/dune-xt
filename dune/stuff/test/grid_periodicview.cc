// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#include "main.hxx"

#if HAVE_DUNE_GRID
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/grid/periodicview.hh>
#include <dune/stuff/grid/provider/interface.hh>
#include <dune/stuff/grid/provider/cube.hh>

using namespace Dune;
using namespace Stuff;


struct PeriodicViewTest : public testing::Test
{
  typedef TESTGRIDTYPE GridType;
  typedef typename GridType::ctype ctype;
  typedef typename GridType::template Codim<0>::Geometry GeometryType;
  typedef Dune::Stuff::Grid::Providers::template Cube<GridType> GridProviderType;
  typedef typename GridType::LeafGridView GridViewType;
  typedef typename GridViewType::IndexSet IndexSet;
  typedef typename GridViewType::template Codim<0>::Geometry::GlobalCoordinate DomainType;
  typedef typename Dune::Stuff::Grid::template PeriodicGridView<GridViewType> PeriodicGridViewType;
  typedef typename PeriodicGridViewType::template Codim<0>::Entity EntityType;
  typedef typename PeriodicGridViewType::template Codim<0>::Iterator EntityIteratorType;
  typedef typename DSG::internal::template PeriodicIntersectionIterator<GridViewType> PeriodicIntersectionIteratorType;
  typedef typename DSG::internal::template PeriodicIntersection<GridViewType> PeriodicIntersectionType;
  typedef typename PeriodicIntersectionType::EntityPointer EntityPointerType;
  typedef typename GridViewType::CollectiveCommunication CollectiveCommunication;
  static const size_t dimDomain = GridViewType::dimension;

  static constexpr int factorial(int n)
  {
    return n > 0 ? n * factorial(n - 1) : 1;
  }

  static void check()
  {
    // create grid and get gridview
    DSC::Configuration grid_config             = DSC_CONFIG.sub("test_grid_periodicview");
    GridProviderType grid_provider             = *(GridProviderType::create(grid_config));
    const std::shared_ptr<const GridType> grid = grid_provider.grid_ptr();
    const GridViewType grid_view               = grid->leafGridView();

    // check whether grid is periodic
    const bool is_nonperiodic        = grid_config["periodic"] == "no";
    const bool is_partially_periodic = grid_config["periodic"] == "partial";

    // create periodic grid_view
    std::bitset<dimDomain> periodic_directions;
    if (is_partially_periodic)
      periodic_directions[0] = 1;
    else if (!is_nonperiodic)
      periodic_directions.set();
    const PeriodicGridViewType periodic_grid_view(grid_view, periodic_directions);

    const bool is_simplex        = DSC::fromString<bool>(grid_config["is_simplex"]);
    const bool is_cube           = !is_simplex;
    const DomainType lower_left  = DSC::fromString<DomainType>(grid_config["lower_left"]);
    const DomainType upper_right = DSC::fromString<DomainType>(grid_config["upper_right"]);

    // check interface
    const GridType& DSC_UNUSED(test_grid) = periodic_grid_view.grid();
    const IndexSet& DSC_UNUSED(test_indexSet) = periodic_grid_view.indexSet();
    const int codim0_size = periodic_grid_view.size(0);
    EXPECT_EQ(grid_view.size(0), codim0_size);
    if (is_cube)
      EXPECT_EQ(std::pow(int(8), dimDomain), codim0_size);
    if (is_simplex)
      EXPECT_EQ(std::pow(int(8), dimDomain) * factorial(dimDomain), codim0_size);
    EXPECT_EQ(grid_view.size(Dune::GeometryType::cube), periodic_grid_view.size(Dune::GeometryType::cube));
    EXPECT_EQ(grid_view.size(Dune::GeometryType::simplex), periodic_grid_view.size(Dune::GeometryType::simplex));
    EXPECT_EQ(grid_view.overlapSize(0), periodic_grid_view.overlapSize(0));
    EXPECT_EQ(grid_view.overlapSize(1), periodic_grid_view.overlapSize(1));
    EXPECT_EQ(grid_view.ghostSize(0), periodic_grid_view.ghostSize(0));
    EXPECT_EQ(grid_view.ghostSize(1), periodic_grid_view.ghostSize(1));
    const CollectiveCommunication& DSC_UNUSED(test_comm) = periodic_grid_view.comm();

    size_t neighbor_count = 0;
    size_t boundary_count = 0;
    size_t periodic_count = 0;
    // iterate over all entities
    const EntityIteratorType it_end = periodic_grid_view.template end<0>();
    for (EntityIteratorType it = periodic_grid_view.template begin<0>(); it != it_end; ++it) {
      const EntityType& entity = *it;
      EXPECT_TRUE(periodic_grid_view.contains(entity));
      // iterate over all intersections on current entity
      const PeriodicIntersectionIteratorType i_it_end = periodic_grid_view.iend(entity);
      for (PeriodicIntersectionIteratorType i_it = periodic_grid_view.ibegin(entity); i_it != i_it_end; ++i_it) {
        const PeriodicIntersectionType& intersection = *i_it;
        if (intersection.neighbor()) {
          ++neighbor_count;
          const EntityType outside = intersection.outside();
          // find corresponding intersection in outside
          const auto index_in_outside                             = intersection.indexInOutside();
          PeriodicIntersectionIteratorType i_it_outside           = periodic_grid_view.ibegin(outside);
          const PeriodicIntersectionIteratorType i_it_outside_end = periodic_grid_view.iend(outside);
          for (; i_it_outside != i_it_outside_end; ++i_it_outside) {
            const PeriodicIntersectionType* outside_intersection = i_it_outside.operator->();
            if (outside_intersection->indexInInside() == index_in_outside) {
              break;
            }
          }
          const PeriodicIntersectionType& intersection_in_outside = *i_it_outside;
          // check outside_intersection coords
          const auto coords_in_outside   = intersection.geometryInOutside().center();
          const auto coords_in_outside_2 = intersection_in_outside.geometryInInside().center();
          EXPECT_TRUE(Dune::Stuff::Common::FloatCmp::eq(coords_in_outside, coords_in_outside_2));
          // check global intersection coords in periodic case
          const auto global_intersection_coords         = intersection.geometry().center();
          const auto global_outside_intersection_coords = intersection_in_outside.geometry().center();
          size_t coord_difference_count                 = 0;
          size_t differing_coordinate;
          for (size_t ii = 0; ii < dimDomain; ++ii) {
            if (Dune::Stuff::Common::FloatCmp::ne(global_outside_intersection_coords[ii],
                                                  global_intersection_coords[ii])) {
              ++coord_difference_count;
              differing_coordinate = ii;
            }
          }
          if (intersection.boundary() && intersection.neighbor()) {
            EXPECT_TRUE(intersection_in_outside.boundary() && intersection_in_outside.neighbor());
            EXPECT_EQ(size_t(1), coord_difference_count);
            EXPECT_TRUE((Dune::Stuff::Common::FloatCmp::eq(global_outside_intersection_coords[differing_coordinate],
                                                           lower_left[differing_coordinate])
                         && Dune::Stuff::Common::FloatCmp::eq(global_intersection_coords[differing_coordinate],
                                                              upper_right[differing_coordinate]))
                        || (Dune::Stuff::Common::FloatCmp::eq(global_outside_intersection_coords[differing_coordinate],
                                                              upper_right[differing_coordinate])
                            && Dune::Stuff::Common::FloatCmp::eq(global_intersection_coords[differing_coordinate],
                                                                 lower_left[differing_coordinate])));
            ++periodic_count;
            if (is_partially_periodic)
              EXPECT_EQ(size_t(0), differing_coordinate);
          }
        }
        if (intersection.boundary() && !intersection.neighbor()) {
          ++boundary_count;
        }
      }
    }

    // the cube/rectangle grid has 2*dimDomain faces
    const size_t num_faces = 2 * dimDomain;
    // on each face, there are 8**(dimDomain-1) intersections. For a simplex grid in 3 dimensions, there are twice as
    // much.
    size_t num_intersections_on_face = std::pow(8, dimDomain - 1);
    // use dimDomain from config here to avoid "code will never be executed" warning
    assert(dimDomain == DSC::fromString<int>(grid_config["dimDomain"]));
    if (is_simplex && DSC::fromString<int>(grid_config["dimDomain"]) == 3)
      num_intersections_on_face *= 2;
    // In a fully periodic grid, all intersections are periodic. In a partially periodic grid, only the intersections on
    // two
    // faces are periodic. In a nonperiodic grid, no intersections are periodic.
    size_t num_periodic_faces = is_partially_periodic ? 2 : num_faces;
    if (is_nonperiodic)
      num_periodic_faces *= 0;
    const size_t expected_num_periodic_intersections = num_periodic_faces * num_intersections_on_face;
    EXPECT_EQ(expected_num_periodic_intersections, periodic_count);
    // The boundary count should be the number of interfaces on the boundary without the periodic interfaces
    const size_t expected_num_boundary_intersections =
        num_faces * num_intersections_on_face - expected_num_periodic_intersections;
    EXPECT_EQ(expected_num_boundary_intersections, boundary_count);
    // neighbor_count should equal the number of intersections without intersections on non-periodic boundaries
    const size_t num_intersections_per_entity = is_cube ? 2 * dimDomain : dimDomain + 1;
    const size_t num_entities = grid_view.size(0);
    EXPECT_EQ(num_entities * num_intersections_per_entity - expected_num_boundary_intersections, neighbor_count);
  } // void check(...)
};

TEST_F(PeriodicViewTest, check_all)
{
  this->check();
}

#else // HAVE_DUNE_GRID

TEST(DISABLED_PeriodicViewTest, check_all)
{
}

#endif // HAVE_DUNE_GRID
