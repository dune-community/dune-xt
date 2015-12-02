// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#include "main.hxx"

#if HAVE_DUNE_GRID
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/common/declaration.hh>
#include <dune/alugrid/3d/alugrid.hh>
#endif // HAVE_DUNE_ALUGRID
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

  static void check(const PeriodicGridViewType periodic_grid_view, const GridViewType& grid_view, const size_t variant)
  {
    const bool is_cube = !IS_SIMPLEX;
    DomainType lower_left(0);
    DomainType upper_right(1);
    if (variant >= 3)
      upper_right = DSC::fromString<DomainType>("[2.0 3.0 1.0 4.0]");
    if (variant >= 6)
      lower_left = DSC::fromString<DomainType>("[-0.5 0.5 0.7 -1.3]");

    const int is_nonperiodic        = variant == 0 || variant == 3 || variant == 6 ? 1 : 0;
    const int is_partially_periodic = variant == 1 || variant == 4 || variant == 7 ? 1 : 0;

    // check interface
    const GridType& DSC_UNUSED(test_grid) = periodic_grid_view.grid();
    const IndexSet& DSC_UNUSED(test_indexSet) = periodic_grid_view.indexSet();
    const int codim0_size = periodic_grid_view.size(0);
    EXPECT_EQ(grid_view.size(0), codim0_size);
    if (is_cube)
      EXPECT_EQ(std::pow(int(8), dimDomain), codim0_size);
    if (IS_SIMPLEX)
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
    if (IS_SIMPLEX && dimDomain == 3)
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

  static void checks_for_all_grids()
  {
    // create grid on unit hypercube
    GridProviderType grid_provider             = *(GridProviderType::create());
    const std::shared_ptr<const GridType> grid = grid_provider.grid_ptr();
    const GridViewType grid_view               = grid->leafGridView();

    // create grid on hyperrectangle
    DSC::Configuration grid_config                            = GridProviderType::default_config();
    grid_config["lower_left"]                                 = "[0.0 0.0 0.0 0.0]";
    grid_config["upper_right"]                                = "[2.0 3.0 1.0 4.0]";
    GridProviderType hyperrectangle_grid_provider             = *(GridProviderType::create(grid_config));
    const std::shared_ptr<const GridType> hyperrectangle_grid = hyperrectangle_grid_provider.grid_ptr();
    const GridViewType hyperrectangle_grid_view               = hyperrectangle_grid->leafGridView();

    std::bitset<dimDomain> periodic_directions;
    // create PeriodicGridViews that actually are not periodic
    const PeriodicGridViewType non_periodic_grid_view(grid_view, periodic_directions);
    const PeriodicGridViewType hr_non_periodic_grid_view(hyperrectangle_grid_view, periodic_directions);
    // create PeriodicGridViews that are periodic only in x-direction
    periodic_directions[0] = 1;
    const PeriodicGridViewType partially_periodic_grid_view(grid_view, periodic_directions);
    const PeriodicGridViewType hr_partially_periodic_grid_view(hyperrectangle_grid_view, periodic_directions);
    // create PeriodicGridViews that are periodic in all directions
    periodic_directions.set();
    const PeriodicGridViewType fully_periodic_grid_view(grid_view, periodic_directions);
    const PeriodicGridViewType hr_fully_periodic_grid_view(hyperrectangle_grid_view, periodic_directions);

    // check
    check(non_periodic_grid_view, grid_view, 0);
    check(partially_periodic_grid_view, grid_view, 1);
    check(fully_periodic_grid_view, grid_view, 2);
    check(hr_non_periodic_grid_view, hyperrectangle_grid_view, 3);
    check(hr_partially_periodic_grid_view, hyperrectangle_grid_view, 4);
    check(hr_fully_periodic_grid_view, hyperrectangle_grid_view, 5);
  } // void checks_for_all_grids()

  static void non_trivial_origin_checks()
  {
    // create grid on hyperrectangle with lower_left != 0 (not possible for YaspGrid)
    DSC::Configuration grid_config             = GridProviderType::default_config();
    grid_config["lower_left"]                  = "[-0.5 0.5 0.7 -1.3]";
    grid_config["upper_right"]                 = "[2.0 3.0 1.0 4.0]";
    GridProviderType grid_provider             = *(GridProviderType::create(grid_config));
    const std::shared_ptr<const GridType> grid = grid_provider.grid_ptr();
    const GridViewType grid_view               = grid->leafGridView();

    std::bitset<dimDomain> periodic_directions;
    // create PeriodicGridViewType that actually is not periodic
    const PeriodicGridViewType non_periodic_grid_view(grid_view, periodic_directions);
    // create PeriodicGridViewType that is periodic only in x-direction
    periodic_directions[0] = 1;
    const PeriodicGridViewType partially_periodic_grid_view(grid_view, periodic_directions);
    // create PeriodicGridViewType that is periodic in all directions
    periodic_directions.set();
    const PeriodicGridViewType fully_periodic_grid_view(grid_view, periodic_directions);

    check(non_periodic_grid_view, grid_view, 6);
    check(partially_periodic_grid_view, grid_view, 7);
    check(fully_periodic_grid_view, grid_view, 8);
  } // void non_trivial_origin_checks()
};

TEST_F(PeriodicViewTest, check_all)
{
  this->checks_for_all_grids();
  this->non_trivial_origin_checks();
}

#else // HAVE_DUNE_GRID

TEST(DISABLED_PeriodicViewTestYaspCube, check_yaspcube)
{
}

#endif // HAVE_DUNE_GRID
