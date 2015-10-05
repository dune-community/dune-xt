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

#define YASPGRIDS                                                                                                      \
  YaspGrid<1, EquidistantOffsetCoordinates<double, 1>>, YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>,          \
      YaspGrid<3, EquidistantOffsetCoordinates<double, 3>>

#if HAVE_ALUGRID
#define ALUCUBEGRIDS ALUGrid<2, 2, cube, nonconforming>, ALUGrid<3, 3, cube, nonconforming>

#define ALUSIMPLEXGRIDS ALUGrid<2, 2, simplex, conforming>, ALUGrid<3, 3, simplex, conforming>
#endif // HAVE_ALUGRID


template <class GridImp>
struct PeriodicViewTestYaspCube : public testing::Test
{
  typedef GridImp GridType;
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

  void check(const PeriodicGridViewType periodic_grid_view, const GridViewType& grid_view, const bool is_simplex,
             const size_t variant)
  {
    const bool is_cube = !is_simplex;
    DomainType lower_left(0);
    DomainType upper_right(1);
    if (variant >= 3)
      upper_right = DSC::fromString<DomainType>("[2.0 3.0 1.0 4.0]");
    if (variant >= 6)
      lower_left = DSC::fromString<DomainType>("[-0.5 0.5 0.7 -1.3]");

    // check interface
    const GridType& DUNE_UNUSED(test_grid) = periodic_grid_view.grid();
    const IndexSet& DUNE_UNUSED(test_indexSet) = periodic_grid_view.indexSet();
    const int codim0_size = periodic_grid_view.size(0);
    EXPECT_EQ(codim0_size, grid_view.size(0));
    if (dimDomain == 1)
      EXPECT_EQ(codim0_size, int(8));
    if (dimDomain == 2 && is_cube)
      EXPECT_EQ(codim0_size, int(8 * 8));
    if (dimDomain == 2 && is_simplex)
      EXPECT_EQ(codim0_size, int(2 * 8 * 8));
    if (dimDomain == 3 && is_cube)
      EXPECT_EQ(codim0_size, int(8 * 8 * 8));
    if (dimDomain == 3 && is_simplex)
      EXPECT_EQ(codim0_size, int(6 * 8 * 8 * 8));
    EXPECT_EQ(periodic_grid_view.size(Dune::GeometryType::cube), grid_view.size(Dune::GeometryType::cube));
    EXPECT_EQ(periodic_grid_view.size(Dune::GeometryType::simplex), grid_view.size(Dune::GeometryType::simplex));
    EXPECT_EQ(periodic_grid_view.overlapSize(0), grid_view.overlapSize(0));
    EXPECT_EQ(periodic_grid_view.overlapSize(1), grid_view.overlapSize(1));
    EXPECT_EQ(periodic_grid_view.ghostSize(0), grid_view.ghostSize(0));
    EXPECT_EQ(periodic_grid_view.ghostSize(1), grid_view.ghostSize(1));
    const CollectiveCommunication& DUNE_UNUSED(test_comm) = periodic_grid_view.comm();

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
            if (variant == 1 || variant == 4 || variant == 7)
              EXPECT_EQ(size_t(0), differing_coordinate);
          }
        }
        if (intersection.boundary() && !intersection.neighbor()) {
          ++boundary_count;
        }
      }
    }

    if (dimDomain == 1) {
      if (variant == 0 || variant == 3 || variant == 6) {
        EXPECT_EQ(size_t(2 * 7), neighbor_count);
        EXPECT_EQ(size_t(2), boundary_count);
        EXPECT_EQ(size_t(0), periodic_count);
      }
      if (variant == 1 || variant == 4 || variant == 7) {
        EXPECT_EQ(size_t(2 * 7 + 2), neighbor_count);
        EXPECT_EQ(size_t(0), boundary_count);
        EXPECT_EQ(size_t(2), periodic_count);
      }
      if (variant == 2 || variant == 5 || variant == 8) {
        EXPECT_EQ(size_t(2 * 7 + 2), neighbor_count);
        EXPECT_EQ(size_t(0), boundary_count);
        EXPECT_EQ(size_t(2), periodic_count);
      }
    }
    if (dimDomain == 2 && is_cube) {
      if (variant == 0 || variant == 3 || variant == 6) {
        EXPECT_EQ(size_t(2 * 2 * 8 * 7), neighbor_count);
        EXPECT_EQ(size_t(2 * 8 * 2), boundary_count);
        EXPECT_EQ(size_t(0), periodic_count);
      }
      if (variant == 1 || variant == 4 || variant == 7) {
        EXPECT_EQ(size_t(2 * 2 * 8 * 7 + 2 * 8), neighbor_count);
        EXPECT_EQ(size_t(2 * 8 * 2 - 2 * 8), boundary_count);
        EXPECT_EQ(size_t(2 * 8), periodic_count);
      }
      if (variant == 2 || variant == 5 || variant == 8) {
        EXPECT_EQ(size_t(2 * 2 * 8 * 7 + 2 * 8 * 2), neighbor_count);
        EXPECT_EQ(size_t(0), boundary_count);
        EXPECT_EQ(size_t(2 * 8 * 2), periodic_count);
      }
    }
    if (dimDomain == 2 && is_simplex) {
      if (variant == 0 || variant == 3 || variant == 6) {
        EXPECT_EQ(size_t(2 * 2 * 8 * 7 + 2 * 64), neighbor_count);
        EXPECT_EQ(size_t(2 * 8 * 2), boundary_count);
        EXPECT_EQ(size_t(0), periodic_count);
      }
      if (variant == 1 || variant == 4 || variant == 7) {
        EXPECT_EQ(size_t(2 * 2 * 8 * 7 + 2 * 8 + 2 * 64), neighbor_count);
        EXPECT_EQ(size_t(2 * 8 * 2 - 2 * 8), boundary_count);
        EXPECT_EQ(size_t(2 * 8), periodic_count);
      }
      if (variant == 2 || variant == 5 || variant == 8) {
        EXPECT_EQ(size_t(2 * 2 * 8 * 7 + 2 * 8 * 2 + 2 * 64), neighbor_count);
        EXPECT_EQ(size_t(0), boundary_count);
        EXPECT_EQ(size_t(2 * 8 * 2), periodic_count);
      }
    }
    if (dimDomain == 3 && is_cube) {
      if (variant == 0 || variant == 3 || variant == 6) {
        EXPECT_EQ(size_t(2 * 3 * 7 * 64), neighbor_count);
        EXPECT_EQ(size_t(6 * 64), boundary_count);
        EXPECT_EQ(size_t(0), periodic_count);
      }
      if (variant == 1 || variant == 4 || variant == 7) {
        EXPECT_EQ(size_t(2 * 3 * 7 * 64 + 2 * 64), neighbor_count);
        EXPECT_EQ(size_t(6 * 64 - 2 * 64), boundary_count);
        EXPECT_EQ(size_t(2 * 64), periodic_count);
      }
      if (variant == 2 || variant == 5 || variant == 8) {
        EXPECT_EQ(size_t(2 * 3 * 7 * 64 + 6 * 64), neighbor_count);
        EXPECT_EQ(size_t(0), boundary_count);
        EXPECT_EQ(size_t(6 * 64), periodic_count);
      }
    }
    if (dimDomain == 3 && is_simplex) {
      if (variant == 0 || variant == 3 || variant == 6) {
        EXPECT_EQ(size_t(2 * 2 * 3 * 7 * 64 + 2 * 6 * 8 * 8 * 8), neighbor_count);
        EXPECT_EQ(size_t(2 * 6 * 64), boundary_count);
        EXPECT_EQ(size_t(0), periodic_count);
      }
      if (variant == 1 || variant == 4 || variant == 7) {
        EXPECT_EQ(size_t(2 * 2 * 3 * 7 * 64 + 2 * 2 * 64 + 2 * 6 * 8 * 8 * 8), neighbor_count);
        EXPECT_EQ(size_t(2 * 6 * 64 - 2 * 2 * 64), boundary_count);
        EXPECT_EQ(size_t(2 * 2 * 64), periodic_count);
      }
      if (variant == 2 || variant == 5 || variant == 8) {
        EXPECT_EQ(size_t(2 * 2 * 3 * 7 * 64 + 2 * 6 * 64 + 2 * 6 * 8 * 8 * 8), neighbor_count);
        EXPECT_EQ(size_t(0), boundary_count);
        EXPECT_EQ(size_t(2 * 6 * 64), periodic_count);
      }
    }
  } // void check(...)

  void checks_for_all_grids(const bool is_simplex)
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
    this->check(non_periodic_grid_view, grid_view, is_simplex, 0);
    this->check(partially_periodic_grid_view, grid_view, is_simplex, 1);
    this->check(fully_periodic_grid_view, grid_view, is_simplex, 2);
    this->check(hr_non_periodic_grid_view, hyperrectangle_grid_view, is_simplex, 3);
    this->check(hr_partially_periodic_grid_view, hyperrectangle_grid_view, is_simplex, 4);
    this->check(hr_fully_periodic_grid_view, hyperrectangle_grid_view, is_simplex, 5);
  } // void checks_for_all_grids(...)

  void non_trivial_origin_checks(const bool is_simplex)
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

    this->check(non_periodic_grid_view, grid_view, is_simplex, 6);
    this->check(partially_periodic_grid_view, grid_view, is_simplex, 7);
    this->check(fully_periodic_grid_view, grid_view, is_simplex, 8);
  } // void additional_checks_for_alu(...)
}; // ... struct PeriodicViewTestYaspCube ...

template <class GridImp>
struct PeriodicViewTestALUCube : public PeriodicViewTestYaspCube<GridImp>
{
  typedef PeriodicViewTestYaspCube<GridImp> BaseType;
  typedef typename BaseType::GridProviderType GridProviderType;
  typedef typename BaseType::GridType GridType;
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::PeriodicGridViewType PeriodicGridViewType;
  static const size_t dimDomain = BaseType::dimDomain;


}; // ... struct PeriodicViewTestALUCube ...

template <class GridImp>
struct PeriodicViewTestALUSimplex : public PeriodicViewTestALUCube<GridImp>
{
};

typedef testing::Types<YASPGRIDS> YaspCubeGridTypes;

TYPED_TEST_CASE(PeriodicViewTestYaspCube, YaspCubeGridTypes);
TYPED_TEST(PeriodicViewTestYaspCube, check_yaspcube)
{
  this->checks_for_all_grids(false);
  this->non_trivial_origin_checks(false);
}

#if HAVE_ALUGRID

typedef testing::Types<ALUCUBEGRIDS> ALUCubeGridTypes;

typedef testing::Types<ALUSIMPLEXGRIDS> ALUSimplexGridTypes;

TYPED_TEST_CASE(PeriodicViewTestALUCube, ALUCubeGridTypes);
TYPED_TEST(PeriodicViewTestALUCube, check_alucube)
{
  this->checks_for_all_grids(false);
  this->non_trivial_origin_checks(false);
}

TYPED_TEST_CASE(PeriodicViewTestALUSimplex, ALUSimplexGridTypes);
TYPED_TEST(PeriodicViewTestALUSimplex, check_alusimplex)
{
  this->checks_for_all_grids(true);
  this->non_trivial_origin_checks(true);
}

#else // HAVE_ALUGRID

TEST(DISABLED_PeriodicViewTestALUCube, check_alucube)
{
}
TEST(DISABLED_PeriodicViewTestALUSimplex, check_alusimplex)
{
}

#endif // HAVE_ALUGRID
#else // HAVE_DUNE_GRID

TEST(DISABLED_PeriodicViewTestYaspCube, check_yaspcube)
{
}

#endif // HAVE_DUNE_GRID
