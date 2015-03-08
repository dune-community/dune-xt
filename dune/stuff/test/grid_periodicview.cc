// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

//# define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
//# define DUNE_STUFF_TEST_MAIN_ENABLE_TIMED_LOGGING 1

#include "main.hxx"

#if HAVE_DUNE_GRID
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/grid/periodicview.hh>
#include <dune/stuff/grid/provider/interface.hh>
#include <dune/stuff/grid/provider/cube.hh>

using namespace Dune;
using namespace Stuff;

#define YASPGRIDS YaspGrid<1>, YaspGrid<2>, YaspGrid<3>

#if HAVE_ALUGRID
#define ALUCUBEGRIDS ALUGrid<2, 2, cube, nonconforming>, ALUGrid<3, 3, cube, nonconforming>

#define ALUSIMPLEXGRIDS ALUGrid<2, 2, simplex, conforming>, ALUGrid<3, 3, simplex, conforming>
#endif // HAVE_ALUGRID


template <class GridImp>
struct PeriodicViewTestCube : public testing::Test
{
  typedef GridImp GridType;
  typedef typename GridType::ctype ctype;
  typedef typename GridType::template Codim<0>::Geometry GeometryType;
  typedef Dune::Stuff::Grid::Providers::template Cube<GridType> GridProviderType;
  typedef typename GridType::LeafGridView GridViewType;
  typedef typename GridViewType::IndexSet IndexSet;
  typedef typename Dune::Stuff::Grid::template PeriodicGridView<GridViewType> PeriodicGridViewType;
  typedef typename PeriodicGridViewType::template Codim<0>::Entity EntityType;
  typedef typename PeriodicGridViewType::template Codim<0>::Iterator EntityIteratorType;
  typedef
      typename Dune::Stuff::Grid::template PeriodicIntersectionIterator<GridViewType> PeriodicIntersectionIteratorType;
  typedef typename Dune::Stuff::Grid::template PeriodicIntersection<GridViewType> PeriodicIntersectionType;
  typedef typename PeriodicIntersectionType::EntityPointer EntityPointerType;
  typedef typename GridViewType::CollectiveCommunication CollectiveCommunication;
  static const size_t dimDomain = GridViewType::dimension;

  void check(const bool is_simplex)
  {
    const bool is_cube                         = !is_simplex;
    GridProviderType grid_provider             = *(GridProviderType::create());
    const std::shared_ptr<const GridType> grid = grid_provider.grid_ptr();
    const GridViewType grid_view               = grid->leafGridView();
    const GridViewType& grid_view_ref          = grid_view;

    std::bitset<dimDomain> periodic_directions;
    // create PeriodicGridviewType that actually is not periodic
    const PeriodicGridViewType non_periodic_grid_view(grid_view_ref, periodic_directions);
    // create PeriodicGridviewType that is periodic only in x-direction
    periodic_directions[0] = 1;
    const PeriodicGridViewType partially_periodic_grid_view(grid_view_ref, periodic_directions);
    // create PeriodicGridviewType that is periodic only in all directions
    periodic_directions.set();
    const PeriodicGridViewType fully_periodic_grid_view(grid_view_ref, periodic_directions);

    std::vector<const PeriodicGridViewType*> periodic_grid_view_vector(3);
    periodic_grid_view_vector[0] = &non_periodic_grid_view;
    periodic_grid_view_vector[1] = &partially_periodic_grid_view;
    periodic_grid_view_vector[2] = &fully_periodic_grid_view;

    for (size_t variant = 0; variant <= 2; ++variant) {
      const PeriodicGridViewType periodic_grid_view = *(periodic_grid_view_vector[variant]);

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
      EXPECT_EQ(periodic_grid_view.size(Dune::cube), grid_view.size(Dune::cube));
      EXPECT_EQ(periodic_grid_view.size(Dune::simplex), grid_view.size(Dune::simplex));
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
            const EntityPointerType outside = intersection.outside();
            // find corresponding intersection in outside
            const auto index_in_outside                             = intersection.indexInOutside();
            PeriodicIntersectionIteratorType i_it_outside           = periodic_grid_view.ibegin(*outside);
            const PeriodicIntersectionIteratorType i_it_outside_end = periodic_grid_view.iend(*outside);
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
            if (intersection.periodic()) {
              EXPECT_TRUE(intersection_in_outside.periodic());
              EXPECT_FALSE(intersection_in_outside.boundary());
              EXPECT_EQ(size_t(1), coord_difference_count);
              EXPECT_TRUE(
                  (Dune::Stuff::Common::FloatCmp::eq(global_outside_intersection_coords[differing_coordinate], ctype(1))
                   && Dune::Stuff::Common::FloatCmp::eq(global_intersection_coords[differing_coordinate], ctype(0)))
                  || (Dune::Stuff::Common::FloatCmp::eq(global_outside_intersection_coords[differing_coordinate],
                                                        ctype(0))
                      && Dune::Stuff::Common::FloatCmp::eq(global_intersection_coords[differing_coordinate],
                                                           ctype(1))));
              ++periodic_count;
              if (variant == 1)
                EXPECT_EQ(size_t(0), differing_coordinate);
            }
          }
          if (intersection.boundary()) {
            ++boundary_count;
            EXPECT_FALSE(intersection.periodic());
            EXPECT_FALSE(intersection.neighbor());
          }
        }
      }

      if (dimDomain == 1) {
        if (variant == 0) {
          EXPECT_EQ(size_t(2 * 7), neighbor_count);
          EXPECT_EQ(size_t(2), boundary_count);
          EXPECT_EQ(size_t(0), periodic_count);
        }
        if (variant == 1) {
          EXPECT_EQ(size_t(2 * 7 + 2), neighbor_count);
          EXPECT_EQ(size_t(0), boundary_count);
          EXPECT_EQ(size_t(2), periodic_count);
        }
        if (variant == 2) {
          EXPECT_EQ(size_t(2 * 7 + 2), neighbor_count);
          EXPECT_EQ(size_t(0), boundary_count);
          EXPECT_EQ(size_t(2), periodic_count);
        }
      }
      if (dimDomain == 2 && is_cube) {
        if (variant == 0) {
          EXPECT_EQ(size_t(2 * 2 * 8 * 7), neighbor_count);
          EXPECT_EQ(size_t(2 * 8 * 2), boundary_count);
          EXPECT_EQ(size_t(0), periodic_count);
        }
        if (variant == 1) {
          EXPECT_EQ(size_t(2 * 2 * 8 * 7 + 2 * 8), neighbor_count);
          EXPECT_EQ(size_t(2 * 8 * 2 - 2 * 8), boundary_count);
          EXPECT_EQ(size_t(2 * 8), periodic_count);
        }
        if (variant == 2) {
          EXPECT_EQ(size_t(2 * 2 * 8 * 7 + 2 * 8 * 2), neighbor_count);
          EXPECT_EQ(size_t(0), boundary_count);
          EXPECT_EQ(size_t(2 * 8 * 2), periodic_count);
        }
      }
      if (dimDomain == 2 && is_simplex) {
        if (variant == 0) {
          EXPECT_EQ(size_t(2 * 2 * 8 * 7 + 2 * 64), neighbor_count);
          EXPECT_EQ(size_t(2 * 8 * 2), boundary_count);
          EXPECT_EQ(size_t(0), periodic_count);
        }
        if (variant == 1) {
          EXPECT_EQ(size_t(2 * 2 * 8 * 7 + 2 * 8 + 2 * 64), neighbor_count);
          EXPECT_EQ(size_t(2 * 8 * 2 - 2 * 8), boundary_count);
          EXPECT_EQ(size_t(2 * 8), periodic_count);
        }
        if (variant == 2) {
          EXPECT_EQ(size_t(2 * 2 * 8 * 7 + 2 * 8 * 2 + 2 * 64), neighbor_count);
          EXPECT_EQ(size_t(0), boundary_count);
          EXPECT_EQ(size_t(2 * 8 * 2), periodic_count);
        }
      }
      if (dimDomain == 3 && is_cube) {
        if (variant == 0) {
          EXPECT_EQ(size_t(2 * 3 * 7 * 64), neighbor_count);
          EXPECT_EQ(size_t(6 * 64), boundary_count);
          EXPECT_EQ(size_t(0), periodic_count);
        }
        if (variant == 1) {
          EXPECT_EQ(size_t(2 * 3 * 7 * 64 + 2 * 64), neighbor_count);
          EXPECT_EQ(size_t(6 * 64 - 2 * 64), boundary_count);
          EXPECT_EQ(size_t(2 * 64), periodic_count);
        }
        if (variant == 2) {
          EXPECT_EQ(size_t(2 * 3 * 7 * 64 + 6 * 64), neighbor_count);
          EXPECT_EQ(size_t(0), boundary_count);
          EXPECT_EQ(size_t(6 * 64), periodic_count);
        }
      }
      if (dimDomain == 3 && is_simplex) {
        if (variant == 0) {
          EXPECT_EQ(size_t(2 * 2 * 3 * 7 * 64 + 2 * 6 * 8 * 8 * 8), neighbor_count);
          EXPECT_EQ(size_t(2 * 6 * 64), boundary_count);
          EXPECT_EQ(size_t(0), periodic_count);
        }
        if (variant == 1) {
          EXPECT_EQ(size_t(2 * 2 * 3 * 7 * 64 + 2 * 2 * 64 + 2 * 6 * 8 * 8 * 8), neighbor_count);
          EXPECT_EQ(size_t(2 * 6 * 64 - 2 * 2 * 64), boundary_count);
          EXPECT_EQ(size_t(2 * 2 * 64), periodic_count);
        }
        if (variant == 2) {
          EXPECT_EQ(size_t(2 * 2 * 3 * 7 * 64 + 2 * 6 * 64 + 2 * 6 * 8 * 8 * 8), neighbor_count);
          EXPECT_EQ(size_t(0), boundary_count);
          EXPECT_EQ(size_t(2 * 6 * 64), periodic_count);
        }
      }
    }
  }
};

template <class GridImp>
struct PeriodicViewTestSimplex : public PeriodicViewTestCube<GridImp>
{
};

typedef testing::Types<YASPGRIDS
#if HAVE_ALUGRID
                       ,
                       ALUCUBEGRIDS
#endif
                       > CubeGridTypes;

typedef testing::Types<
#if HAVE_ALUGRID
    ALUSIMPLEXGRIDS
#endif
    > SimplexGridTypes;


TYPED_TEST_CASE(PeriodicViewTestCube, CubeGridTypes);
TYPED_TEST(PeriodicViewTestCube, checkcube)
{
  this->check(false);
}

TYPED_TEST_CASE(PeriodicViewTestSimplex, SimplexGridTypes);
TYPED_TEST(PeriodicViewTestSimplex, checksimplex)
{
  this->check(true);
}


#else // HAVE_DUNE_GRID

TEST(DISABLED_CubeGridProvider, is_default_creatable)
{
}
TEST(DISABLED_CubeGridProvider, fulfills_const_interface)
{
}
TEST(DISABLED_CubeGridProvider, is_visualizable)
{
}

#endif // HAVE_DUNE_GRID
