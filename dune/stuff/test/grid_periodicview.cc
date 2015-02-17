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

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/grid/periodicview.hh>
#include <dune/stuff/grid/provider/interface.hh>
#include <dune/stuff/grid/provider/cube.hh>

using namespace Dune;
using namespace Stuff;

#define YASPGRIDS YaspGrid<1>, YaspGrid<2>, YaspGrid<3>, YaspGrid<4>

#if HAVE_ALUGRID
#define ALUGRIDS                                                                                                       \
  ALUGrid<2, 2, cube, nonconforming>, ALUGrid<3, 3, cube, nonconforming>, ALUGrid<2, 2, simplex, conforming>,          \
      ALUGrid<3, 3, simplex, conforming>
#endif // HAVE_ALUGRID


template <class GridImp>
struct PeriodicViewTestBase : public testing::Test
{
  typedef typename GridImp GridType;
  typedef Grid::Providers::Cube<GridType> GridProviderType;
  typedef typename GridType::LeafGridView GridViewType;
  typedef typename Grid::PeriodicGridView<GridViewType> PeriodicGridViewType;
  typedef typename PeriodicGridViewType::template Codim<0>::Entity EntityType;
  typedef typename PeriodicGridViewType::template Codim<0>::Iterator EntityIteratorType;
  typedef typename Grid::PeriodicIntersectionIterator<GridViewType> PeriodicIntersectionIteratorType;
  typedef typename Grid::PeriodicIntersection<GridViewType> PeriodicIntersectionType;
  static const size_t dimDomain = GridViewType::dimension;

  void check()
  {

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
    const PeriodicGridViewType periodic_grid_view(grid_view_ref, periodic_directions);

    size_t neighbor_count           = 0;
    size_t boundary_count           = 0;
    const EntityIteratorType it_end = periodic_grid_view.end<0>();
    for (EntityIteratorType it = periodic_grid_view.begin<0>(); it != it_end; ++it) {
      const EntityType& entity                        = *it;
      const PeriodicIntersectionIteratorType i_it_end = periodic_grid_view.iend(entity);
      for (PeriodicIntersectionIteratorType i_it = periodic_grid_view.ibegin(entity); i_it != i_it_end; ++i_it) {
        const PeriodicIntersectionType& intersection = *i_it;
        if (intersection.neighbor())
          ++neighbor_count;
        const auto outside_ptr         = intersection.outside();
        const auto coords_in_outside   = intersection.geometryInOutside().center();
        const auto coords_in_outside_2 = outside_ptr->geometryInInside().center();
        EXPECT_TRUE(Dune::Stuff::Common::FloatCmp::eq(coords_in_outside, coords_in_outside_2));
        const auto index_in_outside           = intersection_indexInOutside();
        const auto global_intersection_coords = intersection.geometry().center();
        const auto global_outside_coords      = outside_ptr->geometry().center();
        const size_t coord_difference         = 0;
        const size_t differing_coordinate;
        for (size_t ii = 0; ii < dimDomain; ++ii) {
          if (global_outside_coords[ii] != global_intersection_coords) {
            ++coord_difference;
            differing_coordinate = ii;
          }
        }
        if (intersection.is_periodic()) {
          EXPECT_EQ(size_t(1), coord_difference);
          EXPECT_TRUE((Dune::Stuff::Common::FloatCmp::eq(global_outside_coords[differing_coordinate], 1)
                       && Dune::Stuff::Common::FloatCmp::eq(global_intersection_coords[differing_coordinate], 0))
                      || (Dune::Stuff::Common::FloatCmp::eq(global_outside_coords[differing_coordinate], 0)
                          && Dune::Stuff::Common::FloatCmp::eq(global_intersection_coords[differing_coordinate], 1)));
        }

        if (intersection.boundary())
          ++boundary_count;
      }
    }


    const GridType& DUNE_UNUSED(test_grid) = periodic_grid_view.grid();
    const IndexSet& DUNE_UNUSED(test_indexSet) = periodic_grid_view.indexSet();
    size_t num_entities = 1;
    for (size_t i = 0; i < dimDomain; ++i) {
      num_entities *= 8;
    }
    EXPECT_EQ(num_entities, periodic_grid_view.size(0));
    EXPECT_EQ(num_entities, periodic_grid_view.size(1));

      int 	size (const GeometryType &type) const
    obtain number of entities with a given geometry type

          template<class EntityType >
  bool 	contains (const EntityType &e) const
    Return true if the given entity is contained in this grid view.
  template<int cd>
  Codim< cd >::Iterator 	begin () const
    obtain begin iterator for this view
  template<int cd>
  Codim< cd >::Iterator 	end () const
    obtain end iterator for this view
  template<int cd, PartitionIteratorType pitype>
  Codim< cd >::template
  Partition< pitype >::Iterator 	begin () const
    obtain begin iterator for this view
  template<int cd, PartitionIteratorType pitype>
  Codim< cd >::template
  Partition< pitype >::Iterator 	end () const
    obtain end iterator for this view
  IntersectionIterator 	ibegin (const typename Codim< 0 >::Entity &entity) const
    obtain begin intersection iterator with respect to this view
  IntersectionIterator 	iend (const typename Codim< 0 >::Entity &entity) const
    obtain end intersection iterator with respect to this view
  const CollectiveCommunication & 	comm () const
    obtain collective communication object
  int 	overlapSize (int codim) const
    Return size of the overlap region for a given codim on the grid view.
  int 	ghostSize (int codim) const
    Return size of the ghost region for a given codim on the grid view.
  template<class DataHandleImp , class DataType >
  void 	communicate (CommDataHandleIF< DataHandleImp, DataType > &data, InterfaceType iftype, CommunicationDirection dir) const
    Communicate data on this view.
  }
};

typedef testing::Types<YASPGRIDS
#if HAVE_ALUGRID
                       ,
                       ALUGRIDS
#endif
                       > GridTypes;

TYPED_TEST_CASE(PeriodicViewTest, GridTypes);
TYPED_TEST(CubeGridProvider, is_default_creatable)
{
  this->is_default_creatable();
}
TYPED_TEST(PeriodicViewTest, fulfills_const_interface)
{
  this->const_interface();
}
TYPED_TEST(PeriodicViewTest, is_visualizable)
{
  this->visualize();
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
