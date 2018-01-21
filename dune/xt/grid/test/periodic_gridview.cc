// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2015 - 2018)
//   Tobias Leibner  (2015 - 2017)

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/type_traits.hh>

#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/view/periodic.hh>
#include <dune/xt/grid/gridprovider.hh>

using namespace Dune::XT;
using namespace Dune::XT::Common;


struct PeriodicViewTest : public testing::Test
{
  static const bool codim_iters_provided = CODIM_ITERS_PROVIDED;
  typedef TESTGRIDTYPE GridType;
  typedef typename GridType::ctype ctype;
  typedef typename GridType::template Codim<0>::Geometry GeometryType;
  typedef typename GridType::LeafGridView GridViewType;
  typedef typename GridViewType::template Codim<0>::Geometry::GlobalCoordinate DomainType;
  typedef typename Dune::XT::Grid::template PeriodicGridView<GridViewType, codim_iters_provided> PeriodicGridViewType;
  typedef typename PeriodicGridViewType::IndexSet IndexSet;
  typedef typename PeriodicGridViewType::template Codim<0>::Entity EntityType;
  typedef typename PeriodicGridViewType::template Codim<0>::Iterator EntityIteratorType;
  typedef typename Dune::XT::Grid::internal::template PeriodicIntersectionIterator<GridViewType>
      PeriodicIntersectionIteratorType;
  typedef typename PeriodicGridViewType::Intersection PeriodicIntersectionType;
  typedef typename GridViewType::CollectiveCommunication CollectiveCommunication;
  static const size_t dimDomain = GridViewType::dimension;

  static constexpr int factorial(int n)
  {
    return n > 0 ? n * factorial(n - 1) : 1;
  }

  static void check()
  {
    // create grid and get gridview
    Configuration grid_config = DXTC_CONFIG.sub("test_grid_periodicview");
    const bool is_simplex = Common::from_string<bool>(grid_config["is_simplex"]);
    const bool is_cube = !is_simplex;
    const DomainType lower_left = Common::from_string<DomainType>(grid_config["lower_left"]);
    const DomainType upper_right = Common::from_string<DomainType>(grid_config["upper_right"]);
    const auto num_elements = Common::from_string<std::array<unsigned int, dimDomain>>(grid_config["num_elements"]);
    for (const auto& elements : num_elements)
      if (elements != num_elements[0])
        DUNE_THROW(Dune::NotImplemented,
                   "Test works only if the grid has the same number of elements in each direction!");
    const int elements_per_direction = num_elements[0];
    std::shared_ptr<GridType> grid;
    if (is_cube)
      grid = Dune::StructuredGridFactory<GridType>::createCubeGrid(lower_left, upper_right, num_elements);
    else
      grid = Dune::StructuredGridFactory<GridType>::createSimplexGrid(lower_left, upper_right, num_elements);
    const GridViewType grid_view = grid->leafGridView();

    // check whether grid is periodic
    const bool is_nonperiodic = grid_config["periodicity"] == "nonperiodic";
    const bool is_partially_periodic = grid_config["periodicity"] == "partially_periodic";

    // create periodic grid_view
    std::bitset<dimDomain> periodic_directions;
    if (is_partially_periodic)
      periodic_directions[0] = 1;
    else if (!is_nonperiodic)
      periodic_directions.set();
    const PeriodicGridViewType periodic_grid_view(grid_view, periodic_directions);

    // check interface
    const GridType& test_grid = periodic_grid_view.grid();
    (void)test_grid;
    const IndexSet& index_set = periodic_grid_view.indexSet();
    const int codim0_size = periodic_grid_view.size(0);
    EXPECT_EQ(grid_view.size(0), codim0_size);
    if (is_cube)
      EXPECT_EQ(std::pow(elements_per_direction, dimDomain), codim0_size);
    if (is_simplex)
      EXPECT_EQ(std::pow(elements_per_direction, dimDomain) * factorial(dimDomain), codim0_size);
    EXPECT_EQ(grid_view.overlapSize(0), periodic_grid_view.overlapSize(0));
    EXPECT_EQ(grid_view.overlapSize(1), periodic_grid_view.overlapSize(1));
    EXPECT_EQ(grid_view.ghostSize(0), periodic_grid_view.ghostSize(0));
    EXPECT_EQ(grid_view.ghostSize(1), periodic_grid_view.ghostSize(1));
    const CollectiveCommunication& test_comm = periodic_grid_view.comm();
    (void)test_comm;


    size_t neighbor_count = 0;
    size_t boundary_count = 0;
    size_t periodic_count = 0;
    // iterate over codim 0 entities
    const EntityIteratorType it_end = periodic_grid_view.template end<0>();
    for (EntityIteratorType it = periodic_grid_view.template begin<0>(); it != it_end; ++it) {
      const EntityType& entity = *it;
      EXPECT_TRUE(periodic_grid_view.contains(entity));
      EXPECT_TRUE(index_set.contains(entity));
      const auto sub_index = index_set.subIndex(entity, 0, 1);
      (void)sub_index;
      // iterate over all intersections on current entity
      const PeriodicIntersectionIteratorType i_it_end = periodic_grid_view.iend(entity);
      for (PeriodicIntersectionIteratorType i_it = periodic_grid_view.ibegin(entity); i_it != i_it_end; ++i_it) {
        const PeriodicIntersectionType& intersection = *i_it;
        if (intersection.neighbor()) {
          ++neighbor_count;
          const EntityType outside = intersection.outside();
          // find corresponding intersection in outside
          const auto index_in_outside = intersection.indexInOutside();
          PeriodicIntersectionIteratorType i_it_outside = periodic_grid_view.ibegin(outside);
          const PeriodicIntersectionIteratorType i_it_outside_end = periodic_grid_view.iend(outside);
          for (; i_it_outside != i_it_outside_end; ++i_it_outside) {
            const PeriodicIntersectionType* outside_intersection = i_it_outside.operator->();
            if (outside_intersection->indexInInside() == index_in_outside) {
              break;
            }
          }
          const PeriodicIntersectionType& intersection_in_outside = *i_it_outside;
          // check outside_intersection coords
          const auto coords_in_outside = intersection.geometryInOutside().center();
          const auto coords_in_outside_2 = intersection_in_outside.geometryInInside().center();
          EXPECT_TRUE(Dune::XT::Common::FloatCmp::eq(coords_in_outside, coords_in_outside_2, 1e-14));
          // check global intersection coords in periodic case
          const auto global_intersection_coords = intersection.geometry().center();
          const auto global_outside_intersection_coords = intersection_in_outside.geometry().center();
          size_t coord_difference_count = 0;
          size_t differing_coordinate;
          for (size_t ii = 0; ii < dimDomain; ++ii) {
            if (Dune::XT::Common::FloatCmp::ne(global_outside_intersection_coords[ii],
                                               global_intersection_coords[ii])) {
              ++coord_difference_count;
              differing_coordinate = ii;
            }
          }
          if (intersection.boundary() && intersection.neighbor()) {
            EXPECT_TRUE(intersection_in_outside.boundary() && intersection_in_outside.neighbor());
            EXPECT_EQ(size_t(1), coord_difference_count);
            EXPECT_TRUE((Dune::XT::Common::FloatCmp::eq(global_outside_intersection_coords[differing_coordinate],
                                                        lower_left[differing_coordinate])
                         && Dune::XT::Common::FloatCmp::eq(global_intersection_coords[differing_coordinate],
                                                           upper_right[differing_coordinate]))
                        || (Dune::XT::Common::FloatCmp::eq(global_outside_intersection_coords[differing_coordinate],
                                                           upper_right[differing_coordinate])
                            && Dune::XT::Common::FloatCmp::eq(global_intersection_coords[differing_coordinate],
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
    /* on each face, there are elements_per_direction**(dimDomain-1) intersections. For a simplex grid in 3 dimensions,
     * there are twice as
     * much. */
    size_t num_intersections_on_face = std::pow(elements_per_direction, dimDomain - 1);
    assert(dimDomain == Common::from_string<int>(grid_config["dimDomain"]));
    const auto domainDim = Common::from_string<int>(grid_config["dimDomain"]);
    if (is_simplex && domainDim == 3) // use dimDomain from config here to avoid "code will never be executed" warning
      num_intersections_on_face *= 2;
    /* In a fully periodic grid, all intersections are periodic. In a partially periodic grid, only the intersections on
     * two faces are periodic. In a nonperiodic grid, no intersections are periodic. */
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

    // the nonperiodic grid has (elements_per_direction-1)**dimDomain inner vertices
    size_t expected_num_vertices = std::pow(elements_per_direction - 1, dimDomain);
    // add number of vertices on faces (codim 1)
    expected_num_vertices += std::pow(elements_per_direction - 1, dimDomain - 1) * (num_faces - num_periodic_faces / 2);
    // add number of vertices on edges (codim 2)
    const size_t num_edges = dimDomain == 1 ? 0 : (dimDomain == 2 ? 4 : 12);
    size_t num_periodic_edges = is_partially_periodic ? num_periodic_faces * std::pow(2, dimDomain - 1) : num_edges;
    if (is_nonperiodic)
      num_periodic_edges = 0;
    expected_num_vertices +=
        dimDomain == 1
            ? 0
            : std::pow(elements_per_direction - 1, dimDomain - 2)
                  * ((num_edges - num_periodic_edges) + num_periodic_edges / (is_partially_periodic ? 2 : 4));
    // add vertices on corners (codim 3) of grid
    if (domainDim == 3)
      expected_num_vertices += is_partially_periodic ? 4 : (is_nonperiodic ? 8 : 1);

    const size_t num_vertices = periodic_grid_view.size(dimDomain);
    EXPECT_EQ(expected_num_vertices, num_vertices);

    // check Codim iterator
    size_t codim_iterator_size = 0;
    const auto codim_it_end = periodic_grid_view.template end<dimDomain>();
    for (auto codim_it = periodic_grid_view.template begin<dimDomain>(); codim_it != codim_it_end; ++codim_it) {
      const auto codim_it_copy DUNE_UNUSED = codim_it;
      const auto& entity DUNE_UNUSED = *codim_it;
      ++codim_iterator_size;
    }
    EXPECT_EQ(num_vertices, codim_iterator_size);

    // check that new indices are zero-starting and consecutive
    std::vector<size_t> index_counter(num_vertices);
    const auto real_codim_it_end = grid_view.template end<dimDomain>();
    for (auto codim_it = grid_view.template begin<dimDomain>(); codim_it != real_codim_it_end; ++codim_it) {
      const auto& entity = *codim_it;
      ++(index_counter[index_set.index(entity)]);
    }
    for (auto& count : index_counter)
      EXPECT_GT(count, 0);


  } // void check(...)
};

TEST_F(PeriodicViewTest, check_all)
{
  this->check();
}
