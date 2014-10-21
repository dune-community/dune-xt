// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "main.hxx"

#if HAVE_DUNE_GRID

#include <iostream>
#include <fstream>
#include <utility>

#include <boost/filesystem.hpp>

#include <dune/common/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/grid/provider/interface.hh>
#include <dune/stuff/grid/provider/cube.hh>

using namespace Dune;
using namespace Dune::Stuff;

typedef testing::Types<SGrid<1, 1>, SGrid<2, 2>, SGrid<3, 3>, SGrid<4, 4>, SGrid<1, 2>, SGrid<2, 3>, SGrid<3, 4>,
                       SGrid<4, 5>, Dune::YaspGrid<1>, Dune::YaspGrid<2>, Dune::YaspGrid<3>> GridTypes;

template <class GridType>
struct GridProviderBaseTest : public testing::Test
{
  typedef Stuff::Grid::ProviderInterface<GridType> GridProviderType;

  template <class GridProviderImp>
  static void can_be_created_by_config()
  {
    std::unique_ptr<GridProviderImp> grid_provider = GridProviderImp::create();
  }

  static void fulfills_interface(GridProviderType& grid_provider)
  {
    using Stuff::Grid::ChoosePartView;
    std::string static_id = grid_provider.static_id();
    const auto& grid = grid_provider.grid();
    typename GridProviderType::template Leaf<ChoosePartView::view>::Type DUNE_UNUSED(leaf_view_auto) =
        grid_provider.template leaf<ChoosePartView::view>();
    typename GridProviderType::LeafGridViewType DUNE_UNUSED(leaf_view) = grid_provider.leaf_view();
#if HAVE_DUNE_FEM
    typename GridProviderType::template Leaf<ChoosePartView::part>::Type DUNE_UNUSED(leaf_part_auto) =
        grid_provider.template leaf<ChoosePartView::part>();
    typename GridProviderType::LeafGridPartType leaf_part = grid_provider.leaf_part();
#endif // HAVE_DUNE_FEM
    for (int level = 0; level <= grid.maxLevel(); ++level) {
      typename GridProviderType::template Level<ChoosePartView::view>::Type DUNE_UNUSED(level_view_auto) =
          grid_provider.template level<ChoosePartView::view>(level);
      typename GridProviderType::LevelGridViewType DUNE_UNUSED(level_view) = grid_provider.level_view(level);
#if HAVE_DUNE_FEM
      typename GridProviderType::template Level<ChoosePartView::part>::Type DUNE_UNUSED(level_part_auto) =
          grid_provider.template level<ChoosePartView::part>(level);
      typename GridProviderType::LevelGridPartType DUNE_UNUSED(level_part) = grid_provider.level_part(level);
#endif // HAVE_DUNE_FEM
      //! TODO: visualization for grid with dim > 3 ??
      if (GridType::dimension <= 3)
        grid_provider.visualize();
    }
  } // ... fulfills_interface()
}; // struct GridProviderBaseTest


template <class GridType>
struct CubeProviderTest : public GridProviderBaseTest<GridType>
{
  typedef GridProviderBaseTest<GridType> BaseType;
  typedef Stuff::Grid::Providers::Cube<GridType> CubeProviderType;

  void can_be_created_by_config() const
  {
    BaseType::template can_be_created_by_config<CubeProviderType>();
  }

  void fulfills_interface() const
  {
    std::string static_id = CubeProviderType::static_id();
    BaseType::fulfills_interface(*(CubeProviderType::create()));
  }
}; // struct CubeProviderTest


TYPED_TEST_CASE(CubeProviderTest, GridTypes);
TYPED_TEST(CubeProviderTest, can_be_created_by_config)
{
  this->can_be_created_by_config();
}

TYPED_TEST_CASE(CubeProviderTest, GridTypes);
TYPED_TEST(CubeProviderTest, fulfills_interface)
{
  this->fulfills_interface();
}

#endif // #if HAVE_DUNE_GRID
