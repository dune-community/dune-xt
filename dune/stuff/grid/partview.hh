// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_GRID_PARTVIEW_HH
#define DUNE_STUFF_GRID_PARTVIEW_HH

#include <memory>


#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/gridpart/levelgridpart.hh>
#endif // HAVE_DUNE_FEM

namespace Dune {
namespace Stuff {
namespace Grid {


enum class ChoosePartView
{
  part,
  view
}; // enum class ChoosePartView


#if HAVE_DUNE_GRID


// forwards
template <class GridType, ChoosePartView type>
struct LevelPartView;

template <class GridType, ChoosePartView type>
struct LeafPartView;


template <class GridType>
struct LevelPartView<GridType, ChoosePartView::view>
{
  typedef typename GridType::LevelGridView Type;

  static std::shared_ptr<Type> create(GridType& grid, const int level)
  {
    assert(level >= 0);
    assert(level <= grid.maxLevel());
    return std::make_shared<Type>(grid.levelGridView(level));
  } // ... create(...)
}; // struct LevelPartView< ..., view >


template <class GridType>
struct LeafPartView<GridType, ChoosePartView::view>
{
  typedef typename GridType::LeafGridView Type;

  static std::shared_ptr<Type> create(GridType& grid)
  {
    return std::make_shared<Type>(grid.leafGridView());
  }
}; // struct LeafPartView< ..., view >


#if HAVE_DUNE_FEM


template <class GridType>
struct LevelPartView<GridType, ChoosePartView::part>
{
  typedef Dune::Fem::LevelGridPart<GridType> Type;

  static std::shared_ptr<Type> create(GridType& grid, const int level)
  {
    assert(level >= 0);
    assert(level <= grid.maxLevel());
    return std::make_shared<Type>(grid, level);
  } // ... create(...)
}; // struct LevelPartView< ..., part >


template <class GridType>
struct LeafPartView<GridType, ChoosePartView::part>
{
  typedef Dune::Fem::LeafGridPart<GridType> Type;

  static std::shared_ptr<Type> create(GridType& grid)
  {
    return std::make_shared<Type>(grid);
  }
}; // struct LeafPartView< ..., part >


#endif // HAVE_DUNE_FEM
#else // HAVE_DUNE_GRID


template <class GridType, ChoosePartView type>
struct LevelPartView
{
  static_assert(AlwaysFalse<GridType>::value, "You are missing dune-grid!");
};

template <class GridType, ChoosePartView type>
struct LeafPartView
{
  static_assert(AlwaysFalse<GridType>::value, "You are missing dune-grid!");
};


#endif // HAVE_DUNE_GRID

} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_PARTVIEW_HH
