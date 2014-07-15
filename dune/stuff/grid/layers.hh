// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_GRID_LAYERS_HH
#define DUNE_STUFF_GRID_LAYERS_HH

#include <memory>

#if HAVE_DUNE_FEM
#include <dune/stuff/common/disable_warnings.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/gridpart/levelgridpart.hh>
#include <dune/stuff/common/reenable_warnings.hh>
#endif // HAVE_DUNE_FEM

#if HAVE_DUNE_GRID_MULTISCALE
#include <dune/grid/multiscale/default.hh>
#endif


namespace Dune {
namespace Stuff {
namespace Grid {


enum class ChoosePartView
{
  view,
  part
}; // enum class ChoosePartView


enum class ChooseLayer
{
  level,
  leaf
#if HAVE_DUNE_GRID_MULTISCALE
  ,
  local,
  local_oversampled
#endif
}; // enum class ChooseLayer


#if HAVE_DUNE_GRID


// forwards
template <class GridType, ChooseLayer layer, ChoosePartView part_view>
struct Layer;

template <class GridType, ChoosePartView type>
struct LevelPartView;

template <class GridType, ChoosePartView type>
struct LeafPartView;

template <class GridType, ChooseLayer type>
struct LayerView;

template <class GridType, ChooseLayer type>
struct LayerPart;


template <class GridType>
struct Layer<GridType, ChooseLayer::level, ChoosePartView::view>
{
  typedef typename GridType::LevelGridView Type;

  static std::shared_ptr<Type> create(GridType& grid, const int level)
  {
    assert(level >= 0);
    assert(level <= grid.maxLevel());
    return std::make_shared<Type>(grid.levelGridView(level));
  } // ... create(...)
}; // struct Layer< ..., level, view >


template <class GridType>
struct Layer<GridType, ChooseLayer::leaf, ChoosePartView::view>
{
  typedef typename GridType::LeafGridView Type;

  static std::shared_ptr<Type> create(GridType& grid, const int /*level*/ = 0)
  {
    return std::make_shared<Type>(grid.leafGridView());
  }
}; // struct Layer< ..., leaf, view >


#if HAVE_DUNE_GRID_MULTISCALE


template <class GridType>
struct Layer<GridType, ChooseLayer::local, ChoosePartView::view>
{
  typedef typename grid::Multiscale::Default<GridType>::LocalGridViewType Type;
};


template <class GridType>
struct Layer<GridType, ChooseLayer::local_oversampled, ChoosePartView::view>
{
  typedef typename grid::Multiscale::Default<GridType>::LocalGridViewType Type;
};


#endif // HAVE_DUNE_GRID_MULTISCALE


template <class GridType>
struct LevelPartView<GridType, ChoosePartView::view>
{
  typedef typename Layer<GridType, ChooseLayer::level, ChoosePartView::view>::Type Type;

  static std::shared_ptr<Type> create(GridType& grid, const int level)
  {
    return Layer<GridType, ChooseLayer::level, ChoosePartView::view>::create(grid, level);
  }
}; // struct LevelPartView< ..., view >


template <class GridType>
struct LeafPartView<GridType, ChoosePartView::view>
{
  typedef typename Layer<GridType, ChooseLayer::leaf, ChoosePartView::view>::Type Type;

  static std::shared_ptr<Type> create(GridType& grid, const int /*level*/ = 0)
  {
    return Layer<GridType, ChooseLayer::leaf, ChoosePartView::view>::create(grid);
  }
}; // struct LeafPartView< ..., view >


template <class GridType>
struct LayerView<GridType, ChooseLayer::level>
{
  typedef typename Layer<GridType, ChooseLayer::level, ChoosePartView::view>::Type Type;

  static std::shared_ptr<Type> create(GridType& grid, const int level)
  {
    return Layer<GridType, ChooseLayer::level, ChoosePartView::view>::create(grid, level);
  }
}; // struct LayerView< ..., level >


template <class GridType>
struct LayerView<GridType, ChooseLayer::leaf>
{
  typedef typename Layer<GridType, ChooseLayer::leaf, ChoosePartView::view>::Type Type;

  static std::shared_ptr<Type> create(GridType& grid, const int /*level*/ = 0)
  {
    return Layer<GridType, ChooseLayer::level, ChoosePartView::view>::create(grid);
  }
}; // struct LayerView< ..., leaf >


#if HAVE_DUNE_FEM


template <class GridType>
struct Layer<GridType, ChooseLayer::level, ChoosePartView::part>
{
  typedef Dune::Fem::LevelGridPart<GridType> Type;

  static std::shared_ptr<Type> create(GridType& grid, const int level)
  {
    assert(level >= 0);
    assert(level <= grid.maxLevel());
    return std::make_shared<Type>(grid, level);
  } // ... create(...)
}; // struct Layer< ..., level, part >


template <class GridType>
struct Layer<GridType, ChooseLayer::leaf, ChoosePartView::part>
{
  typedef Dune::Fem::LeafGridPart<GridType> Type;

  static std::shared_ptr<Type> create(GridType& grid, const int /*level*/ = 0)
  {
    return std::make_shared<Type>(grid);
  }
}; // struct Layer< ..., leaf, part >


#if HAVE_DUNE_GRID_MULTISCALE


template <class GridType>
struct Layer<GridType, ChooseLayer::local, ChoosePartView::part>
{
  typedef typename grid::Multiscale::Default<GridType>::LocalGridPartType Type;
};


template <class GridType>
struct Layer<GridType, ChooseLayer::local_oversampled, ChoosePartView::part>
{
  typedef typename grid::Multiscale::Default<GridType>::LocalGridPartType Type;
};


#endif // HAVE_DUNE_GRID_MULTISCALE


template <class GridType>
struct LevelPartView<GridType, ChoosePartView::part>
{
  typedef typename Layer<GridType, ChooseLayer::level, ChoosePartView::part>::Type Type;

  static std::shared_ptr<Type> create(GridType& grid, const int level)
  {
    return Layer<GridType, ChooseLayer::level, ChoosePartView::part>::create(grid, level);
  }
}; // struct LevelPartView< ..., part >


template <class GridType>
struct LeafPartView<GridType, ChoosePartView::part>
{
  typedef typename Layer<GridType, ChooseLayer::leaf, ChoosePartView::part>::Type Type;

  static std::shared_ptr<Type> create(GridType& grid, const int /*level*/ = 0)
  {
    return Layer<GridType, ChooseLayer::leaf, ChoosePartView::part>::create(grid);
  }
}; // struct LeafPartView< ..., part >


template <class GridType>
struct LayerPart<GridType, ChooseLayer::level>
{
  typedef typename Layer<GridType, ChooseLayer::level, ChoosePartView::part>::Type Type;

  static std::shared_ptr<Type> create(GridType& grid, const int level)
  {
    return Layer<GridType, ChooseLayer::level, ChoosePartView::part>::create(grid, level);
  }
}; // struct LayerPart< ..., level >


template <class GridType>
struct LayerPart<GridType, ChooseLayer::leaf>
{
  typedef typename Layer<GridType, ChooseLayer::leaf, ChoosePartView::part>::Type Type;

  static std::shared_ptr<Type> create(GridType& grid, const int /*level*/ = 0)
  {
    return Layer<GridType, ChooseLayer::leaf, ChoosePartView::part>::create(grid);
  }
}; // struct LayerPart< ..., leaf >


#endif // HAVE_DUNE_FEM
#else // HAVE_DUNE_GRID


template <class GridType, ChooseLayer layer, ChoosePartView type>
struct Layer
{
  static_assert(AlwaysFalse<GridType>::value, "You are missing dune-grid!");
};


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

#endif // DUNE_STUFF_GRID_LAYERS_HH
