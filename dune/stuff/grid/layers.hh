// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_GRID_LAYERS_HH
#define DUNE_STUFF_GRID_LAYERS_HH

#include <memory>

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/gridpart/levelgridpart.hh>
#endif

#if HAVE_DUNE_GRID_MULTISCALE
#include <dune/grid/multiscale/default.hh>
#endif

#include <dune/stuff/common/type_utils.hh>


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


// specializations of LeafPartView

/**
 * \brief Allows to statically create a leaf part or view (view variant).
 */
template <class GridType>
struct LeafPartView<GridType, ChoosePartView::view>
{
  typedef typename Layer<GridType, ChooseLayer::leaf, ChoosePartView::view>::Type Type;

  static Type create(const GridType& grid, const int /*level*/ = 0)
  {
    return grid.leafGridView();
  }
}; // struct LeafPartView< ..., view >


#if HAVE_DUNE_FEM

/**
 * \brief Allows to statically create a leaf part or view (part variant, only from a non-const grid).
 */
template <class GridType>
struct LeafPartView<GridType, ChoosePartView::part>
{
  typedef typename Layer<GridType, ChooseLayer::leaf, ChoosePartView::part>::Type Type;

  static Type create(const GridType& /*grid*/, const int /*level*/ = 0)
  {
    static_assert(AlwaysFalse<GridType>::value,
                  "dune-fem does not allow the creation of a leaf grid part from a non-const grid!");
  }

  static Type create(GridType& grid, const int /*level*/ = 0)
  {
    return Type(grid);
  }
}; // struct LeafPartView< ..., part >

#else // HAVE_DUNE_FEM

template <class GridType>
struct LeafPartView<GridType, ChoosePartView::part>
{
  static_assert(AlwaysFalse<GridType>::value, "You are missing dune-fem!");
};

#endif // HAVE_DUNE_FEM


// specializations of LevelPartView

/**
 * \brief Allows to statically create a level part or view (view variant).
 */
template <class GridType>
struct LevelPartView<GridType, ChoosePartView::view>
{
  typedef typename Layer<GridType, ChooseLayer::level, ChoosePartView::view>::Type Type;

  static Type create(const GridType& grid, const int level)
  {
    assert(level >= 0);
    assert(level <= grid.maxLevel());
    return grid.levelGridView(level);
  }
}; // struct LevelPartView< ..., view >

#if HAVE_DUNE_FEM

/**
 * \brief Allows to statically create a level part or view (part variant, only from a non-const grid).
 */
template <class GridType>
struct LevelPartView<GridType, ChoosePartView::part>
{
  typedef typename Layer<GridType, ChooseLayer::level, ChoosePartView::part>::Type Type;

  static Type create(const GridType& /*grid*/, const int /*level*/)
  {
    static_assert(AlwaysFalse<GridType>::value,
                  "dune-fem does not allow the creation of a level grid part from a non-const grid!");
  }

  static Type create(GridType& grid, const int level)
  {
    assert(level >= 0);
    assert(level <= grid.maxLevel());
    return Type(grid, level);
  }
}; // struct LevelPartView< ..., part >


// specializations of LayerPart

/**
 * \brief Allows to statically create a leaf or level part (leaf variant, only from a non-const grid).
 */
template <class GridType>
struct LayerPart<GridType, ChooseLayer::leaf>
{
  typedef typename Layer<GridType, ChooseLayer::leaf, ChoosePartView::part>::Type Type;

  static Type create(const GridType& /*grid*/, const int /*level*/ = 0)
  {
    static_assert(AlwaysFalse<GridType>::value,
                  "dune-fem does not allow the creation of a leaf grid part from a non-const grid!");
  }

  static Type create(GridType& grid, const int /*level*/ = 0)
  {
    return LeafPartView<GridType, ChoosePartView::part>::create(grid);
  }
}; // struct LayerPart< ..., leaf >

/**
 * \brief Allows to statically create a leaf or level part (level variant, only from a non-const grid).
 */
template <class GridType>
struct LayerPart<GridType, ChooseLayer::level>
{
  typedef typename Layer<GridType, ChooseLayer::level, ChoosePartView::part>::Type Type;

  static Type create(const GridType& /*grid*/, const int /*level*/)
  {
    static_assert(AlwaysFalse<GridType>::value,
                  "dune-fem does not allow the creation of a level grid part from a non-const grid!");
  }

  static Type create(GridType& grid, const int level)
  {
    return LevelPartView<GridType, ChoosePartView::part>::create(grid, level);
  }
}; // struct LayerPart< ..., level >

#else // HAVE_DUNE_FEM

template <class GridType>
struct LevelPartView<GridType, ChoosePartView::part>
{
  static_assert(AlwaysFalse<GridType>::value, "You are missing dune-fem!");
};

template <class GridType>
struct LayerPart<GridType, ChooseLayer::leaf>
{
  static_assert(AlwaysFalse<GridType>::value, "You are missing dune-fem!");
};

template <class GridType>
struct LayerPart<GridType, ChooseLayer::level>
{
  static_assert(AlwaysFalse<GridType>::value, "You are missing dune-fem!");
};

#endif // HAVE_DUNE_FEM


// specializations of LayerView

/**
 * \brief Allows to statically create a leaf or level view (leaf variant).
 */
template <class GridType>
struct LayerView<GridType, ChooseLayer::leaf>
{
  typedef typename Layer<GridType, ChooseLayer::leaf, ChoosePartView::view>::Type Type;

  static Type create(const GridType& grid, const int /*level*/ = 0)
  {
    return LeafPartView<GridType, ChoosePartView::view>::create(grid);
  }
}; // struct LayerView< ..., leaf >

/**
 * \brief Allows to statically create a leaf or level view (level variant).
 */
template <class GridType>
struct LayerView<GridType, ChooseLayer::level>
{
  typedef typename Layer<GridType, ChooseLayer::level, ChoosePartView::view>::Type Type;

  static Type create(const GridType& grid, const int level)
  {
    return LevelPartView<GridType, ChoosePartView::view>::create(grid, level);
  }
}; // struct LayerView< ..., level >


// specializatins of Layer

#if HAVE_DUNE_FEM

/**
 * \brief Allows to statically create a leaf or level part or view (leaf part variant, only from a non-const grid).
 */
template <class GridType>
struct Layer<GridType, ChooseLayer::leaf, ChoosePartView::part>
{
  typedef Dune::Fem::LeafGridPart<GridType> Type;

  static Type create(const GridType& /*grid*/, const int /*level*/ = 0)
  {
    static_assert(AlwaysFalse<GridType>::value,
                  "dune-fem does not allow the creation of a leaf grid part from a non-const grid!");
  }

  static Type create(GridType& grid, const int /*level*/ = 0)
  {
    return LeafPartView<GridType, ChoosePartView::part>::create(grid);
  }
}; // struct Layer< ..., leaf, part >

#else // HAVE_DUNE_FEM

template <class GridType>
struct Layer<GridType, ChooseLayer::leaf, ChoosePartView::part>
{
  static_assert(AlwaysFalse<GridType>::value, "You are missing dune-fem!");
};

#endif // HAVE_DUNE_FEM

/**
 * \brief Allows to statically create a leaf or level part or view (leaf view variant).
 */
template <class GridType>
struct Layer<GridType, ChooseLayer::leaf, ChoosePartView::view>
{
  typedef typename GridType::LeafGridView Type;

  static Type create(const GridType& grid, const int /*level*/ = 0)
  {
    return LeafPartView<GridType, ChoosePartView::view>::create(grid);
  }
}; // struct Layer< ..., leaf, view >

#if HAVE_DUNE_FEM

/**
 * \brief Allows to statically create a leaf or level part or view (level part variant, only from a non-const grid).
 */
template <class GridType>
struct Layer<GridType, ChooseLayer::level, ChoosePartView::part>
{
  typedef Dune::Fem::LevelGridPart<GridType> Type;

  static Type create(const GridType& /*grid*/, const int /*level*/)
  {
    static_assert(AlwaysFalse<GridType>::value,
                  "dune-fem does not allow the creation of a level grid part from a non-const grid!");
  }

  static Type create(GridType& grid, const int level)
  {
    return LevelPartView<GridType, ChoosePartView::part>::create(grid, level);
  }
}; // struct Layer< ..., level, part >

#else // HAVE_DUNE_FEM

template <class GridType>
struct Layer<GridType, ChooseLayer::level, ChoosePartView::part>
{
  static_assert(AlwaysFalse<GridType>::value, "You are missing dune-fem!");
};

#endif // HAVE_DUNE_FEM

/**
 * \brief Allows to statically create a leaf or level part or view (leaf view variant).
 */
template <class GridType>
struct Layer<GridType, ChooseLayer::level, ChoosePartView::view>
{
  typedef typename GridType::LevelGridView Type;

  static Type create(const GridType& grid, const int level)
  {
    return LevelPartView<GridType, ChoosePartView::view>::create(grid, level);
  }
}; // struct Layer< ..., level, view >


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
