// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2014 - 2017, 2019)
//   René Fritze     (2014 - 2019)
//   Tobias Leibner  (2017 - 2018, 2020)

#ifndef DUNE_XT_GRID_LAYERS_HH
#define DUNE_XT_GRID_LAYERS_HH

#include <dune/common/typetraits.hh>

#include <dune/xt/common/fixed_map.hh>
#include <dune/xt/grid/type_traits.hh>

namespace Dune {
namespace XT {
namespace Grid {


// forward
template <class RealGridLayerImp, bool codim_iters_provided>
class PeriodicGridView;


enum class Backends
{
  view
};


enum class Layers
{
  adaptive_leaf,
  leaf,
  level
};

namespace {
const XT::Common::FixedMap<Layers, std::string, 3> layer_names{
    {Layers::adaptive_leaf, "adaptive_leaf"}, {Layers::leaf, "leaf"}, {Layers::level, "level"}};
}

namespace internal {


template <Backends backend>
struct backend_dependent_typename
{
  using type = void;
};


template <Layers layer>
struct layer_dependent_typename
{
  using type = void;
};


} // namespace  internal


/**
 * \brief Allows to statically create a leaf or level part or view (unspecialized variant).
 */
template <class GridType, Layers layer, Backends backend, bool periodic = false>
struct Layer
{
  static_assert(AlwaysFalse<GridType>::value, "No layer available for this combination!");

  using type = void;

  static type create(const GridType& /*grid*/, const int /*level*/ = 0)
  {
    static_assert(AlwaysFalse<GridType>::value, "No layer available for this combination!");
  }

  static type create(GridType& /*grid*/, const int /*level*/ = 0)
  {
    static_assert(AlwaysFalse<GridType>::value, "No layer available for this combination!");
  }
}; // struct Layer


/**
 * \brief Allows to statically create a leaf or level part or view (leaf view variant).
 */
template <class GridType>
struct Layer<GridType, Layers::leaf, Backends::view, false>
{
  using type = typename GridType::LeafGridView;

  static type create(const GridType& grid, const int /*level*/ = 0)
  {
    return grid.leafGridView();
  }
}; // struct Layer<..., leaf, view>


/**
 * \brief Allows to statically create a leaf or level part or view (leaf view variant).
 */
template <class GridType>
struct Layer<GridType, Layers::level, Backends::view, false>
{
  using type = typename GridType::LevelGridView;

  static type create(const GridType& grid, const int level)
  {
    assert(level >= 0);
    assert(level <= grid.maxLevel());
    return grid.levelGridView(level);
  }
}; // struct Layer<..., level, view>


// create periodic grid_layer
template <class GridType, Layers layer, Backends backend>
struct Layer<GridType, layer, backend, true>
{
  using NonPeriodicLayerType = Layer<GridType, layer, backend, false>;
  using type = PeriodicGridView<typename NonPeriodicLayerType::type, false>;

  static type create(const GridType& grid, const int subdomain = 0)
  {
    return type(NonPeriodicLayerType::create(grid, subdomain));
  }

  static type create(GridType& grid, const int subdomain)
  {
    return type(NonPeriodicLayerType::create(grid, subdomain));
  }
}; // struct Layer<..., periodic>


template <class GridLayerType>
struct extract_layer_backend
{
  static constexpr Backends value = Backends::view;
};


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_LAYERS_HH
