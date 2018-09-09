// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014 - 2018)
//   Tobias Leibner  (2017 - 2018)

#ifndef DUNE_XT_GRID_LAYERS_HH
#define DUNE_XT_GRID_LAYERS_HH

#include <dune/common/typetraits.hh>

#include <dune/xt/common/fixed_map.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/view/periodic.hh>
#include <dune/xt/grid/view/subdomain/view.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace DD {


// forward
template <class GridImp>
class SubdomainGrid;


} // namespace DD


enum class Backends
{
  view
};


enum class Layers
{
  adaptive_leaf,
  leaf,
  level,
  dd_subdomain,
  dd_subdomain_oversampled,
  dd_subdomain_boundary,
  dd_subdomain_coupling
};

namespace {
const XT::Common::FixedMap<Layers, std::string, 7> layer_names{
    {Layers::adaptive_leaf, "adaptive_leaf"},
    {Layers::leaf, "leaf"},
    {Layers::level, "level"},
    {Layers::dd_subdomain, "dd_subdomain"},
    {Layers::dd_subdomain_oversampled, "dd_subdomain_oversampled"},
    {Layers::dd_subdomain_boundary, "dd_subdomain_boundary"},
    {Layers::dd_subdomain_coupling, "dd_subdomain_coupling"}};
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
template <class GridType, Layers layer, Backends backend, class DdGridType = int, bool periodic = false>
struct Layer
{
  static_assert(AlwaysFalse<GridType>::value, "No layer available for this combination!");

  using type = void;

  static type
  create(const GridType& /*grid*/, const int /*level*/ = 0, const std::shared_ptr<DdGridType> /*dd_grid*/ = nullptr)
  {
    static_assert(AlwaysFalse<GridType>::value, "No layer available for this combination!");
  }

  static type create(GridType& /*grid*/, const int /*level*/ = 0, std::shared_ptr<DdGridType> /*dd_grid*/ = nullptr)
  {
    static_assert(AlwaysFalse<GridType>::value, "No layer available for this combination!");
  }
}; // struct Layer


/**
 * \brief Allows to statically create a leaf or level part or view (leaf view variant).
 */
template <class GridType, class DdGridType>
struct Layer<GridType, Layers::leaf, Backends::view, DdGridType, false>
{
  using type = typename GridType::LeafGridView;

  static type
  create(const GridType& grid, const int /*level*/ = 0, const std::shared_ptr<DdGridType> /*dd_grid*/ = nullptr)
  {
    return grid.leafGridView();
  }
}; // struct Layer<..., leaf, view>


/**
 * \brief Allows to statically create a leaf or level part or view (leaf view variant).
 */
template <class GridType, class DdGridType>
struct Layer<GridType, Layers::level, Backends::view, DdGridType, false>
{
  using type = typename GridType::LevelGridView;

  static type create(const GridType& grid, const int level, const std::shared_ptr<DdGridType> /*dd_grid*/ = nullptr)
  {
    assert(level >= 0);
    assert(level <= grid.maxLevel());
    return grid.levelGridView(level);
  }
}; // struct Layer<..., level, view>


template <class GridType, class DdGridType>
struct Layer<GridType, Layers::dd_subdomain, Backends::view, DdGridType, false>
{
  using type = typename DD::SubdomainGrid<GridType>::LocalGridViewType;

  static type create(const GridType& /*grid*/,
                     const int subdomain = 0,
                     const std::shared_ptr<DD::SubdomainGrid<GridType>> dd_grid = nullptr)
  {
    static_assert(std::is_same<DdGridType, DD::SubdomainGrid<GridType>>::value,
                  "Only available for DD::SubdomainGrid!");
    return dd_grid->local_grid_view(subdomain, /*oversampling=*/false);
  }

  static type create(GridType& /*grid*/, const int subdomain, std::shared_ptr<DD::SubdomainGrid<GridType>> dd_grid)
  {
    static_assert(std::is_same<DdGridType, DD::SubdomainGrid<GridType>>::value,
                  "Only available for DD::SubdomainGrid!");
    return dd_grid->local_grid_view(subdomain, /*oversampling=*/false);
  }
}; // struct Layer<..., dd_subdomain, view>


template <class GridType, class DdGridType>
struct Layer<GridType, Layers::dd_subdomain_oversampled, Backends::view, DdGridType, false>
{
  using type = typename DD::SubdomainGrid<GridType>::LocalGridViewType;

  static type create(const GridType& /*grid*/,
                     const int subdomain = 0,
                     const std::shared_ptr<DD::SubdomainGrid<GridType>> dd_grid = nullptr)
  {
    static_assert(std::is_same<DdGridType, DD::SubdomainGrid<GridType>>::value,
                  "Only available for DD::SubdomainGrid!");
    return dd_grid->local_grid_view(subdomain, /*oversampling=*/true);
  }

  static type create(GridType& /*grid*/, const int subdomain, std::shared_ptr<DD::SubdomainGrid<GridType>> dd_grid)
  {
    static_assert(std::is_same<DdGridType, DD::SubdomainGrid<GridType>>::value,
                  "Only available for DD::SubdomainGrid!");
    return dd_grid->local_grid_view(subdomain, /*oversampling=*/true);
  }
}; // struct Layer<..., dd_subdomain_oversampled, view>


template <class GridType, class DdGridType>
struct Layer<GridType, Layers::dd_subdomain_coupling, Backends::view, DdGridType, false>
{
  using type = typename DD::SubdomainGrid<GridType>::CouplingGridViewType;

  static type create(const GridType& /*grid*/,
                     const int /*subdomain */ = 0,
                     const std::shared_ptr<DD::SubdomainGrid<GridType>> dd_grid = nullptr)
  {
    static_assert(std::is_same<DdGridType, DD::SubdomainGrid<GridType>>::value,
                  "Only available for DD::SubdomainGrid!");
    DUNE_THROW(NotImplemented, "Only usable to extract the layer type, not the actual layer!");
    return dd_grid->coupling_grid_view(0, 0);
  }

  static type create(GridType& /*grid*/, const int /*subdomain*/, std::shared_ptr<DD::SubdomainGrid<GridType>> dd_grid)
  {
    static_assert(std::is_same<DdGridType, DD::SubdomainGrid<GridType>>::value,
                  "Only available for DD::SubdomainGrid!");
    DUNE_THROW(NotImplemented, "Only usable to extract the layer type, not the actual layer!");
    return dd_grid->coupling_grid_view(0, 0);
  }
}; // struct Layer<..., dd_subdomain_coupling, view>


template <class GridType, class DdGridType>
struct Layer<GridType, Layers::dd_subdomain_boundary, Backends::view, DdGridType, false>
{
  using type = typename DD::SubdomainGrid<GridType>::BoundaryGridViewType;

  static type create(const GridType& /*grid*/,
                     const int subdomain = 0,
                     const std::shared_ptr<DD::SubdomainGrid<GridType>> dd_grid = nullptr)
  {
    static_assert(std::is_same<DdGridType, DD::SubdomainGrid<GridType>>::value,
                  "Only available for DD::SubdomainGrid!");
    return dd_grid->boundary_grid_view(subdomain);
  }

  static type create(GridType& /*grid*/, const int subdomain, std::shared_ptr<DD::SubdomainGrid<GridType>> dd_grid)
  {
    static_assert(std::is_same<DdGridType, DD::SubdomainGrid<GridType>>::value,
                  "Only available for DD::SubdomainGrid!");
    return dd_grid->boundary_grid_view(subdomain);
  }
}; // struct Layer<..., dd_subdomain_boundary, view>


// create periodic grid_layer
template <class GridType, Layers layer, Backends backend, class DdGridType>
struct Layer<GridType, layer, backend, DdGridType, true>
{
  using NonPeriodicLayerType = Layer<GridType, layer, backend, DdGridType, false>;
  using type = XT::Grid::PeriodicGridLayer<typename NonPeriodicLayerType::type>;

  static type create(const GridType& grid,
                     const int subdomain = 0,
                     const std::shared_ptr<DD::SubdomainGrid<GridType>> dd_grid = nullptr)
  {
    return type(NonPeriodicLayerType::create(grid, subdomain, dd_grid));
  }

  static type create(GridType& grid, const int subdomain, std::shared_ptr<DD::SubdomainGrid<GridType>> dd_grid)
  {
    return type(NonPeriodicLayerType::create(grid, subdomain, dd_grid));
  }
}; // struct Layer<..., periodic>


namespace internal {


template <class GL,
          bool view = is_view<GL>::value,
          bool part = is_part<GL>::value,
          bool dd_subdomain = is_dd_subdomain<GL>::value>
struct extract_layer_backend_helper
{
  static_assert(AlwaysFalse<GL>::value, "Please add a specialization for this case!");
};

template <class GL>
struct extract_layer_backend_helper<GL, true, false, false>
{
  static const constexpr Backends value = Backends::view;
};

template <class GL>
struct extract_layer_backend_helper<GL, false, false, true>
{
  static const constexpr Backends value = Backends::view;
};


} // namespace internal


template <class GridLayerType>
struct extract_layer_backend
{
  static const constexpr Backends value = internal::extract_layer_backend_helper<GridLayerType>::value;
};


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_LAYERS_HH
