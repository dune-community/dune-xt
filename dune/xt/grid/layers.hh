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
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_GRID_LAYERS_HH
#define DUNE_XT_GRID_LAYERS_HH

#include <dune/common/typetraits.hh>

#include <dune/xt/grid/type_traits.hh>
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


namespace internal {


template <Backends backend>
struct backend_dependent_typename
{
  typedef void type;
};


template <Layers layer>
struct layer_dependent_typename
{
  typedef void type;
};


} // namespace  internal


/**
 * \brief Allows to statically create a leaf or level part or view (unspecialized variant).
 */
template <class GridType, Layers layer, Backends backend, class DdGridType = int>
struct Layer
{
  static_assert(AlwaysFalse<GridType>::value, "No layer available for this combination!");

  typedef void type;

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
struct Layer<GridType, Layers::leaf, Backends::view, DdGridType>
{
  typedef typename GridType::LeafGridView type;

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
struct Layer<GridType, Layers::level, Backends::view, DdGridType>
{
  typedef typename GridType::LevelGridView type;

  static type create(const GridType& grid, const int level, const std::shared_ptr<DdGridType> /*dd_grid*/ = nullptr)
  {
    assert(level >= 0);
    assert(level <= grid.maxLevel());
    return grid.levelGridView(level);
  }
}; // struct Layer<..., level, view>


template <class GridType, class DdGridType>
struct Layer<GridType, Layers::dd_subdomain, Backends::view, DdGridType>
{
  typedef typename DD::SubdomainGrid<GridType>::LocalGridViewType type;

  static type create(const GridType& /*grid*/,
                     const int /*subdomain*/ = 0,
                     const std::shared_ptr<DD::SubdomainGrid<GridType>> /*dd_grid*/ = nullptr)
  {
    static_assert(AlwaysFalse<GridType>::value,
                  "dune-fem does not allow the creation of grid parts from a const grid!");
  }

  static type create(GridType& /*grid*/, const int subdomain, std::shared_ptr<DD::SubdomainGrid<GridType>> dd_grid)
  {
    static_assert(std::is_same<DdGridType, DD::SubdomainGrid<GridType>>::value,
                  "Only available for DD::SubdomainGrid!");
    return dd_grid->local_grid_view(subdomain, /*oversampling=*/false);
  }
}; // struct Layer<..., dd_subdomain, view>


template <class GridType, class DdGridType>
struct Layer<GridType, Layers::dd_subdomain_oversampled, Backends::view, DdGridType>
{
  typedef typename DD::SubdomainGrid<GridType>::LocalGridViewType type;

  static type create(const GridType& /*grid*/,
                     const int /*subdomain*/ = 0,
                     const std::shared_ptr<DD::SubdomainGrid<GridType>> /*dd_grid*/ = nullptr)
  {
    static_assert(AlwaysFalse<GridType>::value,
                  "dune-fem does not allow the creation of grid parts from a const grid!");
  }

  static type create(GridType& /*grid*/, const int subdomain, std::shared_ptr<DD::SubdomainGrid<GridType>> dd_grid)
  {
    static_assert(std::is_same<DdGridType, DD::SubdomainGrid<GridType>>::value,
                  "Only available for DD::SubdomainGrid!");
    return dd_grid->local_grid_view(subdomain, /*oversampling=*/true);
  }
}; // struct Layer<..., dd_subdomain_oversampled, view>


template <class GridType, class DdGridType>
struct Layer<GridType, Layers::dd_subdomain_coupling, Backends::view, DdGridType>
{
  typedef typename DD::SubdomainGrid<GridType>::CouplingGridViewType type;

  static type create(const GridType& /*grid*/,
                     const int /*subdomain*/ = 0,
                     const std::shared_ptr<DD::SubdomainGrid<GridType>> /*dd_grid*/ = nullptr)
  {
    static_assert(AlwaysFalse<GridType>::value,
                  "dune-fem does not allow the creation of grid parts from a const grid!");
  }

  static type create(GridType& /*grid*/, const int subdomain, std::shared_ptr<DD::SubdomainGrid<GridType>> dd_grid)
  {
    static_assert(std::is_same<DdGridType, DD::SubdomainGrid<GridType>>::value,
                  "Only available for DD::SubdomainGrid!");
    DUNE_THROW(NotImplemented, "Only usable to extract the layer type, not the actual layer!");
    return dd_grid->coupling_grid_view(0, 0);
  }
}; // struct Layer<..., dd_subdomain_coupling, view>


template <class GridType, class DdGridType>
struct Layer<GridType, Layers::dd_subdomain_boundary, Backends::view, DdGridType>
{
  typedef typename DD::SubdomainGrid<GridType>::BoundaryGridViewType type;

  static type create(const GridType& /*grid*/,
                     const int /*subdomain*/ = 0,
                     const std::shared_ptr<DD::SubdomainGrid<GridType>> /*dd_grid*/ = nullptr)
  {
    static_assert(AlwaysFalse<GridType>::value,
                  "dune-fem does not allow the creation of grid parts from a const grid!");
  }

  static type create(GridType& /*grid*/, const int subdomain, std::shared_ptr<DD::SubdomainGrid<GridType>> dd_grid)
  {
    static_assert(std::is_same<DdGridType, DD::SubdomainGrid<GridType>>::value,
                  "Only available for DD::SubdomainGrid!");
    return dd_grid->boundary_grid_view(subdomain);
  }
}; // struct Layer<..., dd_subdomain_boundary, view>


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
