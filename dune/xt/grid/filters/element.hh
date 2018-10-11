// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   Rene Milk       (2014 - 2018)
//   Tobias Leibner  (2015 - 2017)

#ifndef DUNE_XT_GRID_FILTERS_ELEMENT_HH
#define DUNE_XT_GRID_FILTERS_ELEMENT_HH

#include <functional>

#include <dune/xt/common/memory.hh>

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/grid/common/partitionset.hh>

#include "base.hh"

namespace Dune {
namespace XT {
namespace Grid {
namespace ApplyOn {


/**
 * \brief A filter which selects all elements.
 *
 * \sa Walker
 * \sa ElementFilter
 */
template <class GL>
class AllElements : public ElementFilter<GL>
{
  using BaseType = ElementFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;

  explicit AllElements() = default;

  ElementFilter<GridViewType>* copy() const override final
  {
    return new AllElements<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const ElementType& /*element*/) const override final
  {
    return true;
  }
}; // class AllElements


/**
 * \brief A filter which selects no elements.
 *
 * \sa Walker
 * \sa ElementFilter
 */
template <class GL>
class NoElements : public ElementFilter<GL>
{
  using BaseType = ElementFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;

  explicit NoElements() = default;

  ElementFilter<GridViewType>* copy() const override final
  {
    return new NoElements<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const ElementType& /*element*/) const override final
  {
    return false;
  }
}; // class NoElements


/**
 * \brief A filter which selects entities which have a boundary intersection.
 *
 * \sa Walker
 * \sa ElementFilter
 */
template <class GL>
class BoundaryElements : public ElementFilter<GL>
{
  using BaseType = ElementFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;

  explicit BoundaryElements() = default;

  ElementFilter<GridViewType>* copy() const override final
  {
    return new BoundaryElements<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const ElementType& element) const override final
  {
    return element.hasBoundaryIntersections();
  }
}; // class BoundaryElements


/**
 * \brief A filter which selects entities based on a lambda expression.
 *
 * \sa Walker
 * \sa ElementFilter
 */
template <class GL>
class GenericFilteredElements : public ElementFilter<GL>
{
  using BaseType = ElementFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;
  using GenericFunctionType = std::function<bool(const GridViewType&, const ElementType&)>;

  explicit GenericFilteredElements(GenericFunctionType lambda)
    : filter_(lambda)
  {
  }

  ElementFilter<GridViewType>* copy() const override final
  {
    return new GenericFilteredElements<GridViewType>(filter_);
  }

  bool contains(const GridViewType& grid_layer, const ElementType& element) const override final
  {
    return filter_(grid_layer, element);
  }

private:
  const GenericFunctionType filter_;
}; // class GenericFilteredElements


/**
 * \brief A filter which selects entities in a compatible PartitionSet.
 *
 * \sa Walker
 * \sa ElementFilter
 */
template <class GL, class PartitionSetType>
class PartitionSetElements : public ElementFilter<GL>
{
  using BaseType = ElementFilter<GL>;

public:
  using typename BaseType::GridViewType;
  typedef typename BaseType::ElementType ElementType;

  explicit PartitionSetElements() = default;

  ElementFilter<GridViewType>* copy() const override final
  {
    return new PartitionSetElements<GridViewType, PartitionSetType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const ElementType& element) const override final
  {
    return PartitionSetType::contains(element.partitionType());
  }
}; // class PartitionSetElements


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_FILTERS_ELEMENT_HH
