// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2019 - 2020)

#ifndef DUNE_XT_GRID_FILTERS_ELEMENT_HH
#define DUNE_XT_GRID_FILTERS_ELEMENT_HH

#include <functional>

#include <dune/xt/common/memory.hh>

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/grid/common/partitionset.hh>

#include "base.hh"

namespace Dune::XT::Grid::ApplyOn {


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
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;

  explicit AllElements() = default;

  ElementFilter<GridViewType>* copy() const final
  {
    return new AllElements<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const ElementType& /*element*/) const final
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
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;

  explicit NoElements() = default;

  ElementFilter<GridViewType>* copy() const final
  {
    return new NoElements<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const ElementType& /*element*/) const final
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
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;

  explicit BoundaryElements() = default;

  ElementFilter<GridViewType>* copy() const final
  {
    return new BoundaryElements<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const ElementType& element) const final
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
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;
  using GenericFunctionType = std::function<bool(const GridViewType&, const ElementType&)>;

  explicit GenericFilteredElements(GenericFunctionType lambda)
    : filter_(lambda)
  {}

  ElementFilter<GridViewType>* copy() const final
  {
    return new GenericFilteredElements<GridViewType>(filter_);
  }

  bool contains(const GridViewType& grid_layer, const ElementType& element) const final
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
  using ElementType = typename BaseType::ElementType;

  explicit PartitionSetElements() = default;

  ElementFilter<GridViewType>* copy() const final
  {
    return new PartitionSetElements<GridViewType, PartitionSetType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const ElementType& element) const final
  {
    return PartitionSetType::contains(element.partitionType());
  }
}; // class PartitionSetElements


} // namespace Dune::XT::Grid::ApplyOn

#endif // DUNE_XT_GRID_FILTERS_ELEMENT_HH
