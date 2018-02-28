// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   Rene Milk       (2014 - 2018)

#ifndef DUNE_XT_GRID_FUNCTORS_INTERFACES_HH
#define DUNE_XT_GRID_FUNCTORS_INTERFACES_HH

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/grid/type_traits.hh>


namespace Dune {
namespace XT {
namespace Grid {


/**
 * \brief Interface for functors which are applied to elements (codim 0 entities) of a grid layer by the \sa GridWalker.
 * \sa    ElementReturnFunctor
 */
template <class GL>
class ElementFunctor
{
  static_assert(is_layer<GL>::value, "");

protected:
  //! force implementors to use copy() method
  ElementFunctor(const ElementFunctor<GL>&) = default;

public:
  using GridViewType = GL;
  using ElementType = extract_entity_t<GridViewType>;

  using GV = GridViewType;
  using E = ElementType;

  ElementFunctor() = default;
  virtual ~ElementFunctor() = default;

  virtual void prepare()
  {
  }

  virtual void apply_local(const ElementType& element) = 0;

  virtual void finalize()
  {
  }

  virtual ElementFunctor<GridViewType>* copy() = 0;
}; // class ElementFunctor


/**
 * \brief Interface for functors which are applied to (codim 1) intersection of a grid layer by the \sa GridWalker.
 * \sa    ElementFunctor
 */
template <class GL>
class IntersectionFunctor
{
  static_assert(is_layer<GL>::value, "");

protected:
  //! force implementors to use copy() method
  IntersectionFunctor(const IntersectionFunctor<GL>&) = default;

public:
  using GridViewType = GL;
  using ElementType = extract_entity_t<GridViewType>;
  using IntersectionType = extract_intersection_t<GridViewType>;

  using GV = GridViewType;
  using E = ElementType;
  using I = IntersectionType;

  IntersectionFunctor() = default;
  virtual ~IntersectionFunctor() = default;

  virtual void prepare()
  {
  }

  /**
   * \note The meaning of outside_intersection depends on the circumstances. In general, the result of
   *       intersection.outside() is given, but this might differ on periodic or boundary intersections.
   */
  virtual void apply_local(const IntersectionType& intersection,
                           const ElementType& inside_element,
                           const ElementType& outside_element) = 0;

  virtual void finalize()
  {
  }

  virtual IntersectionFunctor<GridViewType>* copy() = 0;
}; // class IntersectionFunctor

/**
 * \brief Interface for functors which are applied to entities and intersections of a grid layer by the \sa GridWalker.
 * \sa    ElementFunctor
 * \sa    IntersectionFunctor
 */
template <class GL>
class ElementAndIntersectionFunctor
{
  static_assert(is_layer<GL>::value, "");

protected:
  //! force implementors to use copy() method
  ElementAndIntersectionFunctor(const ElementAndIntersectionFunctor<GL>&) = default;

public:
  using GridViewType = GL;
  using ElementType = extract_entity_t<GridViewType>;
  using IntersectionType = extract_intersection_t<GridViewType>;

  using GV = GridViewType;
  using E = ElementType;
  using I = IntersectionType;

  ElementAndIntersectionFunctor() = default;
  virtual ~ElementAndIntersectionFunctor() = default;

  virtual void prepare()
  {
  }

  virtual void apply_local(const ElementType& element) = 0;

  virtual void apply_local(const IntersectionType& intersection,
                           const ElementType& inside_element,
                           const ElementType& outside_element) = 0;

  virtual void finalize()
  {
  }

  virtual ElementAndIntersectionFunctor<GL>* copy() = 0;
}; // class ElementAndIntersectionFunctor


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_FUNCTORS_INTERFACES_HH
