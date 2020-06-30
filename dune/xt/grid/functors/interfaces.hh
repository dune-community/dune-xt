// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2018)
//   René Fritze     (2014 - 2019)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_GRID_FUNCTORS_INTERFACES_HH
#define DUNE_XT_GRID_FUNCTORS_INTERFACES_HH

#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/grid/type_traits.hh>


namespace Dune {
namespace XT {
namespace Grid {


template <class GL>
class ElementFunctor;

/**
 * \brief Interface for functors which are applied to elements (codim 0 entities) of a grid layer by the Walker.
 *
 * \sa Walker
 * \sa IntersectionFunctor
 * \sa ElementAndIntersectionFunctor
 */
template <class GL>
class ElementFunctor : public Common::WithLogger<ElementFunctor<GL>>
{
  static_assert(is_layer<GL>::value, "");

protected:
  //! force implementors to use copy() method

public:
  using GridViewType = GL;
  using ElementType = extract_entity_t<GridViewType>;

  using GV = GridViewType;
  using E = ElementType;

  ElementFunctor(const std::string& log_prefix = "", const std::string& log_id = "", const bool log_disabled = true)
    : Common::WithLogger<ElementFunctor<GL>>(
          log_prefix.empty() ? "xt.grid" : log_prefix, log_id.empty() ? "ElementFunctor" : log_id, log_disabled)
  {}

  ElementFunctor(const ElementFunctor<GL>&) = default;

  virtual ~ElementFunctor() = default;

  virtual ElementFunctor<GridViewType>* copy() = 0;

  virtual void prepare() {}

  virtual void apply_local(const ElementType& /*element*/) {}

  virtual void finalize() {}
}; // class ElementFunctor


template <class GL>
class IntersectionFunctor;

/**
 * \brief Interface for functors which are applied to intersections (codim 1 entities) of a grid layer by the Walker.
 *
 * \sa Walker
 * \sa ElementFunctor
 * \sa ElementAndIntersectionFunctor
 */
template <class GL>
class IntersectionFunctor : public Common::WithLogger<IntersectionFunctor<GL>>
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

  IntersectionFunctor(const std::string& log_prefix = "",
                      const std::string& log_id = "",
                      const bool log_disabled = true)
    : Common::WithLogger<IntersectionFunctor<GL>>(
          log_prefix.empty() ? "xt.grid" : log_prefix, log_id.empty() ? "IntersectionFunctor" : log_id, log_disabled)
  {}

  virtual ~IntersectionFunctor() = default;

  virtual IntersectionFunctor<GridViewType>* copy() = 0;

  virtual void prepare() {}

  /**
   * \note The meaning of outside_element depends on the circumstances. If intersection.neighbor() is true, the result
   *       of intersection.outside() is given (the meaning of which is different on inner, periodic or process boundary
   *       intersections). If intersection.neighbor() is false, intersection.inside() is given.
   */
  virtual void apply_local(const IntersectionType& /*intersection*/,
                           const ElementType& /*inside_element*/,
                           const ElementType& /*outside_element*/)
  {}

  virtual void finalize() {}
}; // class IntersectionFunctor


template <class GL>
class ElementAndIntersectionFunctor;

/**
 * \brief Interface for functors which are applied to entities and intersections of a grid layer by the Walker.
 *
 * \sa Walker
 * \sa ElementFunctor
 * \sa IntersectionFunctor
 */
template <class GL>
class ElementAndIntersectionFunctor : public Common::WithLogger<ElementAndIntersectionFunctor<GL>>
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

  ElementAndIntersectionFunctor(const std::string& log_prefix = "",
                                const std::string& log_id = "",
                                const bool log_disabled = true)
    : Common::WithLogger<ElementAndIntersectionFunctor<GL>>(log_prefix.empty() ? "xt.grid" : log_prefix,
                                                            log_id.empty() ? "ElementAndIntersectionFunctor" : log_id,
                                                            log_disabled)
  {}

  virtual ~ElementAndIntersectionFunctor() = default;

  virtual ElementAndIntersectionFunctor<GL>* copy() = 0;

  virtual void prepare() {}

  virtual void apply_local(const ElementType& /*element*/) {}

  virtual void apply_local(const IntersectionType& /*intersection*/,
                           const ElementType& /*inside_element*/,
                           const ElementType& /*outside_element*/)
  {}

  virtual void finalize() {}
}; // class ElementAndIntersectionFunctor


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_FUNCTORS_INTERFACES_HH
