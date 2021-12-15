// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018 - 2020)
//   Tobias Leibner  (2019 - 2020)

#ifndef DUNE_XT_GRID_BOUND_OBJECTS_HH
#define DUNE_XT_GRID_BOUND_OBJECTS_HH

#include <dune/xt/common/timedlogging.hh>

#include "entity.hh"
#include "exceptions.hh"
#include "type_traits.hh"

namespace Dune::XT::Grid {


template <class Element>
class ElementBoundObject /*: Common::EnableDebugLoggingForCtors<ElementBoundObject<Element>>*/
{
  static_assert(is_entity<Element>::value);

  using ThisType = ElementBoundObject;
  //  using DebugLogging = Common::EnableDebugLoggingForCtors<ElementBoundObject<Element>>;

public:
  using ElementType = Element;

  ElementBoundObject()
    : /*DebugLogging("dune.xt.grid", "ElementBoundObject")
    ,*/
    element_(nullptr)

  {}

  ElementBoundObject(const ThisType& other)
    : /*DebugLogging(other)
    ,*/
    element_(nullptr)

  {
    if (other.element_) {
      element_ = std::make_unique<ElementType>(*other.element_);
      is_bound_ = true;
    }
  }

  ElementBoundObject(ThisType&& source) noexcept
    : /*DebugLogging(source)
,*/
    element_(std::move(source.element_))
    , is_bound_(source.is_bound_)
  {}

  virtual ~ElementBoundObject() = default;

  ThisType& operator=(const ThisType& other)
  {
    //    DebugLogging::operator=(other);
    if (&other != this && other.element_) {
      element_ = std::make_unique<ElementType>(other.element_);
      is_bound_ = true;
    }
  }

  ThisType& operator=(ThisType&& source) noexcept
  {
    //    DebugLogging::operator=(source);
    if (&source != this) {
      element_ = source.element_;
      is_bound_ = source.is_bound_;
    }
  }

  /**
   * \attention The returned reference will change as soon as the object is bound to another element!
   * \note      Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  const ElementType& element() const
  {
    DUNE_THROW_IF(!is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    return *element_;
  }

  ThisType& bind(const ElementType& ele)
  {
    if (element_ && ele == *element_)
      return *this;
    element_ = std::make_unique<ElementType>(ele);
    is_bound_ = true;
    this->post_bind(*element_);
    return *this;
  }

protected:
  /**
   * \note Override this function if you need/want to do preparatory work on an element.
   */
  virtual void post_bind(const ElementType& /*ele*/) {}

private:
  std::unique_ptr<ElementType> element_;

protected:
  bool is_bound_{false};
}; // class ElementBoundObject


template <class Intersection>
class IntersectionBoundObject
{
  static_assert(is_intersection<Intersection>::value);

  using ThisType = IntersectionBoundObject;

public:
  using IntersectionType = Intersection;

  IntersectionBoundObject()
    : intersection_(nullptr)

  {}

  IntersectionBoundObject(const ThisType& other)
    : intersection_(nullptr)

  {
    if (other.intersection_) {
      intersection_ = std::make_unique<IntersectionType>(*other.intersection_);
      is_bound_ = true;
    }
  }

  IntersectionBoundObject(ThisType&& source) = default;

  virtual ~IntersectionBoundObject() = default;

  ThisType& operator=(const ThisType& other)
  {
    if (&other != this && other.intersection_) {
      intersection_ = std::make_unique<IntersectionType>(other.intersection_);
      is_bound_ = true;
    }
  }

  ThisType& operator=(ThisType&& source) noexcept
  {
    if (&source != this) {
      intersection_ = source.intersection_;
      is_bound_ = source.is_bound_;
    }
  }

  /**
   * \attention The returned reference will change as soon as the object is bound to another Intersection!
   * \note      Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  const IntersectionType& intersection() const
  {
    DUNE_THROW_IF(!is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    return *intersection_;
  }

  ThisType& bind(const IntersectionType& inter)
  {
    //    if (intersection_ && inter == *intersection_)
    //      return *this;
    intersection_ = std::make_unique<IntersectionType>(inter);
    is_bound_ = true;
    this->post_bind(*intersection_);
    return *this;
  }

protected:
  /**
   * \note Override this function if you need/want to do preparatory work on an intersection.
   */
  virtual void post_bind(const IntersectionType& /*inter*/) {}

private:
  std::unique_ptr<IntersectionType> intersection_;

protected:
  bool is_bound_{false};
}; // class IntersectionBoundObject


} // namespace Dune::XT::Grid

#endif // DUNE_XT_GRID_BOUND_OBJECTS_HH
