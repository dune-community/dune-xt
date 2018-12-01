// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2018)

#ifndef DUNE_XT_GRID_VIEW_SUBDOMAIN_INTERSECTION_WRAPPER_HH
#define DUNE_XT_GRID_VIEW_SUBDOMAIN_INTERSECTION_WRAPPER_HH

#include <boost/numeric/conversion/cast.hpp>

#include <dune/geometry/type.hh>

#include <dune/grid/common/intersection.hh>

#include <dune/xt/common/exceptions.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace internal {


template <class IntersectionIteratorImp, class WrappedIntersectionImp>
class FakeDomainBoundaryIntersection
{
public:
  typedef IntersectionIteratorImp IntersectionIteratorType;
  typedef WrappedIntersectionImp WrappedIntersectionType;
  typedef FakeDomainBoundaryIntersection<IntersectionIteratorType, WrappedIntersectionType> ThisType;
  typedef typename WrappedIntersectionType::Entity Entity;
  typedef typename WrappedIntersectionType::Geometry Geometry;
  typedef typename WrappedIntersectionType::LocalCoordinate LocalCoordinate;
  typedef typename WrappedIntersectionType::GlobalCoordinate GlobalCoordinate;
  typedef typename WrappedIntersectionType::LocalGeometry LocalGeometry;
  typedef typename WrappedIntersectionType::ctype ctype;
  static const int dimension = Entity::dimension;
  static const int dimensionworld = WrappedIntersectionType::dimensionworld;

  FakeDomainBoundaryIntersection(const IntersectionIteratorType& intersectionIterator)
    : intersectionIterator_(intersectionIterator)
    , passThrough_(true)
    , boundary_segment_index_(std::numeric_limits<size_t>::max())
  {}

  void setPassThrough(const bool passThrough)
  {
    passThrough_ = passThrough;
  }

  void setBoundarySegmentIndex(const size_t index)
  {
    boundary_segment_index_ = index;
  }

  bool neighbor() const
  {
    if (passThrough_)
      return intersectionIterator_.getBaseIntersection().neighbor();
    else
      return false;
  }

  bool boundary() const
  {
    if (passThrough_)
      return intersectionIterator_.getBaseIntersection().boundary();
    else
      return true;
  }

  size_t boundarySegmentIndex() const
  {
    if (passThrough_)
      return intersectionIterator_.getBaseIntersection().boundarySegmentIndex();
    else
      return boundary_segment_index_;
  }

  Entity inside() const
  {
    return intersectionIterator_.getBaseIntersection().inside();
  }

  Entity outside() const
  {
    return intersectionIterator_.getBaseIntersection().outside();
  }

  bool conforming() const
  {
    return intersectionIterator_.getBaseIntersection().conforming();
  }

  LocalGeometry geometryInInside() const
  {
    return intersectionIterator_.getBaseIntersection().geometryInInside();
  }

  LocalGeometry geometryInOutside() const
  {
    return intersectionIterator_.getBaseIntersection().geometryInOutside();
  }

  Geometry geometry() const
  {
    return intersectionIterator_.getBaseIntersection().geometry();
  }

  Dune::GeometryType type() const
  {
    return intersectionIterator_.getBaseIntersection().type();
  }

  int indexInInside() const
  {
    return intersectionIterator_.getBaseIntersection().indexInInside();
  }

  int indexInOutside() const
  {
    return intersectionIterator_.getBaseIntersection().indexInOutside();
  }

  GlobalCoordinate outerNormal(const LocalCoordinate& local) const
  {
    return intersectionIterator_.getBaseIntersection().outerNormal(local);
  }

  GlobalCoordinate integrationOuterNormal(const LocalCoordinate& local) const
  {
    return intersectionIterator_.getBaseIntersection().integrationOuterNormal(local);
  }

  GlobalCoordinate unitOuterNormal(const LocalCoordinate& local) const
  {
    return intersectionIterator_.getBaseIntersection().unitOuterNormal(local);
  }

  GlobalCoordinate centerUnitOuterNormal() const
  {
    return intersectionIterator_.getBaseIntersection().centerUnitOuterNormal();
  }

  const typename IntersectionIteratorType::BaseIntersectionType& asBase() const
  {
    return intersectionIterator_.getBaseIntersection();
  }

private:
  const IntersectionIteratorType& intersectionIterator_;
  bool passThrough_;
  size_t boundary_segment_index_;
}; // class FakeDomainBoundaryIntersection


} // namespace internal
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_VIEW_SUBDOMAIN_INTERSECTION_WRAPPER_HH
