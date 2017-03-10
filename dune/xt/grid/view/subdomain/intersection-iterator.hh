// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_GRID_VIEW_SUBDOMAIN_INTERSECTION_ITERATOR_HH
#define DUNE_XT_GRID_VIEW_SUBDOMAIN_INTERSECTION_ITERATOR_HH

#include <map>
#include <vector>
#include <algorithm>

#include <dune/common/shared_ptr.hh>

#include "intersection-wrapper.hh"

namespace Dune {
namespace XT {
namespace Grid {
namespace internal {


template <class GlobalGridPartImp>
class LocalIntersectionIterator : public GlobalGridPartImp::IntersectionIteratorType
{
public:
  typedef GlobalGridPartImp GlobalGridPartType;
  typedef LocalIntersectionIterator<GlobalGridPartType> ThisType;
  typedef typename GlobalGridPartType::IntersectionIteratorType BaseType;
  typedef typename BaseType::Intersection Intersection;
  typedef typename GlobalGridPartType::template Codim<0>::EntityType EntityType;
  typedef typename GlobalGridPartType::IndexSetType::IndexType IndexType;
  typedef std::vector<int> IndexContainerType;

  LocalIntersectionIterator(const GlobalGridPartType& globalGridPart,
                            const EntityType& entity,
                            const IndexContainerType& indexContainer,
                            const bool end = false)
    : BaseType(end ? globalGridPart.iend(entity) : globalGridPart.ibegin(entity))
    , globalGridPart_(globalGridPart)
    , entity_(entity)
    , indexContainer_(indexContainer)
    , workAtAll_(0)
  {
    if (!end) {
      if (indexContainer_.size() > 0) {
        last_ = *(indexContainer_.rbegin());
        end_ = indexContainer_.end();
        ++workAtAll_;
      }
      forward();
    } // if (!end)
  } // LocalIntersectionIterator

  ThisType& operator++()
  {
    if (workAtAll_ > 0) {
      BaseType::operator++();
      forward();
    } else
      BaseType::operator=(globalGridPart_.iend(entity_));
    return *this;
  } // ThisType& operator++()

private:
  //! iterates forward until we find the next intersection of interest
  void forward()
  {
    bool found = false;
    while (!found && (workAtAll_ > 0)) {
      const Intersection& intersection = BaseType::operator*();
      const int index = intersection.indexInInside();
      const auto result = std::find(indexContainer_.begin(), indexContainer_.end(), index);
      if (result != end_) {
        found = true;
        if (*result == last_)
          --workAtAll_;
      } else
        BaseType::operator++();
    } // while (!found && (workAtAll_ > 0))
  } // void forward()

  const GlobalGridPartType& globalGridPart_;
  const EntityType& entity_;
  const IndexContainerType& indexContainer_;
  unsigned int workAtAll_;
  int last_;
  typename IndexContainerType::const_iterator end_;
}; // class LocalIntersectionIterator


template <class GlobalGridPartImp>
class FakeDomainBoundaryIntersectionIterator : public GlobalGridPartImp::IntersectionIteratorType
{
public:
  typedef GlobalGridPartImp GlobalGridPartType;
  typedef FakeDomainBoundaryIntersectionIterator<GlobalGridPartType> ThisType;
  typedef typename GlobalGridPartType::IntersectionIteratorType BaseType;
  typedef typename GlobalGridPartType::template Codim<0>::EntityType EntityType;
  typedef std::map<int, size_t> InfoContainerType;

private:
  typedef typename BaseType::Intersection BaseIntersectionType;

public:
  typedef FakeDomainBoundaryIntersection<ThisType, BaseIntersectionType> Intersection;

  FakeDomainBoundaryIntersectionIterator(const GlobalGridPartType& globalGridPart,
                                         const EntityType& entity,
                                         bool end = false)
    : BaseType(end ? globalGridPart.iend(entity) : globalGridPart.ibegin(entity))
    , passThrough_(true)
    , intersection_(*this)
  {
  }

  FakeDomainBoundaryIntersectionIterator(const GlobalGridPartType& globalGridPart,
                                         const EntityType& entity,
                                         const InfoContainerType infoContainer,
                                         bool end = false)
    : BaseType(end ? globalGridPart.iend(entity) : globalGridPart.ibegin(entity))
    , passThrough_(false)
    , intersection_(*this)
    , infoContainer_(infoContainer)
  {
  }

  const Intersection& operator*() const
  {
    setIntersectionState();
    return intersection_;
  }

  const Intersection* operator->() const
  {
    setIntersectionState();
    return &intersection_;
  }

private:
  friend class FakeDomainBoundaryIntersection<ThisType, BaseIntersectionType>;

  const BaseIntersectionType& getBaseIntersection() const
  {
    return BaseType::operator*();
  }

  void setIntersectionState() const
  {
    // if we are on an entity of interest
    if (passThrough_) {
      intersection_.setPassThrough(true);
    } else {
      const int intersectionIndex = getBaseIntersection().indexInInside();
      // if this intersection is special
      typename InfoContainerType::const_iterator result = infoContainer_.find(intersectionIndex);
      if (result != infoContainer_.end()) {
        intersection_.setPassThrough(false);
        intersection_.setBoundarySegmentIndex(result->second);
      } else {
        intersection_.setPassThrough(true);
      } // if this intersection is special
    } // if we are not on an entity of interest
  } // void setIntersectionState() const

  bool passThrough_;
  mutable Intersection intersection_;
  const InfoContainerType infoContainer_;
}; // class FakeDomainBoundaryIntersectionIterator


} // namespace internal
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_VIEW_SUBDOMAIN_INTERSECTION_ITERATOR_HH
