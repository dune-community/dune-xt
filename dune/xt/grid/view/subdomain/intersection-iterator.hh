// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018)

#ifndef DUNE_XT_GRID_VIEW_SUBDOMAIN_INTERSECTION_ITERATOR_HH
#define DUNE_XT_GRID_VIEW_SUBDOMAIN_INTERSECTION_ITERATOR_HH

#include <map>
#include <vector>
#include <algorithm>

#include <dune/common/shared_ptr.hh>

#include <dune/xt/grid/type_traits.hh>

#include "intersection-wrapper.hh"

namespace Dune {
namespace XT {
namespace Grid {
namespace internal {


template <class GlobalGridViewImp>
class LocalIntersectionIterator : public GlobalGridViewImp::IntersectionIterator
{
public:
  typedef GlobalGridViewImp GlobalGridViewType;
  typedef LocalIntersectionIterator<GlobalGridViewType> ThisType;
  typedef typename GlobalGridViewType::IntersectionIterator BaseType;
  typedef typename BaseType::Intersection Intersection;
  using EntityType = extract_entity_t<GlobalGridViewType>;
  typedef typename GlobalGridViewType::IndexSet::IndexType IndexType;
  typedef std::vector<int> IndexContainerType;

  LocalIntersectionIterator(const GlobalGridViewType& GlobalGridView,
                            const EntityType& entity,
                            const IndexContainerType& indexContainer,
                            const bool end = false)
    : BaseType(end ? GlobalGridView.iend(entity) : GlobalGridView.ibegin(entity))
    , globalGridView_(GlobalGridView)
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
      BaseType::operator=(globalGridView_.iend(entity_));
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

  const GlobalGridViewType& globalGridView_;
  const EntityType& entity_;
  const IndexContainerType& indexContainer_;
  unsigned int workAtAll_;
  int last_;
  typename IndexContainerType::const_iterator end_;
}; // class LocalIntersectionIterator


template <class GlobalGridViewImp>
class FakeDomainBoundaryIntersectionIterator : public GlobalGridViewImp::IntersectionIterator
{
public:
  typedef GlobalGridViewImp GlobalGridViewType;
  typedef FakeDomainBoundaryIntersectionIterator<GlobalGridViewType> ThisType;
  typedef typename GlobalGridViewType::IntersectionIterator BaseType;
  using EntityType = extract_entity_t<GlobalGridViewType>;
  typedef std::map<int, size_t> InfoContainerType;

private:
  typedef typename BaseType::Intersection BaseIntersectionType;

public:
  typedef FakeDomainBoundaryIntersection<ThisType, BaseIntersectionType> IntersectionImp;
  typedef Dune::Intersection<typename GlobalGridViewImp::Grid, IntersectionImp> Intersection;

  FakeDomainBoundaryIntersectionIterator(const GlobalGridViewType& GlobalGridView,
                                         const EntityType& entity,
                                         bool end = false)
    : BaseType(end ? GlobalGridView.iend(entity) : GlobalGridView.ibegin(entity))
    , passThrough_(true)
    , intersection_(IntersectionImp(*this))
  {}

  FakeDomainBoundaryIntersectionIterator(const GlobalGridViewType& GlobalGridView,
                                         const EntityType& entity,
                                         const InfoContainerType infoContainer,
                                         bool end = false)
    : BaseType(end ? GlobalGridView.iend(entity) : GlobalGridView.ibegin(entity))
    , passThrough_(false)
    , intersection_(IntersectionImp(*this))
    , infoContainer_(infoContainer)
  {}

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
      intersection_.impl().setPassThrough(true);
    } else {
      const int intersectionIndex = getBaseIntersection().indexInInside();
      // if this intersection is special
      typename InfoContainerType::const_iterator result = infoContainer_.find(intersectionIndex);
      if (result != infoContainer_.end()) {
        intersection_.impl().setPassThrough(false);
        intersection_.impl().setBoundarySegmentIndex(result->second);
      } else {
        intersection_.impl().setPassThrough(true);
      } // if this intersection is special
    } // if we are not on an entity of interest
  } // void setIntersectionState() const

  bool passThrough_;
  mutable Intersection intersection_;
  const InfoContainerType infoContainer_;
}; // class FakeDomainBoundaryIntersectionIterator


} // namespace internal

/** this is only necessary because iterators get fubared in copying into the default iterator range,
 *  same goes for the intersections(...) functions at the bottom of xt/grid/view/subdomain/view.hh
 *
 */
template <class SubdomainGridViewImp>
class SubdomainGridViewIntersectionRange
{
  typedef typename SubdomainGridViewImp::EntityType EntityType;
  typedef typename SubdomainGridViewImp::IntersectionIterator IntersectionIterator;

public:
  SubdomainGridViewIntersectionRange(const SubdomainGridViewImp& layer, const EntityType& entity)
    : layer_(layer)
    , entity_(entity)
  {}

  IntersectionIterator begin() const
  {
    return layer_.ibegin(entity_);
  }

  IntersectionIterator end() const
  {
    return layer_.iend(entity_);
  }

  const SubdomainGridViewImp& layer_;
  const EntityType& entity_;
}; // class SubdomainGridViewIntersectionRange

} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_VIEW_SUBDOMAIN_INTERSECTION_ITERATOR_HH
