// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GRID_PART_ITERATOR_INTERSECTION_LOCAL_HH
#define DUNE_GRID_PART_ITERATOR_INTERSECTION_LOCAL_HH

#include <map>
#include <vector>
#include <algorithm>

#include <dune/common/shared_ptr.hh>

namespace Dune {

namespace grid {

namespace Part {

namespace Iterator {

namespace Intersection {

template <class GlobalGridPartImp>
class Local : public GlobalGridPartImp::IntersectionIteratorType
{
public:
  typedef GlobalGridPartImp GlobalGridPartType;

  typedef Local<GlobalGridPartType> ThisType;

  typedef typename GlobalGridPartType::IntersectionIteratorType BaseType;

  typedef typename BaseType::Intersection Intersection;

  typedef typename GlobalGridPartType::template Codim<0>::EntityType EntityType;

  typedef typename GlobalGridPartType::IndexSetType::IndexType IndexType;

  typedef std::vector<int> IndexContainerType;

  Local(const GlobalGridPartType& globalGridPart,
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
  } // Local

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
}; // class Local

namespace Wrapper {

template <class GlobalGridPartImp>
class FakeDomainBoundary : public GlobalGridPartImp::IntersectionIteratorType
{
public:
  typedef GlobalGridPartImp GlobalGridPartType;

  typedef FakeDomainBoundary<GlobalGridPartType> ThisType;

  typedef typename GlobalGridPartType::IntersectionIteratorType BaseType;

  typedef typename GlobalGridPartType::template Codim<0>::EntityType EntityType;

  typedef std::map<int, size_t> InfoContainerType;

private:
  typedef typename BaseType::Intersection BaseIntersectionType;

public:
  typedef Dune::grid::Part::Intersection::Wrapper::FakeDomainBoundary<ThisType, BaseIntersectionType> Intersection;

  FakeDomainBoundary(const GlobalGridPartType& globalGridPart, const EntityType& entity, bool end = false)
    : BaseType(end ? globalGridPart.iend(entity) : globalGridPart.ibegin(entity))
    , passThrough_(true)
    , intersection_(*this)
  {
  }

  FakeDomainBoundary(const GlobalGridPartType& globalGridPart,
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
  friend class Dune::grid::Part::Intersection::Wrapper::FakeDomainBoundary<ThisType, BaseIntersectionType>;

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
}; // class FakeDomainBoundary

} // namespace Wrapper

} // namespace Intersection

} // namespace Iterator

} // namespace Part

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_PART_ITERATOR_INTERSECTION_LOCAL_HH
