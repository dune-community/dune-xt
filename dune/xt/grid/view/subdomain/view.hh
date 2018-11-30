// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_XT_GRID_VIEW_SUBDOMAIN_PART_HH
#define DUNE_XT_GRID_VIEW_SUBDOMAIN_PART_HH

#include <map>
#include <set>
#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridview.hh>


#include "entity-iterator.hh"
#include "indexset.hh"
#include "intersection-iterator.hh"
#include "intersection-wrapper.hh"

namespace Dune {
namespace XT {
namespace Grid {


// forwards, needed for the traits below
template <class GlobalGridViewImp>
class SubdomainGridView;

template <class GlobalGridViewImp>
class SubdomainCouplingGridView;

template <class GlobalGridViewImp>
class SubdomainBoundaryGridView;


namespace internal {


template <class GlobalGridViewImp>
class SubdomainGridViewTraits
{
  static_assert(!is_part<GlobalGridViewImp>::value, "only works for real views");

public:
  typedef GlobalGridViewImp GlobalGridViewType;

  typedef SubdomainGridView<GlobalGridViewImp> GridViewType;
  using GridViewImp = GridViewType;
  typedef typename GlobalGridViewType::Grid Grid;
  typedef typename internal::IndexBasedIndexSet<GlobalGridViewType> IndexSet;
  typedef typename GlobalGridViewType::CollectiveCommunication CollectiveCommunication;

  //  static const PartitionIteratorType indexSetPartitionType = GlobalGridViewType::indexSetPartitionType;
  //  static const InterfaceType indexSetInterfaceType = GlobalGridViewType::indexSetInterfaceType;

  typedef internal::FakeDomainBoundaryIntersectionIterator<GlobalGridViewType> IntersectionIterator;
  using Intersection = typename IntersectionIterator::Intersection;

  template <int codim>
  struct Codim : public GlobalGridViewType::template Codim<codim>
  {
    template <PartitionIteratorType pitype>
    struct Partition
    {
      typedef typename internal::IndexBasedEntityIterator<GlobalGridViewType, codim, pitype> IteratorType;
      using Iterator = IteratorType;
    };

    using Iterator = typename Partition<All_Partition>::Iterator;
  };

  static const bool conforming = GlobalGridViewType::Traits::conforming;
}; // class SubdomainGridViewTraits


template <class GlobalGridViewImp>
struct SubdomainCouplingGridViewTraits : public SubdomainGridViewTraits<GlobalGridViewImp>
{
  typedef GlobalGridViewImp GlobalGridViewType;
  typedef SubdomainCouplingGridView<GlobalGridViewImp> GridViewType;
  using GridViewImp = GridViewType;
  typedef internal::LocalIntersectionIterator<GlobalGridViewType> IntersectionIterator;
  using Intersection = typename IntersectionIterator::Intersection;
};


template <class GlobalGridViewImp>
struct SubdomainBoundaryGridViewTraits : public SubdomainGridViewTraits<GlobalGridViewImp>
{
  typedef GlobalGridViewImp GlobalGridViewType;
  typedef SubdomainBoundaryGridView<GlobalGridViewImp> GridViewType;
  using GridViewImp = GridViewType;
  typedef internal::LocalIntersectionIterator<GlobalGridViewType> IntersectionIterator;
  using Intersection = typename IntersectionIterator::Intersection;
}; // class SubdomainBoundaryGridViewTraits


template <class Traits>
class SubdomainGridViewCommon
{
  using ThisType = SubdomainGridViewCommon<Traits>;
  using BaseType = GridView<Traits>;

public:
  typedef typename Traits::Grid Grid;
  typedef typename Traits::CollectiveCommunication CollectiveCommunication;
  typedef typename Traits::GlobalGridViewType GlobalGridViewType;
  typedef typename Traits::IndexSet IndexSet;
  typedef typename Traits::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;
  typedef typename Grid::template Codim<0>::Entity EntityType;

  typedef typename IndexSet::IndexType IndexType;
  typedef std::map<IndexType, IndexType> IndexMapType;
  typedef Dune::GeometryType GeometryType;
  //! container type for the indices
  typedef std::map<GeometryType, IndexMapType> IndexContainerType;
  //! container type for the boundary information
  typedef std::map<IndexType, std::map<int, size_t>> BoundaryInfoContainerType;

public:
  SubdomainGridViewCommon(const std::shared_ptr<const GlobalGridViewType> globalGrdPrt,
                          const std::shared_ptr<const IndexContainerType> indexContainer,
                          const std::shared_ptr<const BoundaryInfoContainerType> boundaryInfoContainer)
    : globalGridView_(globalGrdPrt)
    , indexContainer_(indexContainer)
    , boundaryInfoContainer_(boundaryInfoContainer)
    , indexSet_(std::make_shared<IndexSet>(*globalGridView_, indexContainer_))
  {
  }

  SubdomainGridViewCommon(const ThisType& other) = default;
  SubdomainGridViewCommon(ThisType&& source) = default;

  const IndexSet& indexSet() const
  {
    return *indexSet_;
  }

  const Grid& grid() const
  {
    return globalGridView_->grid();
  }

  const GlobalGridViewType& globalGridView() const
  {
    return *globalGridView_;
  }

  template <int codim>
  typename Traits::template Codim<codim>::Iterator begin() const
  {
    return typename Traits::template Codim<codim>::Iterator(*globalGridView_, indexContainer_);
  }

  template <int codim, PartitionIteratorType pitype>
  typename Traits::template Codim<codim>::template Partition<pitype>::Iterator begin() const
  {
    return
        typename Traits::template Codim<codim>::template Partition<pitype>::Iterator(*globalGridView_, indexContainer_);
  }

  template <int codim>
  typename Traits::template Codim<codim>::Iterator end() const
  {
    return typename Traits::template Codim<codim>::Iterator(*globalGridView_, indexContainer_, true);
  }

  template <int codim, PartitionIteratorType pitype>
  typename Traits::template Codim<codim>::template Partition<pitype>::Iterator end() const
  {
    return typename Traits::template Codim<codim>::template Partition<pitype>::Iterator(
        *globalGridView_, indexContainer_, true);
  }

  int boundaryId(const Intersection& /*intersection*/) const
  {
    DUNE_THROW(Dune::NotImplemented, "Call intersection.boundaryId() instead!");
    return -1;
  }

  int level() const
  {
    return globalGridView_->level();
  }

  template <class DataHandleImp, class DataType>
  void communicate(CommDataHandleIF<DataHandleImp, DataType>& /*data*/,
                   InterfaceType /*iftype*/,
                   CommunicationDirection /*dir*/) const
  {
    DUNE_THROW(Dune::NotImplemented,
               "As long as I am not sure what this does or is used for I will not implement this!");
    //    globalGridView_->communicate(data,iftype,dir);
  }

  const CollectiveCommunication& comm() const
  {
    return grid().comm();
  }

protected:
  const std::shared_ptr<const GlobalGridViewType> globalGridView_;
  const std::shared_ptr<const IndexContainerType> indexContainer_;
  const std::shared_ptr<const BoundaryInfoContainerType> boundaryInfoContainer_;
  const std::shared_ptr<const IndexSet> indexSet_;
};
} // namespace internal


template <class GlobalGridViewImp>
class SubdomainGridView : public internal::SubdomainGridViewCommon<internal::SubdomainGridViewTraits<GlobalGridViewImp>>
{
  using BaseType = internal::SubdomainGridViewCommon<internal::SubdomainGridViewTraits<GlobalGridViewImp>>;
  using ThisType = SubdomainGridView<GlobalGridViewImp>;

public:
  typedef internal::SubdomainGridViewTraits<GlobalGridViewImp> Traits;
  typedef typename Traits::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::GlobalGridViewType GlobalGridViewType;
  typedef typename BaseType::IndexType IndexType;
  typedef typename BaseType::IndexContainerType IndexContainerType;
  typedef typename BaseType::BoundaryInfoContainerType BoundaryInfoContainerType;

  SubdomainGridView(const std::shared_ptr<const GlobalGridViewType> globalGrdPrt,
                    const std::shared_ptr<const IndexContainerType> indexContainer,
                    const std::shared_ptr<const BoundaryInfoContainerType> boundaryInfoContainer)
    : BaseType(globalGrdPrt, indexContainer, boundaryInfoContainer)
  {
  }

  SubdomainGridView(const ThisType& other)
    : BaseType(other)
  {
  }

  SubdomainGridView(ThisType&& source)
    : BaseType(source)
  {
  }

  IntersectionIterator ibegin(const EntityType& ent) const
  {
    const IndexType& globalIndex = BaseType::globalGridView_->indexSet().index(ent);
    const typename BoundaryInfoContainerType::const_iterator result =
        BaseType::boundaryInfoContainer_->find(globalIndex);
    // if this is an entity at the boundary
    if (result != BaseType::boundaryInfoContainer_->end()) {
      // get the information for this entity
      const std::map<int, size_t>& info = result->second;
      // return wrapped iterator
      return IntersectionIterator(*BaseType::globalGridView_, ent, info);
    } else {
      // return iterator which just passes everything thrugh
      return IntersectionIterator(*BaseType::globalGridView_, ent);
    } // if this is an entity at the boundary
  } // IntersectionIteratorType ibegin(const EntityType& entity) const

  IntersectionIterator iend(const EntityType& ent) const
  {
    const IndexType& globalIndex = BaseType::globalGridView_->indexSet().index(ent);
    const typename BoundaryInfoContainerType::const_iterator result =
        BaseType::boundaryInfoContainer_->find(globalIndex);
    // if this is an entity at the boundary
    if (result != BaseType::boundaryInfoContainer_->end()) {
      // get the information for this entity
      const std::map<int, size_t>& info = result->second;
      // return wrapped iterator
      return IntersectionIterator(*BaseType::globalGridView_, ent, info, true);
    } else {
      // return iterator which just passes everything thrugh
      return IntersectionIterator(*BaseType::globalGridView_, ent, true);
    } // if this is an entity at the boundary
  }
}; // class SubdomainGridView


template <class GlobalGridViewImp>
class SubdomainCouplingGridView
    : public internal::SubdomainGridViewCommon<internal::SubdomainCouplingGridViewTraits<GlobalGridViewImp>>
{
  using BaseType = internal::SubdomainGridViewCommon<internal::SubdomainCouplingGridViewTraits<GlobalGridViewImp>>;

public:
  typedef SubdomainCouplingGridView<GlobalGridViewImp> ThisType;
  typedef internal::SubdomainCouplingGridViewTraits<GlobalGridViewImp> Traits;
  typedef typename Traits::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::GlobalGridViewType GlobalGridViewType;
  typedef typename BaseType::IndexType IndexType;
  typedef typename BaseType::IndexContainerType IndexContainerType;
  typedef typename BaseType::BoundaryInfoContainerType BoundaryInfoContainerType;
  typedef SubdomainGridView<GlobalGridViewType> InsideType;
  typedef SubdomainGridView<GlobalGridViewType> OutsideType;
  //! container type for the intersection information
  typedef std::map<IndexType, std::vector<int>> IntersectionInfoContainerType;

  SubdomainCouplingGridView(const std::shared_ptr<const GlobalGridViewType> globalGrdPart,
                            const std::shared_ptr<const IndexContainerType> indexContainer,
                            const std::shared_ptr<const IntersectionInfoContainerType> intersectionContainer,
                            const std::shared_ptr<const InsideType> insd,
                            const std::shared_ptr<const OutsideType> outsd)
    : BaseType(globalGrdPart,
               indexContainer,
               std::shared_ptr<const BoundaryInfoContainerType>(new BoundaryInfoContainerType()))
    , intersectionContainer_(intersectionContainer)
    , inside_(insd)
    , outside_(outsd)
  {
    assert(intersectionContainer_);
    assert(inside_);
    assert(outside_);
  }

  SubdomainCouplingGridView(const ThisType& other)
    : BaseType(other)
    , intersectionContainer_(other.intersectionContainer_)
    , inside_(other.inside_)
    , outside_(other.outside_)
  {
    assert(intersectionContainer_);
    assert(inside_);
    assert(outside_);
  }

  SubdomainCouplingGridView(ThisType&& source)
    : BaseType(source)
    , intersectionContainer_(std::move(source.intersectionContainer_))
    , inside_(std::move(source.inside_))
    , outside_(std::move(source.outside_))
  {
  }

  IntersectionIterator ibegin(const EntityType& ent) const
  {
    assert(intersectionContainer_);
    const IndexType& globalIndex = BaseType::globalGridView().indexSet().index(ent);
    const typename IntersectionInfoContainerType::const_iterator result = intersectionContainer_->find(globalIndex);
    assert(result != intersectionContainer_->end());
    // get the information for this entity
    const auto& info = result->second;
    // return localized iterator
    return IntersectionIterator(BaseType::globalGridView(), ent, info);
  } // IntersectionIteratorType ibegin(const EntityType& entity) const

  IntersectionIterator iend(const EntityType& ent) const
  {
    const IndexType& globalIndex = BaseType::globalGridView().indexSet().index(ent);
    const typename IntersectionInfoContainerType::const_iterator result = intersectionContainer_->find(globalIndex);
    assert(result != intersectionContainer_->end());
    // get the information for this entity
    const auto& info = result->second;
    // return localized iterator
    return IntersectionIterator(BaseType::globalGridView(), ent, info, true);
  } // IntersectionIteratorType iend(const EntityType& entity) const

  std::shared_ptr<const InsideType> inside() const
  {
    return inside_;
  }

  std::shared_ptr<const InsideType> outside() const
  {
    return outside_;
  }

private:
  const std::shared_ptr<const IntersectionInfoContainerType> intersectionContainer_;
  const std::shared_ptr<const InsideType> inside_;
  const std::shared_ptr<const OutsideType> outside_;
}; // class SubdomainCouplingGridView


template <class GlobalGridViewImp>
class SubdomainBoundaryGridView
    : public internal::SubdomainGridViewCommon<internal::SubdomainBoundaryGridViewTraits<GlobalGridViewImp>>
{
  using BaseType = internal::SubdomainGridViewCommon<internal::SubdomainBoundaryGridViewTraits<GlobalGridViewImp>>;

public:
  typedef SubdomainBoundaryGridView<GlobalGridViewImp> ThisType;
  typedef internal::SubdomainBoundaryGridViewTraits<GlobalGridViewImp> Traits;
  typedef typename Traits::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::GlobalGridViewType GlobalGridViewType;
  typedef typename BaseType::IndexType IndexType;
  typedef typename BaseType::IndexContainerType IndexContainerType;
  typedef typename BaseType::BoundaryInfoContainerType BoundaryInfoContainerType;
  typedef SubdomainGridView<GlobalGridViewType> InsideType;
  typedef SubdomainGridView<GlobalGridViewType> OutsideType;
  //! container type for the intersection information
  typedef std::map<IndexType, std::vector<int>> IntersectionInfoContainerType;

  SubdomainBoundaryGridView(const std::shared_ptr<const GlobalGridViewType> globalGrdPart,
                            const std::shared_ptr<const IndexContainerType> indexContainer,
                            const std::shared_ptr<const IntersectionInfoContainerType> intersectionContainer,
                            const std::shared_ptr<const InsideType> insd)
    : BaseType(globalGrdPart,
               indexContainer,
               std::shared_ptr<const BoundaryInfoContainerType>(new BoundaryInfoContainerType()))
    , intersectionContainer_(intersectionContainer)
    , inside_(insd)
  {
  }

  SubdomainBoundaryGridView(const ThisType& other)
    : BaseType(other)
    , intersectionContainer_(other.intersectionContainer_)
    , inside_(other.inside_)
  {
  }

  SubdomainBoundaryGridView(ThisType&& source)
    : BaseType(source)
    , intersectionContainer_(std::move(source.intersectionContainer_))
    , inside_(std::move(source.inside_))
  {
  }

  IntersectionIterator ibegin(const EntityType& ent) const
  {
    const IndexType& globalIndex = BaseType::globalGridView().indexSet().index(ent);
    const typename IntersectionInfoContainerType::const_iterator result = intersectionContainer_->find(globalIndex);
    assert(result != intersectionContainer_->end());
    // get the information for this entity
    const auto& info = result->second;
    // return localized iterator
    return IntersectionIterator(BaseType::globalGridView(), ent, info);
  } // IntersectionIteratorType ibegin(const EntityType& entity) const

  IntersectionIterator iend(const EntityType& ent) const
  {
    const IndexType& globalIndex = BaseType::globalGridView().indexSet().index(ent);
    const typename IntersectionInfoContainerType::const_iterator result = intersectionContainer_->find(globalIndex);
    assert(result != intersectionContainer_->end());
    // get the information for this entity
    const auto& info = result->second;
    // return localized iterator
    return IntersectionIterator(BaseType::globalGridView(), ent, info, true);
  } // IntersectionIteratorType iend(const EntityType& entity) const

  std::shared_ptr<const InsideType> inside() const
  {
    return inside_;
  }

private:
  const std::shared_ptr<const IntersectionInfoContainerType> intersectionContainer_;
  const std::shared_ptr<const InsideType> inside_;
}; // class SubdomainBoundaryGridView

} // namespace Grid
} // namespace XT

//! these are only necessary because iterators get fubared in copying into the default iterator range,
template <class GlobalGridViewImp, class Entity>
inline auto intersections(const GridView<XT::Grid::internal::SubdomainGridViewTraits<GlobalGridViewImp>>& gv,
                          const Entity& e)
    -> XT::Grid::SubdomainGridViewIntersectionRange<XT::Grid::SubdomainGridView<GlobalGridViewImp>>
{
  return XT::Grid::SubdomainGridViewIntersectionRange<XT::Grid::SubdomainGridView<GlobalGridViewImp>>(gv.impl(), e);
}

template <class GlobalGridViewImp, class Entity>
inline auto intersections(const GridView<XT::Grid::internal::SubdomainCouplingGridViewTraits<GlobalGridViewImp>>& gv,
                          const Entity& e)
    -> XT::Grid::SubdomainGridViewIntersectionRange<XT::Grid::SubdomainCouplingGridView<GlobalGridViewImp>>
{
  return XT::Grid::SubdomainGridViewIntersectionRange<XT::Grid::SubdomainCouplingGridView<GlobalGridViewImp>>(gv.impl(),
                                                                                                              e);
}

template <class GlobalGridViewImp, class Entity>
inline auto intersections(const GridView<XT::Grid::internal::SubdomainBoundaryGridViewTraits<GlobalGridViewImp>>& gv,
                          const Entity& e)
    -> XT::Grid::SubdomainGridViewIntersectionRange<XT::Grid::SubdomainBoundaryGridView<GlobalGridViewImp>>
{
  return XT::Grid::SubdomainGridViewIntersectionRange<XT::Grid::SubdomainBoundaryGridView<GlobalGridViewImp>>(gv.impl(),
                                                                                                              e);
}

} // namespace Dune

#endif // DUNE_XT_GRID_VIEW_SUBDOMAIN_PART_HH
