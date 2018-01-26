// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_GRID_VIEW_SUBDOMAIN_PART_HH
#define DUNE_XT_GRID_VIEW_SUBDOMAIN_PART_HH

#include <map>
#include <set>
#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/capabilities.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/common/gridpart.hh>
#endif

#include "entity-iterator.hh"
#include "indexset.hh"
#include "intersection-iterator.hh"
#include "intersection-wrapper.hh"

namespace Dune {
namespace XT {
namespace Grid {


// forwards, needed for the traits below
template <class GlobalGridPartImp>
class SubdomainGridPart;

template <class GlobalGridPartImp>
class SubdomainCouplingGridPart;

template <class GlobalGridPartImp>
class SubdomainBoundaryGridPart;


namespace internal {


template <class GlobalGridPartImp>
class SubdomainGridPartTraits
{
public:
  typedef GlobalGridPartImp GlobalGridPartType;

  typedef SubdomainGridPart<GlobalGridPartImp> GridPartType;
  typedef typename GlobalGridPartType::GridType GridType;
  typedef typename internal::IndexBasedIndexSet<GlobalGridPartType> IndexSetType;
  typedef typename GlobalGridPartType::CollectiveCommunicationType CollectiveCommunicationType;
  typedef typename GlobalGridPartType::TwistUtilityType TwistUtilityType;

  static const PartitionIteratorType indexSetPartitionType = GlobalGridPartType::indexSetPartitionType;
  static const InterfaceType indexSetInterfaceType = GlobalGridPartType::indexSetInterfaceType;

  typedef internal::FakeDomainBoundaryIntersectionIterator<GlobalGridPartType> IntersectionIteratorType;

  template <int codim>
  struct Codim : public GlobalGridPartType::template Codim<codim>
  {
    template <PartitionIteratorType pitype>
    struct Partition
    {
      typedef typename internal::IndexBasedEntityIterator<GlobalGridPartType, codim, pitype> IteratorType;
    };
  };

  static const bool conforming = GlobalGridPartType::Traits::conforming;
}; // class SubdomainGridPartTraits


template <class GlobalGridPartImp>
struct SubdomainCouplingGridPartTraits : public SubdomainGridPartTraits<GlobalGridPartImp>
{
  typedef GlobalGridPartImp GlobalGridPartType;
  typedef SubdomainCouplingGridPart<GlobalGridPartImp> GridPartType;
  typedef internal::LocalIntersectionIterator<GlobalGridPartType> IntersectionIteratorType;
};


template <class GlobalGridPartImp>
struct SubdomainBoundaryGridPartTraits : public SubdomainGridPartTraits<GlobalGridPartImp>
{
  typedef GlobalGridPartImp GlobalGridPartType;
  typedef SubdomainBoundaryGridPart<GlobalGridPartImp> GridPartType;
  typedef internal::LocalIntersectionIterator<GlobalGridPartType> IntersectionIteratorType;
}; // class SubdomainBoundaryGridPartTraits


template <class SubdomainGridPartImp>
class SubdomainGridPartIntersectionRange
{
  typedef typename SubdomainGridPartImp::EntityType EntityType;
  typedef typename SubdomainGridPartImp::IntersectionIteratorType IntersectionIteratorType;

public:
  SubdomainGridPartIntersectionRange(const SubdomainGridPartImp& layer, const EntityType& entity)
    : layer_(layer)
    , entity_(entity)
  {
  }

  IntersectionIteratorType begin() const
  {
    return layer_.ibegin(entity_);
  }

  IntersectionIteratorType end() const
  {
    return layer_.iend(entity_);
  }

  const SubdomainGridPartImp& layer_;
  const EntityType& entity_;
}; // class SubdomainGridPartIntersectionRange

template <class Traits>
class SubdomainGridPartCommon
#if HAVE_DUNE_FEM
    : public Fem::GridPartInterface<Traits>
#endif
{
  using ThisType = SubdomainGridPartCommon<Traits>;

protected:
#if HAVE_DUNE_FEM
  typedef Fem::GridPartInterface<Traits> BaseType;
  typedef BaseType BaseTraits;
#else
  typedef Traits BaseTraits;
  typedef ThisType BaseType;
#endif
public:
  typedef typename Traits::GridType GridType;
  typedef typename Traits::CollectiveCommunicationType CollectiveCommunicationType;
  typedef typename Traits::GlobalGridPartType GlobalGridPartType;
  typedef typename Traits::IndexSetType IndexSetType;
  typedef typename Traits::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename GridType::template Codim<0>::Entity EntityType;

  typedef typename IndexSetType::IndexType IndexType;
  typedef std::map<IndexType, IndexType> IndexMapType;
  typedef Dune::GeometryType GeometryType;
  //! container type for the indices
  typedef std::map<GeometryType, IndexMapType> IndexContainerType;
  //! container type for the boundary information
  typedef std::map<IndexType, std::map<int, size_t>> BoundaryInfoContainerType;

public:
  SubdomainGridPartCommon(const std::shared_ptr<const GlobalGridPartType> globalGrdPrt,
                          const std::shared_ptr<const IndexContainerType> indexContainer,
                          const std::shared_ptr<const BoundaryInfoContainerType> boundaryInfoContainer)
    : globalGridPart_(globalGrdPrt)
    , indexContainer_(indexContainer)
    , boundaryInfoContainer_(boundaryInfoContainer)
    , indexSet_(*globalGridPart_, indexContainer_)
  {
  }

  SubdomainGridPartCommon(const ThisType& other) = default;
  SubdomainGridPartCommon(ThisType&& source) = default;

  const IndexSetType& indexSet() const
  {
    return indexSet_;
  }

  const GridType& grid() const
  {
    return globalGridPart_->grid();
  }

  const GlobalGridPartType& globalGridPart() const
  {
    return *globalGridPart_;
  }

  template <int codim>
  typename BaseTraits::template Codim<codim>::IteratorType begin() const
  {
    return typename BaseTraits::template Codim<codim>::IteratorType(*globalGridPart_, indexContainer_);
  }

  template <int codim, PartitionIteratorType pitype>
  typename BaseTraits::template Codim<codim>::template Partition<pitype>::IteratorType begin() const
  {
    return typename BaseTraits::template Codim<codim>::template Partition<pitype>::IteratorType(*globalGridPart_,
                                                                                                indexContainer_);
  }

  template <int codim>
  typename BaseTraits::template Codim<codim>::IteratorType end() const
  {
    return typename BaseTraits::template Codim<codim>::IteratorType(*globalGridPart_, indexContainer_, true);
  }

  template <int codim, PartitionIteratorType pitype>
  typename BaseTraits::template Codim<codim>::template Partition<pitype>::IteratorType end() const
  {
    return typename BaseTraits::template Codim<codim>::template Partition<pitype>::IteratorType(
        *globalGridPart_, indexContainer_, true);
  }

  int boundaryId(const IntersectionType& intersection) const
  {
    DUNE_THROW(Dune::NotImplemented, "Call intersection.boundaryId() instead!");
    return -1;
  }

  int level() const
  {
    return globalGridPart_->level();
  }

  template <class DataHandleImp, class DataType>
  void communicate(CommDataHandleIF<DataHandleImp, DataType>& /*data*/,
                   InterfaceType /*iftype*/,
                   CommunicationDirection /*dir*/) const
  {
    DUNE_THROW(Dune::NotImplemented,
               "As long as I am not sure what this does or is used for I will not implement this!");
    //    globalGridPart_->communicate(data,iftype,dir);
  }

  const CollectiveCommunicationType& comm() const
  {
    return grid().comm();
  }

protected:
  const std::shared_ptr<const GlobalGridPartType> globalGridPart_;
  const std::shared_ptr<const IndexContainerType> indexContainer_;
  const std::shared_ptr<const BoundaryInfoContainerType> boundaryInfoContainer_;
  const IndexSetType indexSet_;
};
} // namespace internal


template <class GlobalGridPartImp>
class SubdomainGridPart : public internal::SubdomainGridPartCommon<internal::SubdomainGridPartTraits<GlobalGridPartImp>>
{
  using BaseType = internal::SubdomainGridPartCommon<internal::SubdomainGridPartTraits<GlobalGridPartImp>>;
  using ThisType = SubdomainGridPart<GlobalGridPartImp>;

public:
  typedef internal::SubdomainGridPartTraits<GlobalGridPartImp> Traits;
  typedef typename Traits::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::GlobalGridPartType GlobalGridPartType;
  typedef typename BaseType::IndexType IndexType;
  typedef typename BaseType::IndexContainerType IndexContainerType;
  typedef typename BaseType::BoundaryInfoContainerType BoundaryInfoContainerType;

  SubdomainGridPart(const std::shared_ptr<const GlobalGridPartType> globalGrdPrt,
                    const std::shared_ptr<const IndexContainerType> indexContainer,
                    const std::shared_ptr<const BoundaryInfoContainerType> boundaryInfoContainer)
    : BaseType(globalGrdPrt, indexContainer, boundaryInfoContainer)
  {
  }

  SubdomainGridPart(const ThisType& other) = default;

  SubdomainGridPart(ThisType&& source) = default;

  IntersectionIteratorType ibegin(const EntityType& ent) const
  {
    const IndexType& globalIndex = BaseType::globalGridPart_->indexSet().index(ent);
    const typename BoundaryInfoContainerType::const_iterator result =
        BaseType::boundaryInfoContainer_->find(globalIndex);
    // if this is an entity at the boundary
    if (result != BaseType::boundaryInfoContainer_->end()) {
      // get the information for this entity
      const std::map<int, size_t>& info = result->second;
      // return wrapped iterator
      return IntersectionIteratorType(*BaseType::globalGridPart_, ent, info);
    } else {
      // return iterator which just passes everything thrugh
      return IntersectionIteratorType(*BaseType::globalGridPart_, ent);
    } // if this is an entity at the boundary
  } // IntersectionIteratorType ibegin(const EntityType& entity) const

  IntersectionIteratorType iend(const EntityType& ent) const
  {
    const IndexType& globalIndex = BaseType::globalGridPart_->indexSet().index(ent);
    const typename BoundaryInfoContainerType::const_iterator result =
        BaseType::boundaryInfoContainer_->find(globalIndex);
    // if this is an entity at the boundary
    if (result != BaseType::boundaryInfoContainer_->end()) {
      // get the information for this entity
      const std::map<int, size_t>& info = result->second;
      // return wrapped iterator
      return IntersectionIteratorType(*BaseType::globalGridPart_, ent, info, true);
    } else {
      // return iterator which just passes everything thrugh
      return IntersectionIteratorType(*BaseType::globalGridPart_, ent, true);
    } // if this is an entity at the boundary
  }
}; // class SubdomainGridPart


template <class GlobalGridPartImp>
class SubdomainCouplingGridPart
    : public internal::SubdomainGridPartCommon<internal::SubdomainCouplingGridPartTraits<GlobalGridPartImp>>
{
  using BaseType = internal::SubdomainGridPartCommon<internal::SubdomainCouplingGridPartTraits<GlobalGridPartImp>>;

public:
  typedef SubdomainCouplingGridPart<GlobalGridPartImp> ThisType;
  typedef internal::SubdomainCouplingGridPartTraits<GlobalGridPartImp> Traits;
  typedef typename Traits::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::GlobalGridPartType GlobalGridPartType;
  typedef typename BaseType::IndexType IndexType;
  typedef typename BaseType::IndexContainerType IndexContainerType;
  typedef typename BaseType::BoundaryInfoContainerType BoundaryInfoContainerType;
  typedef SubdomainGridPart<GlobalGridPartType> InsideType;
  typedef SubdomainGridPart<GlobalGridPartType> OutsideType;
  //! container type for the intersection information
  typedef std::map<IndexType, std::vector<int>> IntersectionInfoContainerType;

  SubdomainCouplingGridPart(const std::shared_ptr<const GlobalGridPartType> globalGrdPart,
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
  }

  SubdomainCouplingGridPart(const ThisType& other) = default;

  SubdomainCouplingGridPart(ThisType&& source) = default;

  IntersectionIteratorType ibegin(const EntityType& ent) const
  {
    const IndexType& globalIndex = BaseType::globalGridPart().indexSet().index(ent);
    const typename IntersectionInfoContainerType::const_iterator result = intersectionContainer_->find(globalIndex);
    assert(result != intersectionContainer_->end());
    // get the information for this entity
    const auto& info = result->second;
    // return localized iterator
    return IntersectionIteratorType(BaseType::globalGridPart(), ent, info);
  } // IntersectionIteratorType ibegin(const EntityType& entity) const

  IntersectionIteratorType iend(const EntityType& ent) const
  {
    const IndexType& globalIndex = BaseType::globalGridPart().indexSet().index(ent);
    const typename IntersectionInfoContainerType::const_iterator result = intersectionContainer_->find(globalIndex);
    assert(result != intersectionContainer_->end());
    // get the information for this entity
    const auto& info = result->second;
    // return localized iterator
    return IntersectionIteratorType(BaseType::globalGridPart(), ent, info, true);
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
}; // class SubdomainCouplingGridPart


template <class GlobalGridPartImp>
class SubdomainBoundaryGridPart
    : public internal::SubdomainGridPartCommon<internal::SubdomainBoundaryGridPartTraits<GlobalGridPartImp>>
{
  using BaseType = internal::SubdomainGridPartCommon<internal::SubdomainBoundaryGridPartTraits<GlobalGridPartImp>>;

public:
  typedef SubdomainBoundaryGridPart<GlobalGridPartImp> ThisType;
  typedef internal::SubdomainBoundaryGridPartTraits<GlobalGridPartImp> Traits;
  typedef typename Traits::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::GlobalGridPartType GlobalGridPartType;
  typedef typename BaseType::IndexType IndexType;
  typedef typename BaseType::IndexContainerType IndexContainerType;
  typedef typename BaseType::BoundaryInfoContainerType BoundaryInfoContainerType;
  typedef SubdomainGridPart<GlobalGridPartType> InsideType;
  typedef SubdomainGridPart<GlobalGridPartType> OutsideType;
  //! container type for the intersection information
  typedef std::map<IndexType, std::vector<int>> IntersectionInfoContainerType;

  SubdomainBoundaryGridPart(const std::shared_ptr<const GlobalGridPartType> globalGrdPart,
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

  SubdomainBoundaryGridPart(const ThisType& other) = default;

  SubdomainBoundaryGridPart(ThisType&& source) = default;

  IntersectionIteratorType ibegin(const EntityType& ent) const
  {
    const IndexType& globalIndex = BaseType::globalGridPart().indexSet().index(ent);
    const typename IntersectionInfoContainerType::const_iterator result = intersectionContainer_->find(globalIndex);
    assert(result != intersectionContainer_->end());
    // get the information for this entity
    const auto& info = result->second;
    // return localized iterator
    return IntersectionIteratorType(BaseType::globalGridPart(), ent, info);
  } // IntersectionIteratorType ibegin(const EntityType& entity) const

  IntersectionIteratorType iend(const EntityType& ent) const
  {
    const IndexType& globalIndex = BaseType::globalGridPart().indexSet().index(ent);
    const typename IntersectionInfoContainerType::const_iterator result = intersectionContainer_->find(globalIndex);
    assert(result != intersectionContainer_->end());
    // get the information for this entity
    const auto& info = result->second;
    // return localized iterator
    return IntersectionIteratorType(BaseType::globalGridPart(), ent, info, true);
  } // IntersectionIteratorType iend(const EntityType& entity) const

  std::shared_ptr<const InsideType> inside() const
  {
    return inside_;
  }

private:
  const std::shared_ptr<const IntersectionInfoContainerType> intersectionContainer_;
  const std::shared_ptr<const InsideType> inside_;
}; // class SubdomainBoundaryGridPart

} // namespace Grid
} // namespace XT

#if HAVE_DUNE_FEM

namespace Fem {
namespace GridPartCapabilities {


template <class GridPartType>
struct hasGrid<XT::Grid::SubdomainGridPart<GridPartType>>
{
  static const bool v = hasGrid<GridPartType>::v;
};

template <class GridPartType>
struct hasSingleGeometryType<XT::Grid::SubdomainGridPart<GridPartType>>
{
  static const bool v = hasSingleGeometryType<GridPartType>::v;
  static const unsigned int topologyId = hasSingleGeometryType<GridPartType>::topologyId;
};

template <class GridPartType>
struct isCartesian<XT::Grid::SubdomainGridPart<GridPartType>>
{
  static const bool v = isCartesian<GridPartType>::v;
};

template <class GridPartType, int codim>
struct hasEntity<XT::Grid::SubdomainGridPart<GridPartType>, codim>
{
  static const bool v = hasEntity<GridPartType, codim>::v;
};

template <class GridPartType, int codim>
struct canCommunicate<XT::Grid::SubdomainGridPart<GridPartType>, codim>
{
  static const bool v = false;
};

template <class GridPartType>
struct isConforming<XT::Grid::SubdomainGridPart<GridPartType>>
{
  static const bool v = false;
};


template <class GridPartType>
struct hasGrid<XT::Grid::SubdomainCouplingGridPart<GridPartType>>
{
  static const bool v = hasGrid<GridPartType>::v;
};

template <class GridPartType>
struct hasSingleGeometryType<XT::Grid::SubdomainCouplingGridPart<GridPartType>>
{
  static const bool v = hasSingleGeometryType<GridPartType>::v;
  static const unsigned int topologyId = hasSingleGeometryType<GridPartType>::topologyId;
};

template <class GridPartType>
struct isCartesian<XT::Grid::SubdomainCouplingGridPart<GridPartType>>
{
  static const bool v = isCartesian<GridPartType>::v;
};

template <class GridPartType, int codim>
struct hasEntity<XT::Grid::SubdomainCouplingGridPart<GridPartType>, codim>
{
  static const bool v = hasEntity<GridPartType, codim>::v;
};

template <class GridPartType, int codim>
struct canCommunicate<XT::Grid::SubdomainCouplingGridPart<GridPartType>, codim>
{
  static const bool v = false;
};

template <class GridPartType>
struct isConforming<XT::Grid::SubdomainCouplingGridPart<GridPartType>>
{
  static const bool v = false;
};


template <class GridPartType>
struct hasGrid<XT::Grid::SubdomainBoundaryGridPart<GridPartType>>
{
  static const bool v = hasGrid<GridPartType>::v;
};

template <class GridPartType>
struct hasSingleGeometryType<XT::Grid::SubdomainBoundaryGridPart<GridPartType>>
{
  static const bool v = hasSingleGeometryType<GridPartType>::v;
  static const unsigned int topologyId = hasSingleGeometryType<GridPartType>::topologyId;
};

template <class GridPartType>
struct isCartesian<XT::Grid::SubdomainBoundaryGridPart<GridPartType>>
{
  static const bool v = isCartesian<GridPartType>::v;
};

template <class GridPartType, int codim>
struct hasEntity<XT::Grid::SubdomainBoundaryGridPart<GridPartType>, codim>
{
  static const bool v = hasEntity<GridPartType, codim>::v;
};

template <class GridPartType, int codim>
struct canCommunicate<XT::Grid::SubdomainBoundaryGridPart<GridPartType>, codim>
{
  static const bool v = false;
};

template <class GridPartType>
struct isConforming<XT::Grid::SubdomainBoundaryGridPart<GridPartType>>
{
  static const bool v = false;
};


} // namespace GridPartCapabilities
} // namespace Fem

#endif // HAVE_DUNE_FEM


template <typename GP, class Entity>
inline auto intersections(const XT::Grid::SubdomainGridPart<GP>& gv, const Entity& e)
    -> XT::Grid::internal::SubdomainGridPartIntersectionRange<XT::Grid::SubdomainGridPart<GP>>
{
  return XT::Grid::internal::SubdomainGridPartIntersectionRange<XT::Grid::SubdomainGridPart<GP>>(gv, e);
}


template <typename GP, class Entity>
inline auto intersections(const XT::Grid::SubdomainBoundaryGridPart<GP>& gv, const Entity& e)
    -> XT::Grid::internal::SubdomainGridPartIntersectionRange<XT::Grid::SubdomainBoundaryGridPart<GP>>
{
  return XT::Grid::internal::SubdomainGridPartIntersectionRange<XT::Grid::SubdomainBoundaryGridPart<GP>>(gv, e);
}


template <typename GP, class Entity>
inline auto intersections(const XT::Grid::SubdomainCouplingGridPart<GP>& gv, const Entity& e)
    -> XT::Grid::internal::SubdomainGridPartIntersectionRange<XT::Grid::SubdomainCouplingGridPart<GP>>
{
  return XT::Grid::internal::SubdomainGridPartIntersectionRange<XT::Grid::SubdomainCouplingGridPart<GP>>(gv, e);
}


} // namespace Dune

#endif // DUNE_XT_GRID_VIEW_SUBDOMAIN_PART_HH
