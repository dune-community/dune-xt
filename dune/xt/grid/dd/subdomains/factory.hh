// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2017 - 2018)

#ifndef DUNE_XT_GRID_DD_SUBDOMAINS_FACTORY_HH
#define DUNE_XT_GRID_DD_SUBDOMAINS_FACTORY_HH

#include <memory>
#include <vector>
#include <map>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>

#include <dune/xt/common/logging.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/grids.hh>

#include <dune/xt/grid/view/subdomain/view.hh>

#include "grid.hh"

namespace Dune {
namespace XT {
namespace Grid {
namespace DD {
namespace internal {


/// \todo: collect all specializations below into this implementation, differentiate at runtime using std::is_same,
/// allow user to override via template specialization
template <class GridImp>
class NeighborRecursionLevel
{
  static_assert(AlwaysFalse<GridImp>::value, "Please add an appropriate specialization for this GridImp!");

public:
  static size_t compute() = delete;
};

// OneDGrid
template <>
class NeighborRecursionLevel<OneDGrid>
{
public:
  static size_t compute()
  {
    return 1;
  }
};

// YaspGrid
template <class Coordinates>
class NeighborRecursionLevel<YaspGrid<1, Coordinates>>
{
public:
  static size_t compute()
  {
    return 1;
  }
};

template <class Coordinates>
class NeighborRecursionLevel<YaspGrid<2, Coordinates>>
{
public:
  static size_t compute()
  {
    return 2;
  }
};

template <class Coordinates>
class NeighborRecursionLevel<YaspGrid<3, Coordinates>>
{
public:
  static size_t compute()
  {
    return 3;
  }
};

#if HAVE_DUNE_SPGRID
template <class ct, int dim, template <int> class Ref, class Comm>
class NeighborRecursionLevel<SPGrid<ct, dim, Ref, Comm>>
{
public:
  static size_t compute()
  {
    return NeighborRecursionLevel<YaspGrid<dim>>::compute();
  }
};
#endif // HAVE_DUNE_SPGRID

// ALUGrid
#if HAVE_DUNE_ALUGRID

template <Dune::ALUGridRefinementType ref, class Comm>
class NeighborRecursionLevel<Dune::ALUGrid<2, 2, Dune::cube, ref, Comm>>
{
public:
  static size_t compute()
  {
    return 2;
  }
};

template <Dune::ALUGridRefinementType ref, class Comm>
class NeighborRecursionLevel<Dune::ALUGrid<2, 2, Dune::simplex, ref, Comm>>
{
public:
  static size_t compute()
  {
    return 3;
  }
};

template <class Comm, ALUGridRefinementType tp>
class NeighborRecursionLevel<Dune::ALUGrid<3, 3, Dune::simplex, tp, Comm>>
{
public:
  static size_t compute()
  {
    return 9;
  } // just a guess!
};

template <class Comm, ALUGridRefinementType tp>
class NeighborRecursionLevel<Dune::ALUGrid<3, 3, Dune::cube, tp, Comm>>
{
public:
  static size_t compute()
  {
    return 3;
  } // just a guess!
};

#endif // HAVE_DUNE_ALUGRID
#if HAVE_DUNE_UGGRID || HAVE_UG

template <>
class NeighborRecursionLevel<UGGrid<2>>
{
public:
  static size_t compute()
  {
    return 3; // just a guess!
  }
};

template <>
class NeighborRecursionLevel<UGGrid<3>>
{
public:
  static size_t compute()
  {
    return 9; // just a guess!
  }
};

#endif // HAVE_DUNE_UGGRID || HAVE_UG
#if HAVE_ALBERTA

template <>
class NeighborRecursionLevel<AlbertaGrid<1, 1>>
{
public:
  static size_t compute()
  {
    return 1;
  }
};

template <>
class NeighborRecursionLevel<AlbertaGrid<2, 2>>
{
public:
  static size_t compute()
  {
    return 3; // just a guess!
  }
};

template <>
class NeighborRecursionLevel<AlbertaGrid<3, 3>>
{
public:
  static size_t compute()
  {
    return 9; // just a guess!
  }
};

#endif // HAVE_ALBERTA

} // namespace internal


template <class GridImp>
class SubdomainGridFactory
{
public:
  typedef GridImp GridType;
  typedef SubdomainGridFactory<GridType> ThisType;
  static const unsigned int dim = GridType::dimension;
  typedef SubdomainGrid<GridType> DdGridType;
  typedef typename GridType::template Codim<0>::Entity EntityType;

private:
  typedef typename DdGridType::GlobalGridViewType GlobalGridViewType;
  typedef typename DdGridType::LocalGridViewType LocalGridViewType;
  typedef typename DdGridType::BoundaryGridViewType BoundaryGridViewType;
  typedef typename DdGridType::CouplingGridViewType CouplingGridViewType;
  typedef typename GlobalGridViewType::IndexSet::IndexType IndexType;
  typedef Dune::GeometryType GeometryType;
  // i.e. maps a local to a globl index
  typedef std::map<IndexType, IndexType> IndexMapType;
  // i.e. maps a GeometryType to a map of local and global indices
  typedef std::map<GeometryType, IndexMapType> GeometryMapType;
  // i.e. contains a GeometryMap for each subdomain
  typedef std::map<size_t, std::shared_ptr<GeometryMapType>> SubdomainMapType;
  // map type which maps from an entity index (of the global grid parts index set) to a subdomain
  typedef std::map<IndexType, size_t> EntityToSubdomainMapType;
  typedef FieldVector<size_t, dim + 1> CodimSizesType;
  // for the neighbor information between the subdomains
  typedef std::set<size_t> NeighboringSubdomainsSetType;

  template <int c, int d>
  struct Add
  {
    static void subEntities(ThisType& factory,
                            const EntityType& entity,
                            GeometryMapType& geometryMap,
                            CodimSizesType& localCodimSizes)
    {
      // loop over all codim c subentities of the entity
      for (unsigned int i = 0; i < entity.subEntities(c); ++i) {
        const auto codimCentity = entity.template subEntity<c>(i);
        const GeometryType& geometryType = codimCentity.type();
        const IndexType globalIndex = factory.globalGridView_->indexSet().index(codimCentity);
        factory.addGeometryAndIndex(geometryMap, localCodimSizes, geometryType, globalIndex);
      }
      // add all codim c + 1 subentities
      Add<c + 1, d>::subEntities(factory, entity, geometryMap, localCodimSizes);
    } // ... subEntities(...)
  }; // struct Add

public:
  SubdomainGridFactory(const GridType& grid,
                       const size_t boundary_segment_index = std::numeric_limits<size_t>::max() - 42)
    : grid_(stackobject_to_shared_ptr(grid))
    , boundary_segment_index_(boundary_segment_index)
    , prepared_(false)
    , finalized_(false)
    , size_(0)
    , oversampled_(false)
  {
  }

  SubdomainGridFactory(const std::shared_ptr<const GridType> grid,
                       const size_t boundary_segment_index = std::numeric_limits<size_t>::max() - 42)
    : grid_(grid)
    , boundary_segment_index_(boundary_segment_index)
    , prepared_(false)
    , finalized_(false)
    , size_(0)
    , oversampled_(false)
  {
  }

  void prepare()
  {
    if (!prepared_) {
      globalGridView_ = std::make_shared<const GlobalGridViewType>(const_cast<GridType&>(*grid_));
      entityToSubdomainMap_ = std::shared_ptr<EntityToSubdomainMapType>(new EntityToSubdomainMapType());
      prepared_ = true;
    }
  } // ... prepare()

  const std::shared_ptr<const GlobalGridViewType> globalGridView() const
  {
    assert(prepared_ && "Please call prepare() before calling globalGridView()!");
    return globalGridView_;
  }

  void add(const EntityType& entity, const size_t subdomain)
  {
    // prepare
    assert(prepared_ && "Please call prepare() before calling add()!");
    assert(!finalized_ && "Do not call add() after calling finalized()!");
    const IndexType globalIndex = globalGridView_->indexSet().index(entity);
    // add subdomain to this entity index
    typename EntityToSubdomainMapType::iterator indexIt = entityToSubdomainMap_->find(globalIndex);
    if (indexIt == entityToSubdomainMap_->end()) {
      entityToSubdomainMap_->insert(std::pair<IndexType, size_t>(globalIndex, subdomain));
    } else {
      if (indexIt->second != subdomain)
        DUNE_THROW(InvalidStateException, "can not add entity to more than one subdomain!");
    }
    // create geometry map for this subdomain if needed (doing this explicitly (instead of just using insert()) only to
    // increment size)
    if (subdomainToEntityMap_.find(subdomain) == subdomainToEntityMap_.end()) {
      subdomainToEntityMap_.insert(std::pair<size_t, std::shared_ptr<GeometryMapType>>(
          subdomain, std::shared_ptr<GeometryMapType>(new GeometryMapType())));
      ++size_;
    }
    // create local codim sizes for this subdomain
    localCodimSizes_.insert(std::pair<size_t, CodimSizesType>(subdomain, CodimSizesType(0)));
    // add this entity and all subentities to the geometry map (geometry map has to exist, see above)
    assert(subdomainToEntityMap_.find(subdomain) != subdomainToEntityMap_.end() && "This should not happen!");
    GeometryMapType& geometryMap = *(subdomainToEntityMap_.find(subdomain)->second);
    // add geometry and global index of this codim 0 entity
    const GeometryType& geometryType = entity.type();
    CodimSizesType& localCodimSizes = localCodimSizes_.find(subdomain)->second;
    addGeometryAndIndex(geometryMap, localCodimSizes, geometryType, globalIndex);
    // add all remaining codims
    Add<1, dim>::subEntities(*this, entity, geometryMap, localCodimSizes);
  } // ... add(...)

  void finalize(const size_t oversamplingLayers = 0,
                const size_t neighbor_recursion_level = internal::NeighborRecursionLevel<GridType>::compute(),
                bool assert_connected = true)
  {
    assert(prepared_ && "Please call prepare() and add() before calling finalize()!");
    if (finalized_)
      return;

    // prepare
    // init data structures
    // for the subdomains inner boundaries
    //   * to map the local intersection index to the desired fake boundary id
    typedef std::map<int, size_t> IntersectionToBoundaryIdMapType;
    //   * to map the global entity index to one of those maps
    typedef std::map<IndexType, IntersectionToBoundaryIdMapType> EntityToIntersectionInfoMapType;
    //   * to hold one of those maps for each subdomain
    std::vector<std::shared_ptr<EntityToIntersectionInfoMapType>> subdomainInnerBoundaryInfos(size_);
    // for the coupling grid parts
    //   * vector to hold a (neighboring subdomain -> entity) map for each subdomain
    typename std::vector<SubdomainMapType> couplingMaps(size_, SubdomainMapType());
    //   * vector to hold a map of coupling sizes
    std::vector<std::map<size_t, CodimSizesType>> couplingCodimSizeMaps(size_, std::map<size_t, CodimSizesType>());
    //   * set of local coupling intersections
    typedef std::vector<int> IntersectionInfoSetType;
    //   * map to hold the above information for each coupling entity
    typedef std::map<IndexType, IntersectionInfoSetType> EntityToIntersectionSetMapType;
    //   * map to hold the above map for each neighboring subdomain
    typedef std::map<size_t, std::shared_ptr<EntityToIntersectionSetMapType>> CouplingIntersectionMapType;
    //   * vector to hold the above map for each subdomain
    std::vector<CouplingIntersectionMapType> couplingBoundaryInfos(size_, CouplingIntersectionMapType());
    // for the boundary grid parts
    //   * a map to hold the global entity ids for each subdomain
    typename std::map<size_t, std::shared_ptr<GeometryMapType>> boundaryGeometryMapMap;
    //   * a map to hold the codim sizes
    typename std::map<size_t, CodimSizesType> boundaryCodimSizesMap;
    //   * a map to hold the intersection informations
    std::map<size_t, std::shared_ptr<EntityToIntersectionSetMapType>> boundaryInfoMap;
    // for the neighboring information
    neighboringSubdomainSets_ = std::shared_ptr<std::vector<NeighboringSubdomainsSetType>>(
        new std::vector<NeighboringSubdomainsSetType>(size_, NeighboringSubdomainsSetType()));
    std::vector<NeighboringSubdomainsSetType>& neighboringSubdomainSets = *neighboringSubdomainSets_;
    // loop over all subdomains
    //   * to test for consecutive numbering
    //   * to compute the number of codim 0 entities per subdomain
    //   * to initialize data structures
    std::vector<size_t> subdomainSizes(size_, 0);
    for (size_t subdomain = 0; subdomain < size_; ++subdomain) {
      // test if this subdomain exists
      typename SubdomainMapType::iterator subdomainToEntityMapIt = subdomainToEntityMap_.find(subdomain);
      if (subdomainToEntityMapIt == subdomainToEntityMap_.end()) {
        DUNE_THROW(InvalidStateException, "numbering of subdomains has to be consecutive upon calling finalize()!");
      } else {
        // compute number of codim 0 entities
        const GeometryMapType& subdomainGeometryMap = *(subdomainToEntityMapIt->second);
        // loop over all geometry types of this subdomain
        for (typename GeometryMapType::const_iterator subdomainGeometryMapIt = subdomainGeometryMap.begin();
             subdomainGeometryMapIt != subdomainGeometryMap.end();
             ++subdomainGeometryMapIt) {
          // get geometry
          const GeometryType& geometryType = subdomainGeometryMapIt->first;
          // if this is a codim 0 geometry
          if (geometryType.dim() == dim) {
            // get the index map of this geometry type
            const IndexMapType& indexMap = subdomainGeometryMapIt->second;
            subdomainSizes[subdomain] += indexMap.size();
          } // if this is a codim 0 geometry
        } // loop over all geometry types of this subdomain
        // init data structures
        subdomainInnerBoundaryInfos[subdomain] =
            std::shared_ptr<EntityToIntersectionInfoMapType>(new EntityToIntersectionInfoMapType());
      } // test if this subdomain exists
    } // loop over all subdomains
    // walk the global grid part
    //   * to generate the information which sudomains neighbor each other
    for (auto entityIt = globalGridView_->template begin<0>(); entityIt != globalGridView_->template end<0>();
         ++entityIt) {
      // find the subdomains this entity lives in
      const EntityType& entity = *entityIt;
      const IndexType entityGlobalIndex = globalGridView_->indexSet().index(entity);
      const size_t entitySubdomain = getSubdomainOf(entityGlobalIndex);
      // get the set of this subdomains neighbors
      NeighboringSubdomainsSetType& neighborsOfSubdomain = neighboringSubdomainSets[entitySubdomain];
      // get the boundary info map for this subdomain
      EntityToIntersectionInfoMapType& subdomainInnerBoundaryInfo = *(subdomainInnerBoundaryInfos[entitySubdomain]);
      // walk the neighbors
      bool subdomainsEntitiesAreConnected = false;
      for (auto&& intersectionIt = globalGridView_->ibegin(entity); intersectionIt != globalGridView_->iend(entity);
           ++intersectionIt) {
        const auto& intersection = *intersectionIt;
        // check the type of this intersection
        if (intersection.boundary() && !intersection.neighbor()) {
          // get local index of the intersection
          const int intersectionLocalIndex = intersection.indexInInside();
          // for the boundary grid part
          //   * get the maps for this subdomain (and create them, if necessary)
          if (boundaryGeometryMapMap.find(entitySubdomain) == boundaryGeometryMapMap.end())
            boundaryGeometryMapMap.insert(std::pair<size_t, std::shared_ptr<GeometryMapType>>(
                entitySubdomain, std::shared_ptr<GeometryMapType>(new GeometryMapType())));
          GeometryMapType& boundaryGeometryMap = *(boundaryGeometryMapMap[entitySubdomain]);
          if (boundaryCodimSizesMap.find(entitySubdomain) == boundaryCodimSizesMap.end())
            boundaryCodimSizesMap.insert(std::pair<size_t, CodimSizesType>(entitySubdomain, CodimSizesType(0)));
          CodimSizesType& boundaryCodimSizes = boundaryCodimSizesMap[entitySubdomain];
          //   * and add geometry and global index of this codim 0 entity
          const GeometryType& entityGeometryType = entity.type();
          addGeometryAndIndex(boundaryGeometryMap, boundaryCodimSizes, entityGeometryType, entityGlobalIndex);
          //   * and of all remaining codims
          Add<1, dim>::subEntities(*this, entity, boundaryGeometryMap, boundaryCodimSizes);
          // for the intersection information of the boundary grid part
          //   * get the map for this subdomain (and create it, if necessary)
          if (boundaryInfoMap.find(entitySubdomain) == boundaryInfoMap.end())
            boundaryInfoMap.insert(std::pair<size_t, std::shared_ptr<EntityToIntersectionSetMapType>>(
                entitySubdomain,
                std::shared_ptr<EntityToIntersectionSetMapType>(new EntityToIntersectionSetMapType())));
          EntityToIntersectionSetMapType& boundaryInfo = *(boundaryInfoMap[entitySubdomain]);
          //   * and get the entry for this entity
          IntersectionInfoSetType& entityBoundaryInfo = boundaryInfo[entityGlobalIndex];
          //   * and add this local intersection
          entityBoundaryInfo.push_back(intersectionLocalIndex);
        } else if (intersection.neighbor()) {
          // then this entity lies inside the domain
          // and has a neighbor
          const auto neighbor = intersection.outside();
          const IndexType& neighborGlobalIndex = globalGridView_->indexSet().index(neighbor);
          const size_t neighborSubdomain = getSubdomainOf(neighborGlobalIndex);
          // check if neighbor is in another or in the same subdomain
          if (neighborSubdomain != entitySubdomain) {
            // get local index of the intersection
            const int intersectionLocalIndex = intersection.indexInInside();
            // for the neighbor information between the subdomains
            //   * the subdomain of the neighbor is a neighboring subdomain of the entities subdomain
            neighborsOfSubdomain.insert(neighborSubdomain);
            // for the subdomain grid part
            //   * get the boundary info map for this entity
            IntersectionToBoundaryIdMapType& entityInnerBoundaryInfo = subdomainInnerBoundaryInfo[entityGlobalIndex];
            //   * and add the local intersection id and its desired fake boundary id to this entities map
            entityInnerBoundaryInfo.insert(std::pair<int, size_t>(intersectionLocalIndex, boundary_segment_index_));
            // for the coupling grid part
            //   * get the coupling map for this subdomain
            SubdomainMapType& couplingMap = couplingMaps[entitySubdomain];
            //   * create an entry for the neighboring subdomain (if needed)
            if (couplingMap.find(neighborSubdomain) == couplingMap.end())
              couplingMap.insert(std::pair<size_t, std::shared_ptr<GeometryMapType>>(
                  neighborSubdomain, std::shared_ptr<GeometryMapType>(new GeometryMapType)));
            //   * and get it
            GeometryMapType& couplingGeometryMap = *(couplingMap[neighborSubdomain]);
            //   * and add geometry and global index of this codim 0 entity
            const GeometryType& entityGeometryType = entity.type();
            std::map<size_t, CodimSizesType>& couplingCodimSizeMap = couplingCodimSizeMaps[entitySubdomain];
            //   * create codim sizes vector for this neighbor (if needed)
            if (couplingCodimSizeMap.find(neighborSubdomain) == couplingCodimSizeMap.end())
              couplingCodimSizeMap.insert(std::pair<size_t, CodimSizesType>(neighborSubdomain, CodimSizesType(0)));
            //  * and get it
            CodimSizesType& couplingCodimSizes = couplingCodimSizeMap[neighborSubdomain];
            addGeometryAndIndex(couplingGeometryMap, couplingCodimSizes, entityGeometryType, entityGlobalIndex);
            //   * and of all remaining codims
            Add<1, dim>::subEntities(*this, entity, couplingGeometryMap, couplingCodimSizes);
            // for the intersection information of the coupling
            //   * get the map for this subdomain
            CouplingIntersectionMapType& couplingBoundaryInfo = couplingBoundaryInfos[entitySubdomain];
            //   * and create it for this neighbor (if needed)
            if (couplingBoundaryInfo.find(neighborSubdomain) == couplingBoundaryInfo.end())
              couplingBoundaryInfo.insert(std::pair<size_t, std::shared_ptr<EntityToIntersectionSetMapType>>(
                  neighborSubdomain,
                  std::shared_ptr<EntityToIntersectionSetMapType>(new EntityToIntersectionSetMapType())));
            //   * and get it
            EntityToIntersectionSetMapType& couplingBoundaryInfoMap = *(couplingBoundaryInfo[neighborSubdomain]);
            //   * get the entry for this entity
            IntersectionInfoSetType& entityCouplingBoundaryInfo = couplingBoundaryInfoMap[entityGlobalIndex];
            //   * and add this local intersection
            entityCouplingBoundaryInfo.push_back(intersectionLocalIndex);
          } else { // if neighbor is contained in this subdomain
            subdomainsEntitiesAreConnected = true;
          } // check if neighbor is in another subdomain
        } // check the type of this intersection
      } // walk the neighbors
      // check if this entity is connected to the other entities of this subdomain
      if (assert_connected && subdomainSizes[entitySubdomain] != 1 && !subdomainsEntitiesAreConnected) {
        DUNE_THROW(InvalidStateException,
                   "at least one entity of subdomain " << entitySubdomain << " is not connected to entity "
                                                       << entityGlobalIndex
                                                       << " (connected)!");
      } // check if this entity is connected to the other entities of this subdomain
    } // walk the global grid part
    // walk the subdomains
    //   * to create the local grid parts
    localGridParts_ =
        std::make_shared<std::vector<std::shared_ptr<const typename LocalGridViewType::Implementation>>>(size_);
    auto& localGridParts = *localGridParts_;
    for (typename SubdomainMapType::const_iterator subdomainIterator = subdomainToEntityMap_.begin();
         subdomainIterator != subdomainToEntityMap_.end();
         ++subdomainIterator) {
      const size_t subdomain = subdomainIterator->first;
      // for the local grid part
      //   * get the geometry map
      const std::shared_ptr<const GeometryMapType> localGeometryMap = subdomainIterator->second;
      //   * get the boundary info map
      const std::shared_ptr<const EntityToIntersectionInfoMapType> localBoundaryInfo =
          subdomainInnerBoundaryInfos[subdomain];
      //   * and create the local grid part
      localGridParts[subdomain] = std::make_shared<const typename LocalGridViewType::Implementation>(
          globalGridView_, localGeometryMap, localBoundaryInfo);
    } // walk the subdomains

    // walk those subdomains which have a boundary grid part
    //   * to create the boundary grid parts
    boundaryGridParts_ =
        std::make_shared<std::map<size_t, std::shared_ptr<const typename BoundaryGridViewType::Implementation>>>();
    auto& boundaryGridParts = *boundaryGridParts_;
    typename std::map<size_t, CodimSizesType>::const_iterator boundaryCodimSizesMapIt = boundaryCodimSizesMap.begin();
    typename std::map<size_t, std::shared_ptr<EntityToIntersectionSetMapType>>::const_iterator boundaryInfoMapIt =
        boundaryInfoMap.begin();
    for (typename std::map<size_t, std::shared_ptr<GeometryMapType>>::const_iterator
             boundaryGeometryMapMapIt = boundaryGeometryMapMap.begin();
         boundaryGeometryMapMapIt != boundaryGeometryMapMap.end();
         ++boundaryGeometryMapMapIt, ++boundaryCodimSizesMapIt, ++boundaryInfoMapIt) {
      const size_t boundarySubdomain = boundaryGeometryMapMapIt->first;
      assert(boundarySubdomain == boundaryCodimSizesMapIt->first
             && "We should not get here: we are in big trouble, if these maps do not correspond to each other!");
      assert(boundarySubdomain == boundaryInfoMapIt->first
             && "We should not get here: we are in big trouble, if these maps do not correspond to each other!");
      // for the boundary grid part
      //   * get the geometry map
      const std::shared_ptr<const GeometryMapType> boundaryGeometryMap = boundaryGeometryMapMapIt->second;
      //   * get the boundary info map
      const std::shared_ptr<const EntityToIntersectionSetMapType> boundaryBoundaryInfo = boundaryInfoMapIt->second;
      //   * and create the boundary grid part
      boundaryGridParts.emplace(
          boundarySubdomain,
          std::make_shared<const typename BoundaryGridViewType::Implementation>(
              globalGridView_, boundaryGeometryMap, boundaryBoundaryInfo, localGridParts[boundarySubdomain]));
    } // walk those subdomains which have a boundary grid part
    // walk the subdomains
    //   * to create the coupling grid parts
    couplingGridPartsMaps_ =
        std::make_shared<std::vector<std::map<size_t,
                                              std::shared_ptr<const typename CouplingGridViewType::Implementation>>>>(
            size_);
    auto& couplingGridPartsMaps = *couplingGridPartsMaps_;
    for (size_t subdomain = 0; subdomain < couplingMaps.size(); ++subdomain) {
      // get the coupling map for this subdomain
      const SubdomainMapType& couplingMap = couplingMaps[subdomain];
      // get the intersection information map for this subdomain
      const CouplingIntersectionMapType& couplingBoundaryInfo = couplingBoundaryInfos[subdomain];
      // get the target map for this subdomain
      auto& couplingGridPartsMap = couplingGridPartsMaps[subdomain];
      // loop over all neighbors
      for (typename SubdomainMapType::const_iterator neighborIt = couplingMap.begin(); neighborIt != couplingMap.end();
           ++neighborIt) {
        const size_t neighbor = neighborIt->first;
        // get the geometry map
        const std::shared_ptr<const GeometryMapType> couplingGeometryMap = neighborIt->second;
        // get the boundary info map
        typename CouplingIntersectionMapType::const_iterator result = couplingBoundaryInfo.find(neighbor);
        assert(result != couplingBoundaryInfo.end() && "This should not happen (see above)!");
        const std::shared_ptr<const EntityToIntersectionSetMapType> coupling_boundary_info = result->second;
        // and create the coupling grid part
        couplingGridPartsMap.emplace(
            neighbor,
            std::make_shared<const typename CouplingGridViewType::Implementation>(globalGridView_,
                                                                                  couplingGeometryMap,
                                                                                  coupling_boundary_info,
                                                                                  localGridParts[subdomain],
                                                                                  localGridParts[neighbor]));
      } // loop over all neighbors
    } // walk the subdomains

    // create the first layer of oversampling
    if (oversamplingLayers > 0) {
      oversampledLocalGridParts_ = addOneLayerOfOverSampling(
          subdomainToEntityMap_, *localGridParts_, neighbor_recursion_level, subdomainToOversamplingEntitiesMap_);
      oversampled_ = true;
    }
    // and the rest
    for (size_t ii = 1; ii < oversamplingLayers; ++ii) {
      auto tmpOversapledGridParts = oversampledLocalGridParts_;
      oversampledLocalGridParts_ =
          std::make_shared<std::vector<std::shared_ptr<const typename LocalGridViewType::Implementation>>>(size_);
      SubdomainMapType tmpSubdomainToOversamplingEntitiesMap = subdomainToOversamplingEntitiesMap_;
      subdomainToOversamplingEntitiesMap_ = SubdomainMapType();
      oversampledLocalGridParts_ = addOneLayerOfOverSampling(tmpSubdomainToOversamplingEntitiesMap,
                                                             *tmpOversapledGridParts,
                                                             neighbor_recursion_level,
                                                             subdomainToOversamplingEntitiesMap_);
    }

    // done
    finalized_ = true;
  } // ... finalize(...)

  std::shared_ptr<DdGridType> createMsGrid()
  {
    assert(finalized_ && "Please call finalize() before calling createMsGrid()!");
    if (oversampled_)
      return std::make_shared<DdGridType>(grid_,
                                          globalGridView_,
                                          size_,
                                          neighboringSubdomainSets_,
                                          entityToSubdomainMap_,
                                          localGridParts_,
                                          boundaryGridParts_,
                                          couplingGridPartsMaps_,
                                          oversampledLocalGridParts_);
    else
      return std::make_shared<DdGridType>(grid_,
                                          globalGridView_,
                                          size_,
                                          neighboringSubdomainSets_,
                                          entityToSubdomainMap_,
                                          localGridParts_,
                                          boundaryGridParts_,
                                          couplingGridPartsMaps_);
  } // ... createMsGrid(...)

private:
  void addGeometryAndIndex(GeometryMapType& geometryMap,
                           CodimSizesType& localCodimSizes,
                           const GeometryType& geometryType,
                           const IndexType& globalIndex)
  {
    // get the map to this geometry type
    typename GeometryMapType::mapped_type& indexMap = geometryMap[geometryType];
    // add if needed
    if (indexMap.find(globalIndex) == indexMap.end()) {
      const size_t codim = dim - geometryType.dim();
      const IndexType localIndex = boost::numeric_cast<IndexType>(localCodimSizes[codim]);
      indexMap.insert(std::pair<IndexType, IndexType>(globalIndex, localIndex));
      // increase count for this codim
      ++(localCodimSizes[codim]);
    }
  } // ... addGeometryAndIndex(...)

  size_t getSubdomainOf(const IndexType& globalIndex) const
  {
    const typename EntityToSubdomainMapType::const_iterator result = entityToSubdomainMap_->find(globalIndex);
    if (result == entityToSubdomainMap_->end())
      DUNE_THROW(InvalidStateException, "entity " << globalIndex << " not added to any subdomain!");
    return result->second;
  } // ... getSubdomainOf(...)

  std::shared_ptr<std::vector<std::shared_ptr<const typename LocalGridViewType::Implementation>>>
  addOneLayerOfOverSampling(
      const SubdomainMapType& subdomainToEntityMap,
      const std::vector<std::shared_ptr<const typename LocalGridViewType::Implementation>>& localGridParts,
      const size_t neighbor_recursion_level,
      SubdomainMapType& subdomainToOversamplingEntitiesMap)
  {
    // init data structures
    // for the subdomains inner boundaries
    //   * to map the local intersection index to the desired fake boundary id
    typedef std::map<int, size_t> IntersectionToBoundaryIdMapType;
    //   * to map the global entity index to one of those maps
    typedef std::map<IndexType, IntersectionToBoundaryIdMapType> EntityToIntersectionInfoMapType;
    //   * to hold one of those maps for each subdomain
    std::vector<std::shared_ptr<EntityToIntersectionInfoMapType>> oversamplingSubdomainInnerBoundaryInfos(size_);
    // walk the subdomains to create the oversampling
    for (auto subdomainIt = subdomainToEntityMap.begin(); subdomainIt != subdomainToEntityMap.end(); ++subdomainIt) {
      const size_t subdomain = subdomainIt->first;
      const GeometryMapType& geometryMap = *(subdomainIt->second);
      // * therefore, hardcopy the existing map we want to extend,
      std::shared_ptr<GeometryMapType> geometryMapCopy = std::make_shared<GeometryMapType>(geometryMap);
      // * and create an empty local boundary info map for later use
      oversamplingSubdomainInnerBoundaryInfos[subdomain] = std::make_shared<EntityToIntersectionInfoMapType>();
      // * then walk the local grid part to find the local boundary entities
      const LocalGridViewType& localGridPart = *(localGridParts[subdomain]);
      for (auto entityIt = localGridPart.template begin<0>(); entityIt != localGridPart.template end<0>(); ++entityIt) {
        // get the entity index
        const EntityType& entity = *entityIt;
        //        const IndexType entityGlobalIndex = globalGridView_->indexSet().index(entity);
        // lets see if this is a boundary entity of the local grid part
        bool isOnLocalBoundary = false;
        for (auto intersectionIt = localGridPart.ibegin(entity); intersectionIt != localGridPart.iend(entity);
             ++intersectionIt) {
          const auto& intersection = *intersectionIt;
          if (intersection.boundary())
            isOnLocalBoundary = true;
        }
        if (isOnLocalBoundary) {
          // add all the "neighbors"
          // * therefore, iterate over the intersections in the global grid part
          for (auto intersectionIt = globalGridView_->ibegin(entity); intersectionIt != globalGridView_->iend(entity);
               ++intersectionIt) {
            const auto& intersection = *intersectionIt;
            // if this intersection is not on the domain boundary
            if (intersection.neighbor()) {
              // get the neighbor
              const auto neighbor = intersection.outside();
              const IndexType& neighborGlobalIndex = globalGridView_->indexSet().index(neighbor);
              //              const size_t neighborSubdomain = getSubdomainOf(neighborGlobalIndex);
              // if the neighbor is not in the subdomain
              bool neighborIsNotInThisSubdomain = true;
              if (geometryMap.find(neighbor.type()) != geometryMap.end()) {
                if (geometryMap.find(neighbor.type())->second.find(neighborGlobalIndex)
                    != geometryMap.find(neighbor.type())->second.end()) {
                  neighborIsNotInThisSubdomain = false;
                }
              }
              if (neighborIsNotInThisSubdomain) {
                // add him to the oversampling
                // * therefore we can use the old localCodimSizes,
                CodimSizesType& localCodimSizes = localCodimSizes_.find(subdomain)->second;
                // * add the neighbor
                addGeometryAndIndex(*geometryMapCopy, localCodimSizes, neighbor.type(), neighborGlobalIndex);
                // * and all remaining codims entities
                Add<1, dim>::subEntities(*this, neighbor, *geometryMapCopy, localCodimSizes);
                // and also check all its neighbors
                if (neighbor_recursion_level > 0)
                  add_neighbors_neighbors_recursively(
                      entity, neighbor, geometryMap, neighbor_recursion_level, localCodimSizes, *geometryMapCopy);
              } // if the neighbor is not in the subdomain
            } // if this intersection is not on the domain boundary
          } // iterate over the intersections in the global grid part
        } // lets see if this is a boundary entity of the local grid part
      } // then walk the local grid part to find the local boundary entities
      subdomainToOversamplingEntitiesMap.insert(std::make_pair(subdomain, geometryMapCopy));
    } // walk the subdomains to create the oversampling

    // now we need to create the local boundary info for the oversampling, so walk the global grid part
    for (auto entityIt = globalGridView_->template begin<0>(); entityIt != globalGridView_->template end<0>();
         ++entityIt) {
      const auto& entity = *entityIt;
      const IndexType entityIndex = globalGridView_->indexSet().index(entity);
      // now we find all the oversampled subdomains this entity is a part of
      for (auto subdomainToOversamplingEntitiesMapIt : subdomainToOversamplingEntitiesMap) {
        const auto geometryMapIt = subdomainToOversamplingEntitiesMapIt.second->find(entity.type());
        if (geometryMapIt != subdomainToOversamplingEntitiesMapIt.second->end()) {
          const auto& geometryMap = geometryMapIt->second;
          if (geometryMap.find(entityIndex) != geometryMap.end()) {
            // this entity is a part of this subdomain!
            const size_t entitySubdomain = subdomainToOversamplingEntitiesMapIt.first;
            // then walk the neighbors
            for (auto intersectionIt = globalGridView_->ibegin(entity); intersectionIt != globalGridView_->iend(entity);
                 ++intersectionIt) {
              const auto& intersection = *intersectionIt;
              if (intersection.neighbor()) {
                const auto neighbor = intersection.outside();
                const IndexType neighborIndex = globalGridView_->indexSet().index(neighbor);
                // and check, if the neighbor is in the same subdomain
                bool isInSame = false;
                const auto neighborGeometryMapIt = subdomainToOversamplingEntitiesMapIt.second->find(neighbor.type());
                if (neighborGeometryMapIt != subdomainToOversamplingEntitiesMapIt.second->end()) {
                  const auto& neighborGeometryMap = neighborGeometryMapIt->second;
                  if (neighborGeometryMap.find(neighborIndex) != neighborGeometryMap.end()) {
                    isInSame = true;
                  }
                } // and check, if the neighbor is in the same subdomain
                if (!isInSame) {
                  // the neighbor is not part of this oversampled subdomain,
                  // so the entity in question is on the boundary!
                  const int intersectionLocalIndex = intersection.indexInInside();
                  auto& localBoundaryInfo = *(oversamplingSubdomainInnerBoundaryInfos[entitySubdomain]);
                  // get the boundary info map for this entity
                  IntersectionToBoundaryIdMapType& entityBoundaryInfo = localBoundaryInfo[entityIndex];
                  // and add the local intersection id and its desired fake boundary id to this entities map
                  entityBoundaryInfo.insert(std::pair<int, size_t>(intersectionLocalIndex, boundary_segment_index_));
                } // if (!isInSame)
              }
            } // then walk the neighbors
          }
        }
      } // now we find all the oversampled subdomains this entity is a part of
    } // walk the global grid part

    // and create the oversampled local grid parts
    auto oversampledLocalGridPartsRet =
        std::make_shared<std::vector<std::shared_ptr<const typename LocalGridViewType::Implementation>>>(size_);
    auto& oversampledLocalGridParts = *oversampledLocalGridPartsRet;
    for (auto subdomainIterator = subdomainToOversamplingEntitiesMap.begin();
         subdomainIterator != subdomainToOversamplingEntitiesMap.end();
         ++subdomainIterator) {
      // report
      const size_t subdomain = subdomainIterator->first;
      // for the local grid part
      //   * get the geometry map
      const std::shared_ptr<const GeometryMapType> localGeometryMap = subdomainIterator->second;
      //   * get the boundary info map
      const std::shared_ptr<const EntityToIntersectionInfoMapType> localBoundaryInfo =
          oversamplingSubdomainInnerBoundaryInfos[subdomain];
      //   * and create the local grid part
      oversampledLocalGridParts[subdomain] = std::make_shared<const typename LocalGridViewType::Implementation>(
          globalGridView_, localGeometryMap, localBoundaryInfo);
    } // and crete the oversampled local grid parts
    return oversampledLocalGridPartsRet;
  } // ... addOneLayerOfOverSampling(...)

  template <class EntityType, class NeighborType>
  void add_neighbors_neighbors_recursively(const EntityType& entity,
                                           const NeighborType& neighbor,
                                           const GeometryMapType& geometryMap,
                                           size_t recursion_level,
                                           CodimSizesType& localCodimSizes,
                                           GeometryMapType& geometryMapCopy)
  {
    // loop over all the neighbors of the neighbor
    for (auto neighborIntersectionIt = globalGridView_->ibegin(neighbor);
         neighborIntersectionIt != globalGridView_->iend(neighbor);
         ++neighborIntersectionIt) {
      const auto& neighborIntersection = *neighborIntersectionIt;
      if (neighborIntersection.neighbor()) {
        // get the neighbors neighbor
        const auto neighborsNeighbor = neighborIntersection.outside();
        const IndexType neighborsNeighborGlobalIndex = globalGridView_->indexSet().index(neighborsNeighbor);
        bool neighborsNeighborIsNotInThisSubdomain = true;
        if (geometryMap.find(neighborsNeighbor.type()) != geometryMap.end()) {
          if (geometryMap.find(neighborsNeighbor.type())->second.find(neighborsNeighborGlobalIndex)
              != geometryMap.find(neighborsNeighbor.type())->second.end()) {
            neighborsNeighborIsNotInThisSubdomain = false;
          }
        }
        // if the neighbor is not in the subdomain
        if (neighborsNeighborIsNotInThisSubdomain) {
          // check, if he intersects the entity
          const auto& entityGeometry = entity.geometry();
          const auto& neighborsNeighborGeometry = neighborsNeighbor.geometry();
          // * therefore loop over all corners of the entity
          for (int ii = 0; ii < entityGeometry.corners(); ++ii) {
            const auto entityCorner = entityGeometry.corner(ii);
            // then loop over all the corners of the neighbors neighbor
            for (int jj = 0; jj < neighborsNeighborGeometry.corners(); ++jj) {
              const auto neighborsNeighborCorner = neighborsNeighborGeometry.corner(jj);
              // and check for equality
              if (entityCorner == neighborsNeighborCorner) {
                // then add the neighbors neighbor
                addGeometryAndIndex(
                    geometryMapCopy, localCodimSizes, neighborsNeighbor.type(), neighborsNeighborGlobalIndex);
                Add<1, dim>::subEntities(*this, neighborsNeighbor, geometryMapCopy, localCodimSizes);
              }
            }
          }
        } // if the neighbor is not in the subdomain
        // call this function on the neighbours neighbor
        if (recursion_level > 0)
          add_neighbors_neighbors_recursively(
              entity, neighborsNeighbor, geometryMap, --recursion_level, localCodimSizes, geometryMapCopy);
      }
    } // loop over all the neighbors of the neighbor
  } // ... add_neighbors_neighbors_recursively(...)

  // friends
  template <int, int>
  friend struct Add;

  // members
  const std::shared_ptr<const GridType> grid_;
  const size_t boundary_segment_index_;
  bool prepared_;
  bool finalized_;
  size_t size_;
  std::shared_ptr<const GlobalGridViewType> globalGridView_;
  // for the entity <-> subdomain relations
  std::shared_ptr<EntityToSubdomainMapType> entityToSubdomainMap_;
  SubdomainMapType subdomainToEntityMap_;
  SubdomainMapType subdomainToOversamplingEntitiesMap_;
  // for the neighboring information
  std::shared_ptr<std::vector<NeighboringSubdomainsSetType>> neighboringSubdomainSets_;
  // for the local grid parts
  std::map<size_t, CodimSizesType> localCodimSizes_;
  std::shared_ptr<std::vector<std::shared_ptr<const typename LocalGridViewType::Implementation>>> localGridParts_;
  std::shared_ptr<std::vector<std::shared_ptr<const typename LocalGridViewType::Implementation>>>
      oversampledLocalGridParts_;
  // for the boundary grid parts
  std::shared_ptr<std::map<size_t, std::shared_ptr<const typename BoundaryGridViewType::Implementation>>>
      boundaryGridParts_;
  // for the coupling grid parts
  std::shared_ptr<std::vector<std::map<size_t, std::shared_ptr<const typename CouplingGridViewType::Implementation>>>>
      couplingGridPartsMaps_;
  bool oversampled_;
}; // class SubdomainGridFactory


//! specialization to stop the recursion
template <class GridType>
template <int c>
struct SubdomainGridFactory<GridType>::Add<c, c>
{
  static void subEntities(SubdomainGridFactory<GridType>& factory,
                          const typename SubdomainGridFactory<GridType>::EntityType& entity,
                          typename SubdomainGridFactory<GridType>::GeometryMapType& geometryMap,
                          typename SubdomainGridFactory<GridType>::CodimSizesType& localCodimSizes)
  {
    // loop over all codim c subentities of this entity
    for (unsigned int i = 0; i < entity.subEntities(c); ++i) {
      const auto codimCentity = entity.template subEntity<c>(i);
      const SubdomainGridFactory<GridType>::GeometryType& geometryType = codimCentity.type();
      const typename SubdomainGridFactory<GridType>::IndexType globalIndex =
          factory.globalGridView_->indexSet().index(codimCentity);
      factory.addGeometryAndIndex(geometryMap, localCodimSizes, geometryType, globalIndex);
    } // loop over all codim c subentities of this entity
  } // static void subEntities()
}; // struct SubdomainGridFactory< GridType >::Add< c, c >


} // namespace DD
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_DD_SUBDOMAINS_FACTORY_HH
