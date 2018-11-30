// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2017 - 2018)
//   Tobias Leibner  (2017 - 2018)

#ifndef DUNE_XT_GRID_VIEW_SUBDOMAIN_INDEXSET_HH
#define DUNE_XT_GRID_VIEW_SUBDOMAIN_INDEXSET_HH

#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/geometry/type.hh>
#include <dune/xt/grid/type_traits.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace internal {


/**
 *  \brief      Given a Dune::IndexSet and a set of entity indices, provides an index set on those entities only.
 *  \todo       Replace GlobalGridViewImp by Interface!
 *  \todo       Document!
 */
template <class GlobalGridViewImp>
class IndexBasedIndexSet : public Dune::IndexSet<GlobalGridViewImp,
                                                 IndexBasedIndexSet<GlobalGridViewImp>,
                                                 typename extract_index_set_t<GlobalGridViewImp>::IndexType,
                                                 typename extract_index_set_t<GlobalGridViewImp>::Types>
{
public:
  typedef GlobalGridViewImp GlobalGridViewType;
  typedef IndexBasedIndexSet<GlobalGridViewType> ThisType;
  using BaseType = Dune::IndexSet<GlobalGridViewImp,
                                  IndexBasedIndexSet<GlobalGridViewImp>,
                                  typename extract_index_set_t<GlobalGridViewImp>::IndexType,
                                  typename extract_index_set_t<GlobalGridViewImp>::Types>;
  typedef typename GlobalGridViewType::Grid GridType;
  typedef Dune::GeometryType GeometryType;
  typedef typename BaseType::IndexType IndexType;
  static const unsigned int dimension = GridType::dimension;
  typedef std::map<GeometryType, std::map<IndexType, IndexType>> IndexContainerType;

private:
  typedef std::map<IndexType, IndexType> Indices_MapType;

public:
  IndexBasedIndexSet(const GlobalGridViewType& globalGridView,
                     const std::shared_ptr<const IndexContainerType> indexContainer)
    : BaseType()
    , indexContainer_(indexContainer)
    , sizeByCodim_(dimension + 1, IndexType(0))
    , geometryTypesByCodim_(dimension + 1)
    , globalGridView_(globalGridView)
  {
    // get geometry types and compute sizes
    for (typename IndexContainerType::const_iterator iterator = indexContainer_->begin();
         iterator != indexContainer_->end();
         ++iterator) {
      const GeometryType& geometryType = iterator->first;
      const IndexType sz = boost::numeric_cast<IndexType>(iterator->second.size());
      const unsigned int codim = dimension - geometryType.dim();
      geometryTypesByCodim_[codim].push_back(geometryType);
      sizeByGeometryType_.insert(std::pair<GeometryType, IndexType>(geometryType, sz));
      assert(codim <= dimension);
      sizeByCodim_[codim] += sz;
    }
  } // IndexBasedIndexSet(const GridPartType& localGridPart)

  template <int cc>
  IndexType index(const typename GridType::template Codim<cc>::Entity& entity) const
  {
    return getIndex(entity);
  }

  template <class EntityType>
  IndexType index(const EntityType& entity) const
  {
    return getIndex(entity);
  }

  template <int cc>
  IndexType subIndex(const typename GridType::template Codim<cc>::Entity& entity, int i, unsigned int codim) const
  {
    // get the global subindex
    const IndexType& globalSubIndex = globalGridView_.indexSet().template subIndex<cc>(entity, i, codim);
    const int subCodim = cc + codim;
    assert(0 <= subCodim && subCodim <= int(dimension) && "This should not happen, we have a bad codimension");
    const unsigned int subDim = dimension - subCodim;
    // loop over all geometry types
    for (typename IndexContainerType::const_iterator it = indexContainer_->begin(); it != indexContainer_->end();
         ++it) {
      // check if geometry type has correct dimension
      const GeometryType& geometryType = it->first;
      if (geometryType.dim() == subDim) {
        // get the index map
        const std::map<IndexType, IndexType>& indexMap = it->second;
        // search for the global index
        const typename std::map<IndexType, IndexType>::const_iterator result = indexMap.find(globalSubIndex);
        if (result != indexMap.end()) {
          // return the corresponding local index
          const IndexType localSubIndex = result->second;
          return localSubIndex;
        } // search for the global index
      } // check if geometry type has correct dimension
    } // loop over all geometry types
    // if we came this far we did not find it
    std::stringstream msg;
    msg << std::endl
        << "Error in subIndex< " << cc << " >(" << entity.type() << ", " << i << ", " << codim << "):" << std::endl
        << "  subIndex in the global index set is " << globalSubIndex << ", local subIndex could not be found!"
        << std::endl;
    DUNE_THROW(Dune::InvalidStateException, msg.str());
  } // IndexType subIndex(const typename GridType::template Codim< cc >::Entity& entity, int i, unsigned int codim)
  // const

  template <class EntityType>
  IndexType subIndex(const EntityType& entity, int i, unsigned int codim) const
  {
    // get the global subindex
    const IndexType& globalSubIndex = globalGridView_.indexSet().subIndex(entity, i, codim);
    const int subCodim = EntityType::codimension + codim;
    assert(0 <= subCodim && subCodim <= int(dimension) && "This should not happen, we have a bad codimension");
    const unsigned int subDim = dimension - subCodim;
    // loop over all geometry types
    for (typename IndexContainerType::const_iterator it = indexContainer_->begin(); it != indexContainer_->end();
         ++it) {
      // check if geometry type has correct dimension
      const GeometryType& geometryType = it->first;
      if (geometryType.dim() == subDim) {
        // get the index map
        const std::map<IndexType, IndexType>& indexMap = it->second;
        // search for the global index
        const typename std::map<IndexType, IndexType>::const_iterator result = indexMap.find(globalSubIndex);
        if (result != indexMap.end()) {
          // return the corresponding local index
          const IndexType localSubIndex = result->second;
          return localSubIndex;
        } // search for the global index
      } // check if geometry type has correct dimension
    } // loop over all geometry types
    // if we came this far we did not find it
    std::stringstream msg;
    msg << std::endl
        << "Error in subIndex(" << entity.type() << ", " << i << ", " << codim << "):" << std::endl
        << "  subIndex in the global index set is " << globalSubIndex << ", local subIndex could not be found!"
        << std::endl;
    DUNE_THROW(Dune::InvalidStateException, msg.str());
  } // IndexType subIndex(const EntityType& entity, int i, unsigned int codim) const

  const std::vector<GeometryType>& geomTypes(int codim) const
  {
    return geometryTypesByCodim_[codim];
  }

  //! note here the container is returned by value opposed to a reference in geomTypes !
  const std::vector<GeometryType> types(int codim) const
  {
    return geomTypes(codim);
  }


  IndexType size(GeometryType type) const
  {
    assert(sizeByGeometryType_.find(type) != sizeByGeometryType_.end());
    return sizeByGeometryType_.find(type)->second;
  }

  IndexType size(int codim) const
  {
    assert(0 <= codim);
    assert(codim <= int(dimension));
    return sizeByCodim_[codim];
  }

  template <class EntityType>
  bool contains(const EntityType& entity) const
  {
    // check if we have an index map for this GeometryType
    const GeometryType& geometryType = entity.type();
    const typename IndexContainerType::const_iterator indexMap = indexContainer_->find(geometryType);
    if (indexMap != indexContainer_->end()) {
      // check if this entity is listen in the map
      const IndexType globalIndex = globalGridView_.indexSet().index(entity);
      if (indexMap->second.find(globalIndex) != indexMap->second.end())
        return true;
      else
        return false;
    } else
      return false;
  } // bool contains(const EntityType& entity) const

private:
  //! \todo This method can be optimized by the use of assert()s instead of if()s and throw()s!
  template <class EntityType>
  IndexType getIndex(const EntityType& entity) const
  {
    // check if we have an index map for this GeometryType
    const GeometryType& geometryType = entity.type();
    const typename IndexContainerType::const_iterator indexMap = indexContainer_->find(geometryType);
    if (indexMap != indexContainer_->end()) {
      // check if this entity is listen in the map
      const IndexType globalIndex = globalGridView_.indexSet().index(entity);
      const typename Indices_MapType::const_iterator indexPair = indexMap->second.find(globalIndex);
      if (indexPair != indexMap->second.end())
        return indexPair->second;
      else
        DUNE_THROW(Dune::InvalidStateException, "Given entity not contained in index set!");
    } else
      DUNE_THROW(Dune::InvalidStateException, "Given geometry not contained in index set!");
  } // IndexType getIndex(const EntityType& entity) const

  IndexType findLocalIndex(const IndexType& globalIndex) const
  {
    // find the map with this global index and return the local one
    for (typename IndexContainerType::const_iterator it = indexContainer_->begin(); it != indexContainer_->end();
         ++it) {
      const Indices_MapType& indicesMap = it->second;
      const typename Indices_MapType::const_iterator indexPair = indicesMap.find(globalIndex);
      if (indexPair != indicesMap.end())
        return indexPair->second;
      else {
        std::cerr << std::endl
                  << "Error in indexContainer_[" << it->first << "]: globalIndex " << globalIndex << " not found!"
                  << std::endl;
        std::cerr << "  ";
        for (typename Indices_MapType::const_iterator indexIt = indicesMap.begin(); indexIt != indicesMap.end();
             ++indexIt) {
          std::cerr << "(" << indexIt->first << "/" << indexIt->second << "), ";
        }
        std::cerr << std::endl;
      }
    } // find the map with this global index and return the local one
    // we should not come this far
    DUNE_THROW(Dune::InvalidStateException, "Given entity not contained in index set!");
  } // IndexType findLocalIndex(const IndexType& globalIndex) const

  const Dune::shared_ptr<const IndexContainerType> indexContainer_;
  std::vector<IndexType> sizeByCodim_;
  std::vector<std::vector<GeometryType>> geometryTypesByCodim_;
  std::map<GeometryType, IndexType> sizeByGeometryType_;
  const GlobalGridViewType globalGridView_;
}; // class IndexBasedIndexSet


} // namespace internal
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_VIEW_SUBDOMAIN_INDEXSET_HH
