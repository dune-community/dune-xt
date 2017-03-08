// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_GRID_DD_SUBDOMAINS_GRID_HH
#define DUNE_XT_GRID_DD_SUBDOMAINS_GRID_HH

#include <vector>
#include <set>
#include <map>
#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/leafgridpart.hh>
#endif

#include <dune/xt/common/color.hh>
#include <dune/xt/common/type_traits.hh>

#include <dune/xt/grid/view/subdomain/part.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace DD {

#if HAVE_DUNE_FEM


/**
 *  \brief      Implementation of domain decomposition by subdomain grid parts.
 *
 *              The basic assumption is, that each entity of the global grid part is added to at most one local grid
 *              part. Overlapping or the like is then handled afterwards by adding additional overlapping regions to
 *              an existing local grid part (not yet implemented).
 *  \attention  Works only for one GeometryType per Codim (I think)! Problem is, that indices are unique per
 * GeometryType, not per Codim.
 *  \todo       Resolve the above Issue (should be easy, compute the size as a sum over all GeometryTypes for a given
 * codim)!
 *
 *  \todo       Giving the local and global gridparts as shared pointers is quite misleading, since it is not
 *              guaranteed that the underlying grid parts and the grid will exist forever. So we should change those to
 *              reference imho.
 */
template <class GridImp>
class SubdomainGrid
{
public:
  typedef GridImp GridType;
  typedef SubdomainGrid<GridType> ThisType;
  static const unsigned int dim = GridType::dimension;
  static const unsigned int dimension = GridType::dimension;
  typedef typename GridType::ctype ctype;
  typedef Fem::LeafGridPart<GridType> GlobalGridPartType;
  typedef SubdomainGridPart<GlobalGridPartType> LocalGridPartType;
  typedef SubdomainBoundaryGridPart<GlobalGridPartType> BoundaryGridPartType;
  typedef SubdomainCouplingGridPart<GlobalGridPartType> CouplingGridPartType;
  typedef typename GlobalGridPartType::template Codim<0>::EntityType EntityType;
  typedef typename GlobalGridPartType::IndexSetType::IndexType IndexType;
  typedef std::set<size_t> NeighborSetType;
  //! map type which maps from an entity index (of the global grid parts index set) to a subdomain
  typedef std::map<IndexType, size_t> EntityToSubdomainMapType;

  SubdomainGrid(
      const std::shared_ptr<const GridType> grd,
      const std::shared_ptr<const GlobalGridPartType> globalGrdPrt,
      const size_t sz,
      const std::shared_ptr<const std::vector<NeighborSetType>> neighboringSets,
      const std::shared_ptr<const EntityToSubdomainMapType> entityToSubdMap,
      const std::shared_ptr<const std::vector<std::shared_ptr<const LocalGridPartType>>> localGridParts,
      const std::shared_ptr<const std::map<size_t, std::shared_ptr<const BoundaryGridPartType>>> boundaryGridParts,
      const std::shared_ptr<const std::vector<std::map<size_t, std::shared_ptr<const CouplingGridPartType>>>>
          couplingGridPartsMaps)
    : grid_(grd)
    , globalGridPart_(globalGrdPrt)
    , size_(sz)
    , neighboringSetsPtr_(neighboringSets)
    , entityToSubdomainMap_(entityToSubdMap)
    , localGridParts_(localGridParts)
    , boundaryGridParts_(boundaryGridParts)
    , couplingGridPartsMaps_(couplingGridPartsMaps)
    , oversampling_(false)
  {
    // check for correct sizes
    std::stringstream msg;
    bool error = false;
    msg << "Error:" << std::endl;
    if (localGridParts_->size() != size_) {
      msg << "  - 'localGridParts' has wrong size (is " << localGridParts_->size() << ", should be " << size_ << ")!"
          << std::endl;
      error = true;
    }
    if (couplingGridPartsMaps_->size() != size_) {
      msg << "  - 'couplingGridPartsMaps' has wrong size (is " << couplingGridPartsMaps_->size() << ", should be "
          << size_ << ")!" << std::endl;
      error = true;
    }
    if (error)
      DUNE_THROW(InvalidStateException, msg.str());
  } // SubdomainGrid()

  SubdomainGrid(
      const std::shared_ptr<const GridType> grd,
      const std::shared_ptr<const GlobalGridPartType> globalGrdPrt,
      const size_t sz,
      const std::shared_ptr<const std::vector<NeighborSetType>> neighboringSets,
      const std::shared_ptr<const EntityToSubdomainMapType> entityToSubdMap,
      const std::shared_ptr<const std::vector<std::shared_ptr<const LocalGridPartType>>> localGridParts,
      const std::shared_ptr<const std::map<size_t, std::shared_ptr<const BoundaryGridPartType>>> boundaryGridParts,
      const std::shared_ptr<const std::vector<std::map<size_t, std::shared_ptr<const CouplingGridPartType>>>>
          couplingGridPartsMaps,
      const std::shared_ptr<const std::vector<std::shared_ptr<const LocalGridPartType>>> oversampledLocalGridParts)
    : grid_(grd)
    , globalGridPart_(globalGrdPrt)
    , size_(sz)
    , neighboringSetsPtr_(neighboringSets)
    , entityToSubdomainMap_(entityToSubdMap)
    , localGridParts_(localGridParts)
    , boundaryGridParts_(boundaryGridParts)
    , couplingGridPartsMaps_(couplingGridPartsMaps)
    , oversampling_(true)
    , oversampledLocalGridParts_(oversampledLocalGridParts)
  {
    // check for correct sizes
    std::stringstream msg;
    bool error = false;
    msg << "Error:" << std::endl;
    if (localGridParts_->size() != size_) {
      msg << "  - 'localGridParts' has wrong size (is " << localGridParts_->size() << ", should be " << size_ << ")!"
          << std::endl;
      error = true;
    }
    if (couplingGridPartsMaps_->size() != size_) {
      msg << "  - 'couplingGridPartsMaps' has wrong size (is " << couplingGridPartsMaps_->size() << ", should be "
          << size_ << ")!" << std::endl;
      error = true;
    }
    if (error)
      DUNE_THROW(InvalidStateException, msg.str());
  } // SubdomainGrid()

  SubdomainGrid(const ThisType& other) = default;

  SubdomainGrid(ThisType&& source) = default;

  const std::shared_ptr<const GridType>& grid() const
  {
    return grid_;
  }

  const GlobalGridPartType& globalGridPart() const
  {
    return *globalGridPart_;
  }

  size_t size() const
  {
    return size_;
  }

  bool oversampling() const
  {
    return oversampling_;
  }

  LocalGridPartType localGridPart(const size_t subdomain, const bool ovrsmplng = false) const
  {
    assert(subdomain < size_);
    if (!ovrsmplng) {
      const std::vector<std::shared_ptr<const LocalGridPartType>>& localGridParts = *localGridParts_;
      return *(localGridParts[subdomain]);
    } else {
      if (!oversampling_)
        DUNE_THROW(InvalidStateException,
                   "\n" << Common::color_string_red("ERROR:")
                        << " oversampled local gridpart requested from a grid without oversampling!");
      const std::vector<std::shared_ptr<const LocalGridPartType>>& oversampledLocalGridParts =
          *oversampledLocalGridParts_;
      return *(oversampledLocalGridParts[subdomain]);
    }
  } // ... localGridPart(...)

  bool boundary(const size_t subdomain) const
  {
    assert(subdomain < size_);
    const std::map<size_t, std::shared_ptr<const BoundaryGridPartType>>& boundaryGridParts = *boundaryGridParts_;
    return (boundaryGridParts.find(subdomain) != boundaryGridParts.end());
  }

  BoundaryGridPartType boundaryGridPart(const size_t subdomain) const
  {
    assert(subdomain < size_);
    const std::map<size_t, std::shared_ptr<const BoundaryGridPartType>>& boundaryGridParts = *boundaryGridParts_;
    typename std::map<size_t, std::shared_ptr<const BoundaryGridPartType>>::const_iterator result =
        boundaryGridParts.find(subdomain);
    assert(result != boundaryGridParts.end()
           && "Only call boundaryGridPart(subdomain), if boundary(subdomain) is true!");
    return *(result->second);
  }

  CouplingGridPartType couplingGridPart(const size_t subdomain, const size_t neighbor) const
  {
    assert(subdomain < size_);
    assert(neighbor < size_);
    const std::vector<std::map<size_t, std::shared_ptr<const CouplingGridPartType>>>& couplingGridPartsMaps =
        *couplingGridPartsMaps_;
    const std::map<size_t, std::shared_ptr<const CouplingGridPartType>>& couplingGridPartsMap =
        couplingGridPartsMaps[subdomain];
    const typename std::map<size_t, std::shared_ptr<const CouplingGridPartType>>::const_iterator result =
        couplingGridPartsMap.find(neighbor);
    if (result == couplingGridPartsMap.end()) {
      std::stringstream msg;
      msg << "subdomain " << neighbor << " is not a neighbor of subdomain " << subdomain << "!";
      DUNE_THROW(InvalidStateException, msg.str());
    }
    return *(result->second);
  } // ... couplingGridPart(...)

  const std::shared_ptr<const EntityToSubdomainMapType>& entityToSubdomainMap() const
  {
    return entityToSubdomainMap_;
  }

  const NeighborSetType& neighborsOf(const size_t subdomain) const
  {
    assert(subdomain < size_);
    const std::vector<NeighborSetType>& neighboringSets = *neighboringSetsPtr_;
    return neighboringSets[subdomain];
  }
  size_t subdomainOf(const IndexType& globalIndex) const
  {
    const typename EntityToSubdomainMapType::const_iterator result = entityToSubdomainMap_->find(globalIndex);
    if (result == entityToSubdomainMap_->end()) {
      std::stringstream msg;
      msg << "missing information for entity " << globalIndex << "in entityToSubdomainMap_!";
      DUNE_THROW(InvalidStateException, msg.str());
    }
    return result->second;
  } // ... getSubdomainOf(...)

  size_t subdomainOf(const EntityType& entity) const
  {
    return subdomainOf(globalGridPart_->indexSet().index(entity));
  }

private:
  const std::shared_ptr<const GridType> grid_;
  const std::shared_ptr<const GlobalGridPartType> globalGridPart_;
  const size_t size_;
  const std::shared_ptr<const std::vector<NeighborSetType>> neighboringSetsPtr_;
  const std::shared_ptr<const EntityToSubdomainMapType> entityToSubdomainMap_;
  const std::shared_ptr<const std::vector<std::shared_ptr<const LocalGridPartType>>> localGridParts_;
  const std::shared_ptr<const std::map<size_t, std::shared_ptr<const BoundaryGridPartType>>> boundaryGridParts_;
  const std::shared_ptr<const std::vector<std::map<size_t, std::shared_ptr<const CouplingGridPartType>>>>
      couplingGridPartsMaps_;
  bool oversampling_;
  const std::shared_ptr<const std::vector<std::shared_ptr<const LocalGridPartType>>> oversampledLocalGridParts_;
}; // class SubdomainGrid


#else // HAVE_DUNE_FEM


template <class GridImp>
class SubdomainGrid
{
  static_assert(AlwaysFalse<GridImp>::value, "You are missing dune-fem!");
};


#endif // HAVE_DUNE_FEM

} // namespace DD
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_DD_SUBDOMAINS_GRID_HH
