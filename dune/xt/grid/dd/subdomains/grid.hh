// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

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
#include <dune/xt/grid/type_traits.hh>
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
  typedef typename GlobalGridPartType::GridViewType GlobalGridViewType;
  typedef SubdomainGridPart<GlobalGridPartType> LocalGridPartType;
  typedef Fem::GridPart2GridView<LocalGridPartType> LocalGridViewType;
  typedef SubdomainBoundaryGridPart<GlobalGridPartType> BoundaryGridPartType;
  typedef typename BoundaryGridPartType::GridViewType BoundaryGridViewType;
  typedef SubdomainCouplingGridPart<GlobalGridPartType> CouplingGridPartType;
  typedef typename CouplingGridPartType::GridViewType CouplingGridViewType;
  using EntityType = XT::Grid::extract_entity_t<GlobalGridPartType>;
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
    , global_grid_view_(std::make_shared<GlobalGridViewType>(*globalGridPart_))
    , size_(sz)
    , neighboringSetsPtr_(neighboringSets)
    , entityToSubdomainMap_(entityToSubdMap)
    , localGridParts_(localGridParts)
    , boundaryGridParts_(boundaryGridParts)
    , couplingGridPartsMaps_(couplingGridPartsMaps)
    , oversampling_(false)
    , oversampledLocalGridParts_(nullptr)
    , local_grid_views_(new std::vector<std::shared_ptr<const LocalGridViewType>>(localGridParts_->size()))
    , oversampled_local_grid_views_(nullptr)
    , boundary_grid_views_(new std::map<size_t, std::shared_ptr<const BoundaryGridViewType>>())
    , coupling_grid_views_maps_(new std::vector<std::map<size_t, std::shared_ptr<const CouplingGridViewType>>>(
          couplingGridPartsMaps_->size()))
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
    build_views();
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
    , global_grid_view_(std::make_shared<GlobalGridViewType>(*globalGridPart_))
    , size_(sz)
    , neighboringSetsPtr_(neighboringSets)
    , entityToSubdomainMap_(entityToSubdMap)
    , localGridParts_(localGridParts)
    , boundaryGridParts_(boundaryGridParts)
    , couplingGridPartsMaps_(couplingGridPartsMaps)
    , oversampling_(true)
    , oversampledLocalGridParts_(oversampledLocalGridParts)
    , local_grid_views_(new std::vector<std::shared_ptr<const LocalGridViewType>>(localGridParts_->size()))
    , oversampled_local_grid_views_(
          new std::vector<std::shared_ptr<const LocalGridViewType>>(oversampledLocalGridParts_->size()))
    , boundary_grid_views_(new std::map<size_t, std::shared_ptr<const BoundaryGridViewType>>())
    , coupling_grid_views_maps_(new std::vector<std::map<size_t, std::shared_ptr<const CouplingGridViewType>>>(
          couplingGridPartsMaps_->size()))
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
    if (oversampledLocalGridParts_->size() != size_) {
      msg << "  - 'oversampledLocalGridParts' has wrong size (is " << oversampledLocalGridParts_->size()
          << ", should be " << size_ << ")!" << std::endl;
      error = true;
    }
    if (couplingGridPartsMaps_->size() != size_) {
      msg << "  - 'couplingGridPartsMaps' has wrong size (is " << couplingGridPartsMaps_->size() << ", should be "
          << size_ << ")!" << std::endl;
      error = true;
    }
    if (error)
      DUNE_THROW(InvalidStateException, msg.str());
    build_views();
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

  const GlobalGridViewType& global_grid_view() const
  {
    return *global_grid_view_;
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
    if (!ovrsmplng)
      return *((*localGridParts_)[subdomain]);
    else {
      if (!oversampling_)
        DUNE_THROW(InvalidStateException,
                   "\n" << Common::color_string_red("ERROR:")
                        << " oversampled local grid part requested from a grid without oversampling!");
      return *((*oversampledLocalGridParts_)[subdomain]);
    }
  } // ... localGridPart(...)

  LocalGridViewType local_grid_view(const size_t subdomain, const bool ovrsmplng = false) const
  {
    assert(subdomain < size_);
    if (!ovrsmplng)
      return *((*local_grid_views_)[subdomain]);
    else {
      if (!oversampling_)
        DUNE_THROW(InvalidStateException,
                   "\n" << Common::color_string_red("ERROR:")
                        << " oversampled local grid part requested from a grid without oversampling!");
      return *((*oversampled_local_grid_views_)[subdomain]);
    }
  } // ... local_grid_view(...)

  bool boundary(const size_t subdomain) const
  {
    assert(subdomain < size_);
    const std::map<size_t, std::shared_ptr<const BoundaryGridPartType>>& boundaryGridParts = *boundaryGridParts_;
    return (boundaryGridParts.find(subdomain) != boundaryGridParts.end());
  }

  BoundaryGridPartType boundaryGridPart(const size_t subdomain) const
  {
    assert(subdomain < size_);
    const auto result = boundaryGridParts_->find(subdomain);
    assert(result != boundaryGridParts_->end()
           && "Only call boundaryGridPart(subdomain), if boundary(subdomain) is true!");
    return *(result->second);
  }

  BoundaryGridViewType boundary_grid_view(const size_t subdomain) const
  {
    // See for instance test_dd_subdomains_cube_2d, when adapting
    // boundary_parts_contain_only_boundary_entities_and_intersections to use these views!
    DUNE_THROW(InvalidStateException,
               "The intersections of the boundary grid view are currently broken, use boundaryGridPart() instead!");
    assert(subdomain < size_);
    const auto result = boundary_grid_views_->find(subdomain);
    assert(result != boundary_grid_views_->end()
           && "Only call boundary_grid_view(subdomain), if boundary(subdomain) is true!");
    return *(result->second);
  }

  CouplingGridPartType couplingGridPart(const size_t subdomain, const size_t neighbor) const
  {
    assert(subdomain < size_);
    assert(neighbor < size_);
    const auto& couplingGridPartsMap = (*couplingGridPartsMaps_)[subdomain];
    const auto result = couplingGridPartsMap.find(neighbor);
    if (result == couplingGridPartsMap.end()) {
      DUNE_THROW(InvalidStateException,
                 "subdomain " << neighbor << " is not a neighbor of subdomain " << subdomain << "!");
    }
    return *(result->second);
  } // ... couplingGridPart(...)

  CouplingGridViewType coupling_grid_view(const size_t subdomain, const size_t neighbor) const
  {
    // See for instance test_dd_subdomains_cube_2d, when adapting
    // coupling_parts_contain_only_inner_entities_and_intersections to use these views!
    DUNE_THROW(InvalidStateException,
               "The intersections of the coupling grid view are currently broken, use couplingGridPart() instead!");
    assert(subdomain < size_);
    assert(neighbor < size_);
    const auto& coupling_grid_views_map = (*coupling_grid_views_maps_)[subdomain];
    const auto result = coupling_grid_views_map.find(neighbor);
    if (result == coupling_grid_views_map.end()) {
      DUNE_THROW(InvalidStateException,
                 "subdomain " << neighbor << " is not a neighbor of subdomain " << subdomain << "!");
    }
    return *(result->second);
  }

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
  void build_views()
  {
    assert(local_grid_views_->size() == localGridParts_->size());
    for (size_t ii = 0; ii < localGridParts_->size(); ++ii)
      (*local_grid_views_)[ii] = std::make_shared<LocalGridViewType>(*((*localGridParts_)[ii]));
    if (oversampling_) {
      assert(oversampled_local_grid_views_->size() == oversampledLocalGridParts_->size());
      for (size_t ii = 0; ii < oversampledLocalGridParts_->size(); ++ii)
        (*oversampled_local_grid_views_)[ii] =
            std::make_shared<LocalGridViewType>(*((*oversampledLocalGridParts_)[ii]));
    }
    for (const auto& key_value_pair : *boundaryGridParts_)
      boundary_grid_views_->insert(
          std::make_pair(key_value_pair.first, std::make_shared<BoundaryGridViewType>(*key_value_pair.second)));
    assert(coupling_grid_views_maps_->size() == couplingGridPartsMaps_->size());
    for (size_t ii = 0; ii < couplingGridPartsMaps_->size(); ++ii) {
      (*coupling_grid_views_maps_)[ii] = std::map<size_t, std::shared_ptr<const CouplingGridViewType>>();
      for (const auto& key_value_pair : (*couplingGridPartsMaps_)[ii])
        (*coupling_grid_views_maps_)[ii].insert(
            std::make_pair(key_value_pair.first, std::make_shared<CouplingGridViewType>(*key_value_pair.second)));
    }
  } // ... build_views(...)

  const std::shared_ptr<const GridType> grid_;
  const std::shared_ptr<const GlobalGridPartType> globalGridPart_;
  const std::shared_ptr<const GlobalGridViewType> global_grid_view_;
  const size_t size_;
  const std::shared_ptr<const std::vector<NeighborSetType>> neighboringSetsPtr_;
  const std::shared_ptr<const EntityToSubdomainMapType> entityToSubdomainMap_;
  const std::shared_ptr<const std::vector<std::shared_ptr<const LocalGridPartType>>> localGridParts_;
  const std::shared_ptr<const std::map<size_t, std::shared_ptr<const BoundaryGridPartType>>> boundaryGridParts_;
  const std::shared_ptr<const std::vector<std::map<size_t, std::shared_ptr<const CouplingGridPartType>>>>
      couplingGridPartsMaps_;
  bool oversampling_;
  const std::shared_ptr<const std::vector<std::shared_ptr<const LocalGridPartType>>> oversampledLocalGridParts_;
  std::shared_ptr<std::vector<std::shared_ptr<const LocalGridViewType>>> local_grid_views_;
  std::shared_ptr<std::vector<std::shared_ptr<const LocalGridViewType>>> oversampled_local_grid_views_;
  std::shared_ptr<std::map<size_t, std::shared_ptr<const BoundaryGridViewType>>> boundary_grid_views_;
  std::shared_ptr<std::vector<std::map<size_t, std::shared_ptr<const CouplingGridViewType>>>> coupling_grid_views_maps_;
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
