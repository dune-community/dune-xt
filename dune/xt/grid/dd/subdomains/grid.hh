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
#include <dune/grid/common/gridview.hh>

#include <dune/xt/common/color.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/view/subdomain/view.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace DD {

/**
 *  \brief      Implementation of domain decomposition by subdomain grid views.
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
  typedef typename GridType::LeafGridView GlobalGridViewType;
  typedef GridView<XT::Grid::internal::SubdomainGridViewTraits<GlobalGridViewType>> LocalGridViewType;
  typedef GridView<XT::Grid::internal::SubdomainBoundaryGridViewTraits<GlobalGridViewType>> BoundaryGridViewType;
  typedef GridView<XT::Grid::internal::SubdomainCouplingGridViewTraits<GlobalGridViewType>> CouplingGridViewType;
  using EntityType = XT::Grid::extract_entity_t<GlobalGridViewType>;
  typedef typename GlobalGridViewType::IndexSet::IndexType IndexType;
  typedef std::set<size_t> NeighborSetType;
  //! map type which maps from an entity index (of the global grid parts index set) to a subdomain
  typedef std::map<IndexType, size_t> EntityToSubdomainMapType;

  SubdomainGrid(
      const std::shared_ptr<const GridType> grd,
      const std::shared_ptr<const GlobalGridViewType> globalGrdPrt,
      const size_t sz,
      const std::shared_ptr<const std::vector<NeighborSetType>> neighboringSets,
      const std::shared_ptr<const EntityToSubdomainMapType> entityToSubdMap,
      const std::shared_ptr<const std::vector<std::shared_ptr<const typename LocalGridViewType::Implementation>>>
          localGridViews,
      const std::shared_ptr<const std::map<size_t,
                                           std::shared_ptr<const typename BoundaryGridViewType::Implementation>>>
          boundaryGridViews,
      const std::
          shared_ptr<const std::vector<std::map<size_t,
                                                std::shared_ptr<const typename CouplingGridViewType::Implementation>>>>
              couplingGridViewsMaps)
    : grid_(grd)
    , globalGridView_(globalGrdPrt)
    , size_(sz)
    , neighboringSetsPtr_(neighboringSets)
    , entityToSubdomainMap_(entityToSubdMap)
    , localGridViews_(localGridViews)
    , boundaryGridViews_(boundaryGridViews)
    , couplingGridViewsMaps_(couplingGridViewsMaps)
    , oversampling_(false)
    , oversampledLocalGridViews_(nullptr)
  {
    // check for correct sizes
    std::stringstream msg;
    bool error = false;
    msg << "Error:" << std::endl;
    if (localGridViews_->size() != size_) {
      msg << "  - 'localGridViews' has wrong size (is " << localGridViews_->size() << ", should be " << size_ << ")!"
          << std::endl;
      error = true;
    }
    if (couplingGridViewsMaps_->size() != size_) {
      msg << "  - 'couplingGridViewsMaps' has wrong size (is " << couplingGridViewsMaps_->size() << ", should be "
          << size_ << ")!" << std::endl;
      error = true;
    }
    if (error)
      DUNE_THROW(InvalidStateException, msg.str());
  } // SubdomainGrid()

  SubdomainGrid(
      const std::shared_ptr<const GridType> grd,
      const std::shared_ptr<const GlobalGridViewType> globalGrdPrt,
      const size_t sz,
      const std::shared_ptr<const std::vector<NeighborSetType>> neighboringSets,
      const std::shared_ptr<const EntityToSubdomainMapType> entityToSubdMap,
      const std::shared_ptr<const std::vector<std::shared_ptr<const typename LocalGridViewType::Implementation>>>
          localGridViews,
      const std::shared_ptr<const std::map<size_t,
                                           std::shared_ptr<const typename BoundaryGridViewType::Implementation>>>
          boundaryGridViews,
      const std::
          shared_ptr<const std::vector<std::map<size_t,
                                                std::shared_ptr<const typename CouplingGridViewType::Implementation>>>>
              couplingGridViewsMaps,
      const std::shared_ptr<const std::vector<std::shared_ptr<const typename LocalGridViewType::Implementation>>>
          oversampledLocalGridViews)
    : grid_(grd)
    , globalGridView_(globalGrdPrt)
    , size_(sz)
    , neighboringSetsPtr_(neighboringSets)
    , entityToSubdomainMap_(entityToSubdMap)
    , localGridViews_(localGridViews)
    , boundaryGridViews_(boundaryGridViews)
    , couplingGridViewsMaps_(couplingGridViewsMaps)
    , oversampling_(true)
    , oversampledLocalGridViews_(oversampledLocalGridViews)
  {
    // check for correct sizes
    std::stringstream msg;
    bool error = false;
    msg << "Error:" << std::endl;
    if (localGridViews_->size() != size_) {
      msg << "  - 'localGridViews' has wrong size (is " << localGridViews_->size() << ", should be " << size_ << ")!"
          << std::endl;
      error = true;
    }
    if (oversampledLocalGridViews_->size() != size_) {
      msg << "  - 'oversampledLocalGridViews' has wrong size (is " << oversampledLocalGridViews_->size()
          << ", should be " << size_ << ")!" << std::endl;
      error = true;
    }
    if (couplingGridViewsMaps_->size() != size_) {
      msg << "  - 'couplingGridViewsMaps' has wrong size (is " << couplingGridViewsMaps_->size() << ", should be "
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

  const GlobalGridViewType& globalGridView() const
  {
    return *globalGridView_;
  }

  const GlobalGridViewType& global_grid_view() const
  {
    return *globalGridView_;
  }

  size_t size() const
  {
    return size_;
  }

  bool oversampling() const
  {
    return oversampling_;
  }

  LocalGridViewType local_grid_view(const size_t subdomain, const bool ovrsmplng = false) const
  {
    assert(subdomain < size_);
    if (!ovrsmplng)
      return *((*localGridViews_)[subdomain]);
    else {
      if (!oversampling_)
        DUNE_THROW(InvalidStateException,
                   "\n" << Common::color_string_red("ERROR:")
                        << " oversampled local grid part requested from a grid without oversampling!");
      return *((*oversampledLocalGridViews_)[subdomain]);
    }
  } // ... local_grid_view(...)

  bool boundary(const size_t subdomain) const
  {
    assert(subdomain < size_);
    const auto& boundaryGridViews = *boundaryGridViews_;
    return (boundaryGridViews.find(subdomain) != boundaryGridViews.end());
  }

  BoundaryGridViewType boundary_grid_view(const size_t subdomain) const
  {
    // See for instance test_dd_subdomains_cube_2d, when adapting
    // boundary_parts_contain_only_boundary_entities_and_intersections to use these views!
    assert(subdomain < size_);
    const auto result = boundaryGridViews_->find(subdomain);
    assert(result != boundaryGridViews_->end()
           && "Only call boundary_grid_view(subdomain), if boundary(subdomain) is true!");
    return *(result->second);
  }

  CouplingGridViewType coupling_grid_view(const size_t subdomain, const size_t neighbor) const
  {
    // See for instance test_dd_subdomains_cube_2d, when adapting
    // coupling_parts_contain_only_inner_entities_and_intersections to use these views!
    assert(subdomain < size_);
    assert(neighbor < size_);
    const auto& couplingGridViews = (*couplingGridViewsMaps_)[subdomain];
    const auto result = couplingGridViews.find(neighbor);
    if (result == couplingGridViews.end()) {
      DUNE_THROW(InvalidStateException,
                 "subdomain " << neighbor << " is not a neighbor of subdomain " << subdomain << "!");
    }
    assert(result->second);
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
    return subdomainOf(globalGridView_->indexSet().index(entity));
  }

  const std::shared_ptr<const GridType> grid_;
  const std::shared_ptr<const GlobalGridViewType> globalGridView_;
  const size_t size_;
  const std::shared_ptr<const std::vector<NeighborSetType>> neighboringSetsPtr_;
  const std::shared_ptr<const EntityToSubdomainMapType> entityToSubdomainMap_;
  const std::shared_ptr<const std::vector<std::shared_ptr<const typename LocalGridViewType::Implementation>>>
      localGridViews_;
  const std::shared_ptr<const std::map<size_t, std::shared_ptr<const typename BoundaryGridViewType::Implementation>>>
      boundaryGridViews_;
  const std::
      shared_ptr<const std::vector<std::map<size_t,
                                            std::shared_ptr<const typename CouplingGridViewType::Implementation>>>>
          couplingGridViewsMaps_;
  bool oversampling_;
  const std::shared_ptr<const std::vector<std::shared_ptr<const typename LocalGridViewType::Implementation>>>
      oversampledLocalGridViews_;
}; // class SubdomainGrid


} // namespace DD
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_DD_SUBDOMAINS_GRID_HH
