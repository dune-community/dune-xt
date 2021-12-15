// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018, 2020)
//   René Fritze     (2018 - 2020)
//   Tobias Leibner  (2017 - 2020)

#ifndef DUNE_XT_GRID_DD_GLUED_HH
#define DUNE_XT_GRID_DD_GLUED_HH

#include <array>
#include <algorithm>
#include <memory>
#include <map>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/grid/common/rangegenerators.hh>

#if HAVE_DUNE_GRID_GLUE
#  include <dune/grid-glue/extractors/codim1extractor.hh>
#  include <dune/grid-glue/gridglue.hh>
#  include <dune/grid-glue/merging/contactmerge.hh>
#endif // HAVE_DUNE_GRID_GLUE

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/functions/base/function-as-grid-function.hh>
#include <dune/xt/functions/base/visualization.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/grid-function.hh>
#include <dune/xt/functions/interfaces/function.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/search.hh>
#include <dune/xt/grid/type_traits.hh>
#include <utility>

namespace Dune::XT::Grid::DD {
namespace Exceptions {


class intersection_orientation_is_broken : public Dune::InvalidStateException
{};


} // namespace Exceptions


#if HAVE_DUNE_GRID_GLUE


template <typename P0, typename P1>
size_t check_for_broken_coupling_intersections(
    const GridGlue::GridGlue<P0, P1>& glue,
    const typename P0::ctype& tolerance = 10 * XT::Common::FloatCmp::DefaultEpsilon<typename P0::ctype>::value())
{
  const auto inside_grid_view = glue.template gridView<0>();
  size_t failures = 0;
  // walk the coupling
  const auto coupling_intersection_it_end = glue.template iend<0>();
  for (auto coupling_intersection_it = glue.template ibegin<0>();
       coupling_intersection_it != coupling_intersection_it_end;
       ++coupling_intersection_it) {
    const auto& coupling_intersection = *coupling_intersection_it;
    const auto coupling_intersection_normal = coupling_intersection.centerUnitOuterNormal();
    const auto local_entity = coupling_intersection.inside();
    typename std::remove_const<decltype(coupling_intersection_normal)>::type local_intersection_normal(0.);
    // find the intersection of the local inside entity that corresponds to the coupling intersection
    size_t found = 0;
    for (auto&& local_intersection : intersections(inside_grid_view, local_entity)) {
      // the coupling intersection may be smaller than the local intersection, so check if all of the corners of the
      // coupling intersection lie within this local intersection
      int corners_inside = 0;
      for (auto ii : XT::Common::value_range(coupling_intersection.geometry().corners()))
        if (XT::Grid::contains(local_intersection, coupling_intersection.geometry().corner(ii)))
          ++corners_inside;
      if (corners_inside == coupling_intersection.geometry().corners()) {
        // this is the one
        ++found;
        local_intersection_normal = local_intersection.centerUnitOuterNormal();
      }
    }
    if (found != 1)
      DUNE_THROW(InvalidStateException,
                 "This should not happen!\n"
                     << "There were " << found << " local intersections which contain the coupling intersection, "
                     << "and there must not be more than one!");
    // now the expected normal is local_intersection_normal
    // and we would like coupling_intersection_normal to point in the same direction
    // since they have unit length, they should be identical
    if ((local_intersection_normal - coupling_intersection_normal).infinity_norm() > tolerance)
      ++failures;
  }
  return failures;
} // ... check_for_broken_coupling_intersections(...)


// forward
template <class MacroGridType, class LocalGridType, Layers layer = Layers::level>
class GluedVTKWriter;


template <class MacroGridImp, class LocalGridImp, Layers layer = Layers::level>
class Glued
{
  template <class G, bool anything = true>
  struct allowed_macro_grid
  {
    static constexpr bool value = true;
  };

  template <class G, bool anything = true>
  struct allowed_local_grid
  {
    static constexpr bool value = true;
  };

#  if HAVE_DUNE_ALUGRID
  template <class Comm, bool anything>
  struct allowed_local_grid<ALUGrid<3, 3, simplex, conforming, Comm>, anything>
  {
    static constexpr bool value = false;
  };

  template <class Comm, bool anything>
  struct allowed_local_grid<ALUGrid<3, 3, cube, conforming, Comm>, anything>
  {
    static constexpr bool value = false;
  };
#  endif // HAVE_DUNE_ALUGRID

#  if HAVE_MPI && (HAVE_DUNE_UGGRID || HAVE_UG)
  // UGGrid does not support multiple parallel instances in parallel and we have no means yet to create multiple
  // sequential grids once MPI was found.
  template <bool anything>
  struct allowed_local_grid<UGGrid<2>, anything>
  {
    static constexpr bool value = false;
  };

  template <bool anything>
  struct allowed_local_grid<UGGrid<3>, anything>
  {
    static constexpr bool value = false;
  };
#  endif

  static_assert(allowed_macro_grid<MacroGridImp>::value,
                "This macro grid is known to fail, enable on your onw risk by disabling this check!");
  static_assert(allowed_local_grid<LocalGridImp>::value,
                "This local grid is known to fail, enable on your onw risk by disabling this check!");

public:
  using MacroGridType = MacroGridImp;
  using MacroGridProviderType = XT::Grid::GridProvider<MacroGridType>;
  using MacroGridViewType = typename MacroGridProviderType::LeafGridViewType;
  using MacroEntityType = typename MacroGridType::template Codim<0>::Entity;

  using LocalGridType = LocalGridImp;
  using LocalGridProviderType = XT::Grid::GridProvider<LocalGridType>;
  using MicroGridViewType = typename LocalGridType::LevelGridView;
  using MicroEntityType = XT::Grid::extract_entity_t<MicroGridViewType>;

  using ctype = typename MacroGridType::ctype;
  static constexpr size_t dimDomain = MacroGridType::dimension;
  static constexpr size_t dimWorld = MacroGridType::dimensionworld;


  using LocalViewType = typename Layer<LocalGridType, layer, Backends::view>::type;

private:
  using LocalExtractorType = GridGlue::Codim1Extractor<LocalViewType>;

public:
  using GlueType = GridGlue::GridGlue<LocalExtractorType, LocalExtractorType>;

private:
  template <class GridView, class MacroIntersectionType>
  class CouplingFaceDescriptor
  {
    using LocalEntityType = XT::Grid::extract_entity_t<GridView>;
    using ctype = typename GridView::ctype;

  public:
    CouplingFaceDescriptor(const MacroIntersectionType& macro_intersection)
      : macro_intersection_(macro_intersection)
    {}

    bool operator()(const LocalEntityType& element, unsigned int face) const
    {
      const auto local_intersection_ptr = element.template subEntity<1>(face);
#  if DUNE_VERSION_GTE(DUNE_GRID, 2, 4)
      const auto& local_intersection = local_intersection_ptr;
#  else
      const auto& local_intersection = *local_intersection_ptr;
#  endif
      const auto& local_intersection_geometry = local_intersection.geometry();
      // Check if all corners of the local intersection lie within the macro intersection.
      for (auto ii : XT::Common::value_range(local_intersection_geometry.corners()))
        if (!XT::Grid::contains(macro_intersection_, local_intersection_geometry.corner(ii)))
          return false;
      return true;
    } // ... contains(...)

  private:
    const MacroIntersectionType& macro_intersection_;
  }; // CouplingFaceDescriptor

  template <class GV, class MI>
  static CouplingFaceDescriptor<GV, MI> create_descriptor(const GV& /*gv*/, const MI& mi)
  {
    return CouplingFaceDescriptor<GV, MI>(mi);
  }

public:
  Glued(MacroGridProviderType& macro_grid_provider,
        const size_t num_local_refinements = 0,
        const bool prepare_glues = false,
        const bool allow_for_broken_orientation_of_coupling_intersections = false,
        const ctype& allowed_overlap = 10 * XT::Common::FloatCmp::DefaultEpsilon<ctype>::value())
    : macro_grid_(macro_grid_provider)
    , allowed_overlap_(allowed_overlap)
    , macro_leaf_view_(macro_grid_.leaf_view())
    , macro_leaf_view_size_(macro_leaf_view_.indexSet().size(0))
    , local_grids_(macro_leaf_view_.indexSet().size(0), nullptr)
    , glues_(macro_leaf_view_.indexSet().size(0))
  {
    setup_local_grids();
    if (num_local_refinements > 0)
      for (auto& local_grid_provider : local_grids_) {
        assert(local_grid_provider);
        local_grid_provider->grid().globalRefine(boost::numeric_cast<int>(num_local_refinements));
      }
    if (prepare_glues)
      setup_glues(allow_for_broken_orientation_of_coupling_intersections);
  } // Glued(...)

  Glued(MacroGridProviderType& macro_grid_provider,
        const std::array<unsigned int, dimDomain>& num_elements_per_subdomain,
        const size_t num_local_refinements = 0,
        const bool prepare_glues = false,
        const bool allow_for_broken_orientation_of_coupling_intersections = false,
        const ctype& allowed_overlap = 10 * XT::Common::FloatCmp::DefaultEpsilon<ctype>::value())
    : macro_grid_(macro_grid_provider)
    , allowed_overlap_(allowed_overlap)
    , macro_leaf_view_(macro_grid_.leaf_view())
    , macro_leaf_view_size_(macro_leaf_view_.indexSet().size(0))
    , local_grids_(macro_leaf_view_.indexSet().size(0), nullptr)
    , glues_(macro_leaf_view_.indexSet().size(0))
  {
    setup_local_grids(num_elements_per_subdomain);
    if (num_local_refinements > 0)
      for (auto& local_grid_provider : local_grids_) {
        assert(local_grid_provider);
        local_grid_provider->grid().globalRefine(boost::numeric_cast<int>(num_local_refinements));
      }
    if (prepare_glues)
      setup_glues(allow_for_broken_orientation_of_coupling_intersections);
  } // Glued(...)

  const MacroGridViewType& macro_grid_view() const
  {
    assert_macro_grid_state();
    return macro_leaf_view_;
  }

  size_t num_subdomains() const
  {
    assert_macro_grid_state();
    return macro_leaf_view_size_;
  }

  size_t subdomain(const MacroEntityType& macro_entity) const
  {
    assert_macro_grid_state();
    assert(macro_leaf_view_.indexSet().contains(macro_entity));
    return macro_leaf_view_.indexSet().index(macro_entity);
  }

  bool boundary(const MacroEntityType& macro_entity) const
  {
    assert_macro_grid_state();
    assert(macro_leaf_view_.indexSet().contains(macro_entity));
    return macro_entity.hasBoundaryIntersections();
  }

  MicroGridViewType global_grid_view()
  {
    auto logger = XT::Common::TimedLogger().get("grid-multiscale.glued.global_grid_view");
    logger.warn() << "Requiring inefficient access to global micro grid!" << std::endl;
    assert_macro_grid_state();
    prepare_global_grid();
    return global_grid_->level_view(global_grid_->grid().maxLevel());
  }

  const std::vector<std::vector<size_t>>& local_to_global_indices()
  {
    auto logger = XT::Common::TimedLogger().get("grid-multiscale.glued.local_to_global_indices");
    logger.warn() << "Requiring inefficient access to global micro grid!" << std::endl;
    assert_macro_grid_state();
    prepare_global_grid();
    return *local_to_global_indices_;
  }

  MicroEntityType local_to_global_entity(const size_t subd, const MicroEntityType& local_entity)
  {
    return local_to_global_entity(
        subd, extract_local_view<layer>()(*local_grids_[subd], max_local_level(subd)).indexSet().index(local_entity));
  }

  MicroEntityType local_to_global_entity(const size_t subd, const size_t local_entity_index)
  {
    auto logger = XT::Common::TimedLogger().get("grid-multiscale.glued.local_to_global_entity");
    logger.warn() << "Requiring inefficient access to global micro grid!" << std::endl;
    assert_macro_grid_state();
    prepare_global_grid();
    const size_t global_index_of_local_entity = local_to_global_indices_->operator[](subd)[local_entity_index];
    const auto global_micro_grid_view = global_grid_view();
    for (auto&& entity : elements(global_micro_grid_view)) {
      if (global_micro_grid_view.indexSet().index(entity) == global_index_of_local_entity)
        return entity;
    }
    DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
               "subdomain: " << subd << "\n"
                             << "local_entity_index: " << local_entity_index << "\n"
                             << "global_index_of_local_entity: " << global_index_of_local_entity);
    return *global_micro_grid_view.template begin<0>();
  } // ... local_to_global_entity(...)

  MicroEntityType global_to_local_entity(const MicroEntityType& micro_entity)
  {
    assert_macro_grid_state();
    prepare_global_grid();
    return global_to_local_entity(global_grid_view().indexSet().index(micro_entity));
  }

  MicroEntityType global_to_local_entity(const size_t micro_entity_index)
  {
    auto logger = XT::Common::TimedLogger().get("grid-multiscale.glued.global_to_local_entity");
    logger.warn() << "Requiring inefficient access to global micro grid!" << std::endl;
    assert_macro_grid_state();
    prepare_global_grid();
    auto subdomain_and_local_entity_index = global_to_local_indices_->operator[](micro_entity_index);
    const auto subd = subdomain_and_local_entity_index.first;
    const auto local_index_of_global_entity = subdomain_and_local_entity_index.second;
    const auto local_grid_view = extract_local_view<layer>()(*local_grids_[subd], max_local_level(subd));
    for (auto&& entity : elements(local_grid_view)) {
      if (local_grid_view.indexSet().index(entity) == local_index_of_global_entity)
        return entity;
    }
    DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
               "micro_entity_index: " << micro_entity_index
                                      << "\n"
                                         "subdomain: "
                                      << subd << "\n"
                                      << "local_index_of_global_entity: " << local_index_of_global_entity);
    return *local_grid_view.template begin<0>();
  } // ... global_to_local_entity(...)

  const std::vector<std::pair<size_t, size_t>>& global_to_local_indices()
  {
    auto logger = XT::Common::TimedLogger().get("grid-multiscale.glued.global_to_local_indices");
    logger.warn() << "Requiring inefficient access to global micro grid!" << std::endl;
    assert_macro_grid_state();
    prepare_global_grid();
    assert(global_to_local_indices_);
    return *global_to_local_indices_;
  }

  const LocalGridProviderType& local_grid(const MacroEntityType& macro_entity) const
  {
    assert_macro_grid_state();
    assert(macro_leaf_view_.indexSet().contains(macro_entity));
    return local_grid(macro_leaf_view_.indexSet().index(macro_entity));
  }

  LocalGridProviderType& local_grid(const MacroEntityType& macro_entity)
  {
    assert_macro_grid_state();
    assert(macro_leaf_view_.indexSet().contains(macro_entity));
    return local_grid(macro_leaf_view_.indexSet().index(macro_entity));
  }

  LocalGridProviderType& local_grid(const size_t macro_entity_index)
  {
    assert_macro_grid_state();
    assert(macro_entity_index < macro_leaf_view_size_);
    return *(local_grids_.at(macro_entity_index));
  }

  const LocalGridProviderType& local_grid(const size_t macro_entity_index) const
  {
    assert_macro_grid_state();
    assert(macro_entity_index < macro_leaf_view_size_);
    return *(local_grids_.at(macro_entity_index));
  }

  /**
   * \brief Returns (and creates, if it does not exist) the coupling glue between the local grid view of level
   *        local_level_macro_entity on macro_entity and the local grid view of level local_level_macro_neighbor on
   * macro_neighbor.
   * \note  Access is not implemented efficiently. This could be improved by using std::map::find instead of
   *        std::map::operator[].
   */
  const GlueType& coupling(const MacroEntityType& macro_entity,
                           const int local_level_macro_entity,
                           const MacroEntityType& macro_neighbor,
                           const int local_level_macro_neighbor,
                           const bool allow_for_broken_orientation_of_coupling_intersections = false)
  {
    assert_macro_grid_state();
    const auto& macro_index_set = macro_leaf_view_.indexSet();
    assert(macro_index_set.contains(macro_entity));
    assert(macro_index_set.contains(macro_neighbor));
    if (local_level_macro_entity > max_local_level(macro_entity))
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                 "max_local_level(macro_entity): " << max_local_level(macro_entity) << "\n"
                                                   << "   local_level_macro_entity:      " << local_level_macro_entity);
    if (local_level_macro_neighbor > max_local_level(macro_neighbor))
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                 "max_local_level(macro_neighbor): " << max_local_level(macro_neighbor) << "\n"
                                                     << "   local_level_macro_neighbor:      "
                                                     << local_level_macro_neighbor);
    int local_entity_grid_level = local_level_macro_entity;
    int local_neighbor_grid_level = local_level_macro_neighbor;
    // in case of local level views, the local_level_macro... have to be non-negative
    if (layer == Layers::level) {
      if (local_level_macro_entity < 0)
        DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                   "local_level_macro_entity has to be non-negative (is " << local_level_macro_entity << ")!");
      if (local_level_macro_neighbor < 0)
        DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                   "local_level_macro_neighbor has to be non-negative (is " << local_level_macro_neighbor << ")!");
    } else if (layer == Layers::leaf) {
      // in case of local leaf views, it has to be -1
      if (local_level_macro_entity != -1)
        DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                   "local_level_macro_entity has to be -1 (is " << local_level_macro_entity << ")!");
      if (local_level_macro_neighbor != -1)
        DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                   "local_level_macro_neighbor has to be -1 (is " << local_level_macro_neighbor << ")!");
      // since the local leaf views may have been adapted, we use their respective sizes as keys for the glue map
      local_entity_grid_level *= local_grid(macro_entity).leaf_view().indexSet().size(0);
      local_neighbor_grid_level *= local_grid(macro_neighbor).leaf_view().indexSet().size(0);
    }
    const auto entity_index = macro_index_set.index(macro_entity);
    const auto neighbor_index = macro_index_set.index(macro_neighbor);
    auto& glues_for_this_entity_neighbor = glues_[entity_index][neighbor_index];
    if (glues_for_this_entity_neighbor[local_entity_grid_level][local_neighbor_grid_level] == nullptr) {
      if (layer == Layers::leaf) {
        // if the glue for the current leaf sizes does not exist, but for other leaf sizes, then those are not valid any
        // more and we can get rid of them
        if (glues_for_this_entity_neighbor.size() > 0)
          glues_for_this_entity_neighbor.clear();
      }
      // find the corresponding macro intersection ...
      for (auto&& macro_intersection : intersections(macro_leaf_view_, macro_entity)) {
        if (macro_intersection.neighbor() && !macro_intersection.boundary()) {
          const auto real_neighbor_ptr = macro_intersection.outside();
#  if DUNE_VERSION_GTE(DUNE_GRID, 2, 4)
          const auto& real_neighbor = real_neighbor_ptr;
#  else
          const auto& real_neighbor = *real_neighbor_ptr;
#  endif
          if (macro_index_set.index(real_neighbor) == neighbor_index)
            glues_for_this_entity_neighbor[local_entity_grid_level][local_neighbor_grid_level] =
                create_glue(macro_entity,
                            macro_neighbor,
                            macro_intersection,
                            (layer == Layers::level) ? local_entity_grid_level : -1,
                            (layer == Layers::level) ? local_neighbor_grid_level : -1);
        }
      } // ... find the corresponding macro intersection
    }

    const auto& glue = *(glues_for_this_entity_neighbor[local_entity_grid_level][local_neighbor_grid_level]);
    if (!allow_for_broken_orientation_of_coupling_intersections) {
      const size_t brocken_intersections = check_for_broken_coupling_intersections(glue);
      if (brocken_intersections > 0)
        DUNE_THROW(Exceptions::intersection_orientation_is_broken,
                   "The coupling glue between the grid views of\n"
                       << "     level " << local_entity_grid_level << " on macro entity   "
                       << macro_leaf_view_.indexSet().index(macro_entity) << " and\n"
                       << "     level " << local_neighbor_grid_level << " on macro neighbor "
                       << macro_leaf_view_.indexSet().index(macro_neighbor) << "\n"
                       << "   contains\n"
                       << "     " << brocken_intersections << "/" << glue.size()
                       << " intersections with wrong orientation!");
    }
    return glue;
  } // ... coupling(...)

  const GlueType& coupling(const MacroEntityType& macro_entity,
                           const MacroEntityType& macro_neighbor,
                           const bool allow_for_broken_orientation_of_coupling_intersections = false)
  {
    return coupling(macro_entity,
                    max_local_level(macro_entity),
                    macro_neighbor,
                    max_local_level(macro_neighbor),
                    allow_for_broken_orientation_of_coupling_intersections);
  }

  int max_local_level(const MacroEntityType& macro_entity) const
  {
    assert_macro_grid_state();
    if (layer == Layers::leaf)
      return -1;
    return local_grid(macro_entity).grid().maxLevel();
  }

  int max_local_level(const size_t macro_entity_index) const
  {
    assert_macro_grid_state();
    if (layer == Layers::leaf)
      return -1;
    return local_grid(macro_entity_index).grid().maxLevel();
  }

  const std::vector<std::pair<MicroEntityType, std::vector<int>>>&
  local_boundary_entities(const MacroEntityType& macro_entity, const int local_level)
  {
    assert_macro_grid_state();
    assert(macro_leaf_view_.indexSet().contains(macro_entity));
    const size_t macro_entity_index = macro_leaf_view_.indexSet().index(macro_entity);
    if (local_level > max_local_level(macro_entity))
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                 "macro_entity_index: " << macro_entity_index << "\n"
                                        << "local_level: " << local_level << "\n"
                                        << "max_local_level(macro_entity): " << max_local_level(macro_entity));
    auto& local_level_to_boundary_entity_ptrs_with_local_intersections =
        macro_entity_to_local_level_to_boundary_entity_ptrs_with_local_intersections_[macro_entity_index];
    auto& boundary_entity_ptrs_with_local_intersections =
        local_level_to_boundary_entity_ptrs_with_local_intersections[local_level];
    //    logger.debug() << "macro_entity: " << macro_entity_index << std::endl;
    if (boundary(macro_entity) && boundary_entity_ptrs_with_local_intersections.empty()) {
      // create the container, therefore
      const auto local_leaf_view = local_grids_[macro_entity_index]->leaf_view();
      // * walk the local grid (manually, to have access to the entity pointer)
      for (auto&& local_entity : elements(local_leaf_view)) {
        //        logger.debug() << "local_entity: " << local_leaf_view.indexSet().index(local_entity) << " ";
        if (local_entity.hasBoundaryIntersections()) {
          //          logger.debug() << "(boundary entity)" << std::endl;
          std::vector<int> local_boundary_intersections;
          // This entity has intersections on the local grid boundary, those could either be the domain boundary (which
          // we are looking for) or a boundary to another local grid (which we are not looking for). To find out
          // * walk the intersections
          for (auto&& local_intersection : intersections(local_leaf_view, local_entity)) {
            if (local_intersection.boundary() && !local_intersection.neighbor()) {
              //              logger.debug() << "local_intersection: " << local_intersection.indexInInside() << ":" <<
              //              std::endl;
              const auto local_intersection_geometry = local_intersection.geometry();
              const auto num_corners = boost::numeric_cast<size_t>(local_intersection_geometry.corners());
              // ** Check if all corners of the intersection lie on the domain boundary (aka the boundary intersection
              // of
              //    the macro entity this local grid belongs to. Therefore
              //    *** walk the intersections of the macro entity
              for (auto&& macro_intersection : intersections(macro_leaf_view_, macro_entity)) {
                if (macro_intersection.boundary() && !macro_intersection.neighbor()) {
                  // This macro intersection lies on the domain boundary, check if the local intersection is contained.
                  size_t corners_lie_on_boundary = 0;
                  for (size_t ii = 0; ii < num_corners; ++ii) {
                    if (XT::Grid::contains(macro_intersection,
                                           local_intersection_geometry.corner(boost::numeric_cast<int>(ii))))
                      ++corners_lie_on_boundary;
                  }
                  if (corners_lie_on_boundary == num_corners) {
                    // add the information to the container
                    local_boundary_intersections.push_back(local_intersection.indexInInside());
                  } // add the information to the container
                }
              } //    *** walk the intersections of the macro entity
            }
          } // * walk the intersections
          if (!local_boundary_intersections.empty()) {
            // add this local entity and its local intersections to the container
            boundary_entity_ptrs_with_local_intersections.emplace_back(local_entity, local_boundary_intersections);
          }
        } /*else
          logger.debug() << "(inner entity)" << std::endl;*/
      } // * walk the local grid
    } // create the container
    return boundary_entity_ptrs_with_local_intersections;
  } // ... local_boundary_entities(...)

  template <class... Args>
  void visualize(const std::string& filename = "grid.multiscale.glued", Args&&... args)
  {
    assert_macro_grid_state();
    auto logger = XT::Common::TimedLogger().get("grid-multiscale.glued.visualize");
    macro_grid_.visualize(filename + ".macro", std::forward<Args>(args)...);
    GluedVTKWriter<MacroGridType, LocalGridType, layer> vtk_writer(*this);
    const auto& macro_index_set = macro_leaf_view_.indexSet();
    std::vector<std::vector<double>> subdomain_visualization(macro_index_set.size(0));
    std::vector<std::vector<double>> boundary_visualization(macro_index_set.size(0));
    std::vector<std::vector<double>> inside_outside_coupling_visualization(macro_index_set.size(0));
    std::vector<std::vector<double>> outside_inside_coupling_visualization(macro_index_set.size(0));
    // walk the macro grid
    for (auto&& macro_entity : elements(macro_leaf_view_)) {
      const size_t macro_entity_index = macro_index_set.index(macro_entity);
      logger.debug() << "macro_entity: " << macro_entity_index << " ";
      const auto local_level = max_local_level(macro_entity);
      const auto local_grid_view = extract_local_view<layer>()(*local_grids_[macro_entity_index], local_level);
      subdomain_visualization[macro_entity_index] =
          std::vector<double>(local_grid_view.indexSet().size(0), macro_entity_index);
      boundary_visualization[macro_entity_index] = std::vector<double>(local_grid_view.indexSet().size(0), -1);
      if (inside_outside_coupling_visualization[macro_entity_index].empty())
        inside_outside_coupling_visualization[macro_entity_index] =
            std::vector<double>(local_grid_view.indexSet().size(0), -1);
      if (outside_inside_coupling_visualization[macro_entity_index].empty())
        outside_inside_coupling_visualization[macro_entity_index] =
            std::vector<double>(local_grid_view.indexSet().size(0), -1);
      // local boundary entities
      if (boundary(macro_entity)) {
        logger.debug() << "(boundary entity)" << std::endl;
        const auto& boundary_entities_with_local_intersections = local_boundary_entities(macro_entity, local_level);
        logger.debug() << "  " << boundary_entities_with_local_intersections.size() << "/"
                       << local_grid_view.indexSet().size(0) << " boundary entities" << std::endl;
        for (const auto& element : boundary_entities_with_local_intersections) {
          const auto& local_entity = element.first;
          const auto& local_intersections = element.second;
          if (!local_intersections.empty()) {
            const size_t local_entity_index = local_grid_view.indexSet().index(local_entity);
            boundary_visualization[macro_entity_index][local_entity_index] = static_cast<double>(macro_entity_index);
          }
        }
      } else
        logger.debug() << "(inner entity)" << std::endl;
      for (auto&& macro_intersection : intersections(macro_leaf_view_, macro_entity)) {
        if (!macro_intersection.boundary() && macro_intersection.neighbor()) {
          const auto macro_neighbor_ptr = macro_intersection.outside();
#  if DUNE_VERSION_GTE(DUNE_GRID, 2, 4)
          const auto& macro_neighbor = macro_neighbor_ptr;
#  else
          const auto& macro_neighbor = *macro_neighbor_ptr;
#  endif
          const size_t macro_neighbor_index = macro_leaf_view_.indexSet().index(macro_neighbor);
          const auto local_neighbor_level = max_local_level(macro_neighbor);
          const auto local_neighbor_grid_view =
              extract_local_view<layer>()(*local_grids_[macro_neighbor_index], local_neighbor_level);
          if (inside_outside_coupling_visualization[macro_neighbor_index].empty())
            inside_outside_coupling_visualization[macro_neighbor_index] =
                std::vector<double>(local_neighbor_grid_view.indexSet().size(0), -1);
          if (outside_inside_coupling_visualization[macro_neighbor_index].empty())
            outside_inside_coupling_visualization[macro_neighbor_index] =
                std::vector<double>(local_neighbor_grid_view.indexSet().size(0), -1);
          // walk the coupling, where this is the inside
          size_t num_coupling_intersections = 0;
          const auto& in_out_coupling_glue = coupling(macro_entity,
                                                      local_level,
                                                      macro_neighbor,
                                                      local_neighbor_level,
                                                      /*allow_for_broken_orientation_of_coupling_intersections=*/true);
          const auto in_out_coupling_intersection_it_end = in_out_coupling_glue.template iend<0>();
          for (auto in_out_coupling_intersection_it = in_out_coupling_glue.template ibegin<0>();
               in_out_coupling_intersection_it != in_out_coupling_intersection_it_end;
               ++in_out_coupling_intersection_it) {
            ++num_coupling_intersections;
            const auto& coupling_intersection = *in_out_coupling_intersection_it;
            const auto local_entity = coupling_intersection.inside();
            const size_t local_entity_index = local_grid_view.indexSet().index(local_entity);
            inside_outside_coupling_visualization[macro_entity_index][local_entity_index] =
                static_cast<double>(macro_entity_index);
            const auto local_neighbor = coupling_intersection.outside();
            const size_t local_neighbor_index = local_neighbor_grid_view.indexSet().index(local_neighbor);
            inside_outside_coupling_visualization[macro_neighbor_index][local_neighbor_index] =
                static_cast<double>(macro_neighbor_index);
          }
          // walk the coupling, where this is the outside
          size_t out_in_num_coupling_intersections = 0;
          const auto& out_in_coupling_glue = coupling(macro_neighbor,
                                                      local_neighbor_level,
                                                      macro_entity,
                                                      local_level,
                                                      /*allow_for_broken_orientation_of_coupling_intersections=*/true);
          const auto out_in_coupling_intersection_it_end = out_in_coupling_glue.template iend<0>();
          for (auto out_in_coupling_intersection_it = out_in_coupling_glue.template ibegin<0>();
               out_in_coupling_intersection_it != out_in_coupling_intersection_it_end;
               ++out_in_coupling_intersection_it) {
            ++out_in_num_coupling_intersections;
            const auto& coupling_intersection = *out_in_coupling_intersection_it;
            const auto local_entity = coupling_intersection.inside();
            const size_t local_entity_index = local_grid_view.indexSet().index(local_entity);
            outside_inside_coupling_visualization[macro_neighbor_index][local_entity_index] =
                static_cast<double>(macro_neighbor_index);
            const auto local_neighbor = coupling_intersection.outside();
            const size_t local_neighbor_index = local_neighbor_grid_view.indexSet().index(local_neighbor);
            outside_inside_coupling_visualization[macro_entity_index][local_neighbor_index] =
                static_cast<double>(macro_entity_index);
          }
          if (num_coupling_intersections != out_in_num_coupling_intersections)
            DUNE_THROW(XT::Common::Exceptions::internal_error,
                       "The coupling glue is broken!\n"
                           << "macro entity (local level):   " << macro_entity_index << " (" << local_level << ")\n"
                           << "macro neighbor (local level): " << macro_neighbor_index << " (" << local_neighbor_level
                           << ")");
          logger.debug() << "  " << num_coupling_intersections << " coupling intersections with neighbor "
                         << macro_neighbor_index << std::endl;
        }
      }
    } // walk the macro grid
    vtk_writer.addCellData(subdomain_visualization, "subdomains");
    vtk_writer.addCellData(boundary_visualization, "local boundary entities");
    vtk_writer.addCellData(inside_outside_coupling_visualization, "local coupling entities (inside/outside)");
    vtk_writer.addCellData(outside_inside_coupling_visualization, "local coupling entities (outside/inside)");
    vtk_writer.write(filename, VTK::appendedraw);
  } // ... visualize(...)

  void write_global_visualization(const std::string& filename_prefix,
                                  const XT::Functions::FunctionInterface<dimDomain>& func,
                                  const std::list<int>& subdomains)
  {
    GluedVTKWriter<MacroGridType, LocalGridType, layer> vtk_writer(*this);
    for (size_t subdomain = 0; subdomain < macro_leaf_view_size_; ++subdomain) {
      if (std::find(subdomains.begin(), subdomains.end(), subdomain) != subdomains.end()) {
        auto visualization_adapter = std::make_shared<XT::Functions::VisualizationAdapter<LocalViewType, 1, 1, double>>(
            XT::Functions::make_grid_function<typename LocalViewType::template Codim<0>::Entity>(func), func.name());
        vtk_writer.addVertexData(subdomain, visualization_adapter);
      }
    }
    vtk_writer.write(filename_prefix, VTK::appendedraw);
    vtk_writer.clear();
  } // ... write_global_visualization(...)

private:
  template <class MacroEntityType>
  static std::shared_ptr<LocalGridProviderType> create_grid_of_simplex(const MacroEntityType& macro_entity)
  {
    try {
      GridFactory<LocalGridType> subdomain_factory;
      const auto num_vertices = macro_entity.subEntities(dimDomain);
      std::vector<unsigned int> vertex_ids(num_vertices, 0);
      for (auto&& local_vertex_id : XT::Common::value_range(num_vertices)) {
        const auto vertex = macro_entity.template subEntity<dimDomain>(local_vertex_id).geometry().center();
        subdomain_factory.insertVertex(vertex);
        vertex_ids[local_vertex_id] = local_vertex_id;
      }
      subdomain_factory.insertElement(macro_entity.type(), vertex_ids);
      return std::make_shared<XT::Grid::GridProvider<LocalGridType>>(subdomain_factory.createGrid());
    } catch (GridError& ee) {
      DUNE_THROW(GridError,
                 "It was not possible to create a grid for this simplex with the given GridType!\n\n"
                     << "GridType: " << XT::Common::Typename<LocalGridType>::value() << "\n\n"
                     << "This was the original error: " << ee.what());
    }
  } // ... create_grid_of_simplex(...)

  template <class MacroEntityType>
  static std::shared_ptr<LocalGridProviderType>
  create_grid_of_cube(const MacroEntityType& macro_entity, const std::array<unsigned int, dimDomain>& num_elements)
  {
    const auto num_vertices = macro_entity.subEntities(dimDomain);
    FieldVector<ctype, dimDomain> lower_left(std::numeric_limits<ctype>::max());
    FieldVector<ctype, dimDomain> upper_right(std::numeric_limits<ctype>::min());
    for (auto&& local_vertex_id : XT::Common::value_range(num_vertices)) {
      const auto vertex = macro_entity.template subEntity<dimDomain>(local_vertex_id).geometry().center();
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        lower_left[dd] = std::min(lower_left[dd], vertex[dd]);
        upper_right[dd] = std::max(upper_right[dd], vertex[dd]);
      }
    }
    return std::make_shared<XT::Grid::GridProvider<LocalGridType>>(
        XT::Grid::make_cube_grid<LocalGridType>(lower_left, upper_right, num_elements));
  } // ... create_grid_of_cube(...)

  template <class MacroEntityType>
  static std::shared_ptr<LocalGridProviderType> create_grid_of_cube(const MacroEntityType& macro_entity)
  {
    std::array<unsigned int, dimDomain> num_elements;
    std::fill(num_elements.begin(), num_elements.end(), 1);
    return create_grid_of_cube(macro_entity, num_elements);
  }

  void setup_local_grids()
  {
    const auto& macro_index_set = macro_leaf_view_.indexSet();
    for (auto&& macro_entity : elements(macro_leaf_view_)) {
      auto macro_entity_index = macro_index_set.index(macro_entity);
      if (macro_entity.type().isSimplex())
        local_grids_[macro_entity_index] = create_grid_of_simplex(macro_entity);
      else if (macro_entity.type().isCube())
        local_grids_[macro_entity_index] = create_grid_of_cube(macro_entity);
      else
        DUNE_THROW(GridError, "Unknown entity.type() encountered: " << macro_entity.type());
    }
  } // ... setup_local_grids()

  void setup_local_grids(const std::array<unsigned int, dimDomain>& num_elements_per_subdomain)
  {
    const auto& macro_index_set = macro_leaf_view_.indexSet();
    for (auto&& macro_entity : elements(macro_leaf_view_)) {
      DUNE_THROW_IF(!macro_entity.type().isCube(),
                    GridError,
                    "Prescribing number of elements for local subdomain grids only suitable for cubic subdomains!");
      auto macro_entity_index = macro_index_set.index(macro_entity);
      local_grids_[macro_entity_index] = create_grid_of_cube(macro_entity, num_elements_per_subdomain);
    }
  } // ... setup_local_grids()

  template <Layers l, bool anything = true>
  struct extract_local_view
  {
    template <class G>
    auto operator()(G& g, int /*lv*/) -> decltype(g.leaf_view())
    {
      return g.leaf_view();
    }
  };

  template <bool anything>
  struct extract_local_view<Layers::level, anything>
  {
    template <class G>
    auto operator()(G& g, int lv) -> decltype(g.level_view(lv))
    {
      return g.level_view(lv);
    }
  };

  template <class MacroIntersectionType>
  std::shared_ptr<GlueType> create_glue(const MacroEntityType& macro_entity,
                                        const MacroEntityType& macro_neighbor,
                                        const MacroIntersectionType& macro_intersection,
                                        const int local_entity_level,
                                        const int local_neighbor_level) const
  {
    assert(local_entity_level >= -1);
    assert(local_neighbor_level >= -1);
    assert(local_entity_level <= max_local_level(macro_entity));
    assert(local_neighbor_level <= max_local_level(macro_neighbor));
    const auto& local_entity_grid = local_grid(macro_entity);
    const auto& local_neighbor_grid = local_grid(macro_neighbor);
    auto local_entity_view = extract_local_view<layer>()(local_entity_grid, local_entity_level);
    auto local_neighbor_view = extract_local_view<layer>()(local_neighbor_grid, local_neighbor_level);
    // create descriptors, these can be discarded after creating the extractors
    auto entity_descriptor = create_descriptor(local_entity_view, macro_intersection);
    auto neighbor_descriptor = create_descriptor(local_neighbor_view, macro_intersection);
    // create extractors and merger as shared_ptr, so glue will handle memory
    auto entity_extractor = std::make_shared<LocalExtractorType>(local_entity_view, entity_descriptor);
    auto neighbor_extractor = std::make_shared<LocalExtractorType>(local_neighbor_view, neighbor_descriptor);
    auto contact_merger = std::make_shared<GridGlue::ContactMerge<dimWorld, ctype>>(allowed_overlap_);
    // create glue
    auto glue = std::make_shared<GlueType>(entity_extractor, neighbor_extractor, contact_merger);
    glue->build();
    if (glue->size() == 0)
      DUNE_THROW(GridError,
                 "Something went wrong, the coupling glue is empty!\n"
                     << "   macro_entity " << macro_leaf_view_.indexSet().index(macro_entity) << "\n"
                     << "   local_entity_level " << local_entity_level << "\n"
                     << "   macro_neighbor " << macro_leaf_view_.indexSet().index(macro_neighbor) << "\n"
                     << "   local_neighbor_level " << local_neighbor_level << "\n");
    return glue;
  } // ... create_glue(...)

  void setup_glues(const bool allow_for_broken_orientation_of_coupling_intersections = false)
  {
    const auto& macro_index_set = macro_leaf_view_.indexSet();
    for (auto&& macro_entity : elements(macro_leaf_view_)) {
      const auto macro_entity_index = macro_index_set.index(macro_entity);
      auto& entity_glues = glues_[macro_entity_index];
      // walk the neighbors ...
      for (auto&& macro_intersection : intersections(macro_leaf_view_, macro_entity)) {
        if (macro_intersection.neighbor() && !macro_intersection.boundary()) {
          const auto macro_neighbor_ptr = macro_intersection.outside();
#  if DUNE_VERSION_GTE(DUNE_GRID, 2, 4)
          const auto& macro_neighbor = macro_neighbor_ptr;
#  else
          const auto& macro_neighbor = *macro_neighbor_ptr;
#  endif
          const auto macro_neighbor_index = macro_index_set.index(macro_neighbor);
          for (auto local_entity_level :
               XT::Common::value_range(local_grids_[macro_entity_index]->grid().maxLevel() + 1))
            for (auto local_neighbor_level :
                 XT::Common::value_range(local_grids_[macro_neighbor_index]->grid().maxLevel() + 1)) {
              auto glue = create_glue(
                  macro_entity, macro_neighbor, macro_intersection, local_entity_level, local_neighbor_level);
              if (!allow_for_broken_orientation_of_coupling_intersections) {
                const size_t brocken_intersections = check_for_broken_coupling_intersections(*glue);
                if (brocken_intersections > 0)
                  DUNE_THROW(Exceptions::intersection_orientation_is_broken,
                             "The coupling glue between the grid views of\n"
                                 << "  level " << local_entity_level << " on macro entity   "
                                 << macro_leaf_view_.indexSet().index(macro_entity) << " and\n"
                                 << "  level " << local_neighbor_level << " on macro neighbor "
                                 << macro_leaf_view_.indexSet().index(macro_neighbor) << "\n"
                                 << "contains\n"
                                 << "  " << brocken_intersections << "/" << glue->size() << " intersections with wrong"
                                 << "orientation!");
              }
              entity_glues[macro_neighbor_index][local_entity_level][local_neighbor_level] = glue;
            }
        }
      } // ... walk the neighbors
    }
  } // ... setup_glues(...)

  void prepare_global_grid()
  {
    if (global_grid_)
      return;
    const auto& macro_index_set = macro_leaf_view_.indexSet();
    std::vector<FieldVector<ctype, dimDomain>> vertices;
    std::vector<std::vector<std::vector<unsigned int>>> entity_to_vertex_ids(local_grids_.size());
    std::vector<std::vector<GeometryType>> geometry_types(local_grids_.size());
    // walk the grid for the first time
    for (auto&& macro_entity : elements(macro_leaf_view_)) {
      const auto local_leaf_view = local_grid(macro_entity).leaf_view();
      const auto& local_index_set = local_leaf_view.indexSet();
      const size_t macro_index = macro_index_set.index(macro_entity);
      entity_to_vertex_ids[macro_index] = std::vector<std::vector<unsigned int>>(local_index_set.size(0));
      geometry_types[macro_index] = std::vector<GeometryType>(local_index_set.size(0));
      for (auto&& micro_entity : elements(local_leaf_view)) {
        const size_t micro_index = local_index_set.index(micro_entity);
        const auto num_vertices = micro_entity.
#  if DUNE_VERSION_GTE(DUNE_GRID, 2, 4)

                                  subEntities(dimDomain);
#  else
                                  template count<dimDomain>();
#  endif
        entity_to_vertex_ids[macro_index][micro_index] = std::vector<unsigned int>(num_vertices);
        geometry_types[macro_index][micro_index] = micro_entity.type();
        for (unsigned int local_vertex_id = 0; local_vertex_id < num_vertices; ++local_vertex_id) {
          const unsigned int global_vertex_id = find_insert_vertex(vertices,
                                                                   micro_entity
                                                                       .template subEntity<dimDomain>(local_vertex_id)
#  if DUNE_VERSION_GTE(DUNE_GRID, 2, 4)
                                                                       .
#  else
                                                                       ->
#  endif
                                                                   geometry()
                                                                       .center());
          entity_to_vertex_ids[macro_index][micro_index][local_vertex_id] = global_vertex_id;
        }
      } // walk the local grid
    } // walk the macro grid
    GridFactory<LocalGridType> global_factory;
    for (const auto& vertex : vertices)
      global_factory.insertVertex(vertex);
    size_t II = 0;
    size_t JJ = 0;
    try {
      for (size_t ii = 0; ii < local_grids_.size(); ++ii, II = ii)
        for (size_t jj = 0; jj < entity_to_vertex_ids[ii].size(); ++jj, JJ = jj)
          global_factory.insertElement(geometry_types[ii][jj], entity_to_vertex_ids[ii][jj]);
    } catch (GridError& ee) {
      DUNE_THROW(GridError,
                 "It was not possible to insert an element into the grid factory!\n\n"
                     << "GridType: " << XT::Common::Typename<LocalGridType>::value() << "\n"
                     << "GeometryType: " << geometry_types[II][JJ] << "\n"
                     << "\n"
                     << "This was the original error: " << ee.what());
    } // try
    global_grid_ = std::make_unique<XT::Grid::GridProvider<LocalGridType>>(global_factory.createGrid());

    // build maps of indices, relating the local to the global grid
    const auto global_view = global_grid_->leaf_view();
    local_to_global_indices_ = std::make_unique<std::vector<std::vector<size_t>>>(macro_leaf_view_.indexSet().size(0));
    auto& local_to_global_inds = *local_to_global_indices_;
    global_to_local_indices_ = std::make_unique<std::vector<std::pair<size_t, size_t>>>(global_view.indexSet().size(0));
    auto& global_to_local_inds = *global_to_local_indices_;
    // therefore
    // * create a search on the global view: if all is fine, we only need to walk that one once
    std::vector<FieldVector<ctype, dimDomain>> local_entity_center{FieldVector<ctype, dimDomain>(0.0)};
    auto global_search = XT::Grid::make_entity_in_level_search(global_view);
    // * walk the macro grid
    for (auto&& macro_entity : elements(macro_leaf_view_)) {
      const size_t subd = macro_leaf_view_.indexSet().index(macro_entity);
      const auto local_leaf_view = local_grid(macro_entity).leaf_view();
      const auto& local_index_set = local_leaf_view.indexSet();
      local_to_global_inds[subd] = std::vector<size_t>(local_index_set.size(0));
      // * walk the local grid
      for (auto&& local_entity : elements(local_leaf_view)) {
        const size_t local_entity_index = local_index_set.index(local_entity);
        local_entity_center[0] = local_entity.geometry().center();
        const auto global_entity_ptr_unique_ptrs = global_search(local_entity_center);
        // the search has to be successfull, since the global grid has been constructed to exactly contain each local
        // entity
        assert(global_entity_ptr_unique_ptrs.size() == 1);
        const auto& global_entity_ptr_unique_ptr = global_entity_ptr_unique_ptrs.at(0);
        assert(global_entity_ptr_unique_ptr);
        const auto& global_entity_ptr = *global_entity_ptr_unique_ptr;
        const auto& global_entity = *global_entity_ptr;
        const size_t global_entity_index = global_view.indexSet().index(global_entity);
        // store information
        local_to_global_inds[subd][local_entity_index] = global_entity_index;
        global_to_local_inds[global_entity_index] = {subd, local_entity_index};
      } // * walk the local grid
    } // * walk the macro grid
  } // ... prepare_global_grid(...)

  size_t find_insert_vertex(std::vector<FieldVector<ctype, dimDomain>>& vertices,
                            FieldVector<ctype, dimDomain>&& vertex) const
  {
    // check if vertex is already contained
    for (size_t ii = 0; ii < vertices.size(); ++ii)
      if (XT::Common::FloatCmp::eq(vertex, vertices[ii]))
        return ii;
    // if not, add it
    vertices.emplace_back(std::move(vertex));
    return vertices.size() - 1;
  } // ... find_insert_vertex(...)

  void assert_macro_grid_state() const
  {
    if (macro_leaf_view_.indexSet().size(0) != macro_leaf_view_size_)
      DUNE_THROW(InvalidStateException,
                 "The size of the macro leaf grid view has changed (from "
                     << macro_leaf_view_size_ << " to " << macro_leaf_view_.indexSet().size(0)
                     << "), which invalidates this object!\n"
                     << "Do not adapt the macro grid after instantiating this object!");
  } // void assert_macro_grid_state()

  MacroGridProviderType& macro_grid_;
  const ctype allowed_overlap_;
  MacroGridViewType macro_leaf_view_;
  const size_t macro_leaf_view_size_;
  std::vector<std::shared_ptr<LocalGridProviderType>> local_grids_;
  std::vector<std::map<size_t, std::map<int, std::map<int, std::shared_ptr<GlueType>>>>> glues_;
  std::map<size_t, std::map<int, std::vector<std::pair<MicroEntityType, std::vector<int>>>>>
      macro_entity_to_local_level_to_boundary_entity_ptrs_with_local_intersections_;
  std::unique_ptr<LocalGridProviderType> global_grid_;
  std::unique_ptr<std::vector<std::vector<size_t>>> local_to_global_indices_;
  std::unique_ptr<std::vector<std::pair<size_t, size_t>>> global_to_local_indices_;
}; // class Glued


template <class MacroGridType, class LocalGridType, Layers layer>
class GluedVTKWriter
{
  using LocalGridViewType = typename Layer<LocalGridType, layer, Backends::view>::type;

  // we only need this class to access a protected pwrite method of VTKWriter
  class LocalVTKWriter : public VTKWriter<LocalGridViewType>
  {
    using BaseType = VTKWriter<LocalGridViewType>;

  public:
    LocalVTKWriter(const LocalGridViewType& local_grid_view, const size_t subdomain, const size_t num_subdomains)
      : BaseType(local_grid_view)
      , commRank_(boost::numeric_cast<int>(subdomain))
      , commSize_(boost::numeric_cast<int>(num_subdomains))
    {}

    void write_locally(const std::string& name, VTK::OutputType ot)
    {
      BaseType::pwrite(name, "", "", ot, commRank_, commSize_);
    }

  private:
    const int commRank_;
    const int commSize_;
  }; // class LocalVTKWriter

public:
  using GluedGridType = Glued<MacroGridType, LocalGridType, layer>;

  GluedVTKWriter(const GluedGridType& glued_grid, const int local_level = -1)
    : glued_grid_(glued_grid)
    , local_levels_(glued_grid_.num_subdomains(), -1)
  {
    // set each local level to its respective max
    if (local_level < 0)
      for (auto&& macro_entity : elements(glued_grid_.macro_grid_view()))
        local_levels_[glued_grid_.subdomain(macro_entity)] = glued_grid_.max_local_level(macro_entity);
    prepare_local_vtk_writers();
  } // GluedVTKWriter(...)

  GluedVTKWriter(const GluedGridType& glued_grid, std::vector<int> local_levels)
    : glued_grid_(glued_grid)
    , local_levels_(std::move(local_levels))
  {
    for (size_t ss = 0; ss < glued_grid_.num_subdomains(); ++ss)
      if (local_levels_[ss] < 0)
        local_levels_[ss] = glued_grid_.max_local_level(ss);
    prepare_local_vtk_writers();
  }

  template <class V>
  void addCellData(const size_t subdomain, const std::vector<V>& vector, const std::string& name, const int ncomps = 1)
  {
    DUNE_THROW_IF(subdomain > glued_grid_.num_subdomains(),
                  Common::Exceptions::index_out_of_range,
                  "subdomain = " << subdomain
                                 << "\n   glued_grid_.num_subdomains() = " << glued_grid_.num_subdomains());
    local_vtk_writers_[subdomain]->addCellData(vector, name, ncomps);
  } // ... addCellData(...)

  template <class V>
  void addCellData(const std::vector<std::vector<V>>& vectors, const std::string& name, const int ncomps = 1)
  {
    DUNE_THROW_IF(vectors.size() != glued_grid_.num_subdomains(),
                  Common::Exceptions::shapes_do_not_match,
                  "vectors.size() =  " << vectors.size()
                                       << "\n   glued_grid_.num_subdomains(): " << glued_grid_.num_subdomains());
    for (size_t ss = 0; ss < glued_grid_.num_subdomains(); ++ss)
      local_vtk_writers_[ss]->addCellData(vectors[ss], name, ncomps);
  } // ... addCellData(...)

  template <class VTKFunctionType>
  void addCellData(const size_t subdomain, const std::shared_ptr<VTKFunctionType>& function, const std::string& name)
  {
    DUNE_THROW_IF(subdomain > glued_grid_.num_subdomains(),
                  Common::Exceptions::index_out_of_range,
                  "subdomain = " << subdomain
                                 << "\n   glued_grid_.num_subdomains() = " << glued_grid_.num_subdomains());
    local_vtk_writers_[subdomain]->addCellData(function, name);
  } // ... addCellData(...)

  template <class VTKFunctionType>
  void addCellData(const std::vector<std::shared_ptr<VTKFunctionType>>& functions, const std::string& name)
  {
    DUNE_THROW_IF(functions.size() != glued_grid_.num_subdomains(),
                  Common::Exceptions::shapes_do_not_match,
                  "functions.size() =  " << functions.size()
                                         << "\n   glued_grid_.num_subdomains(): " << glued_grid_.num_subdomains());
    for (size_t ss = 0; ss < glued_grid_.num_subdomains(); ++ss)
      local_vtk_writers_[ss]->addCellData(functions[ss], name);
  } // ... addCellData(...)

  template <class V>
  void
  addVertexData(const size_t subdomain, const std::vector<V>& vector, const std::string& name, const int ncomps = 1)
  {
    DUNE_THROW_IF(subdomain > glued_grid_.num_subdomains(),
                  Common::Exceptions::index_out_of_range,
                  "subdomain = " << subdomain
                                 << "\n   glued_grid_.num_subdomains() = " << glued_grid_.num_subdomains());
    local_vtk_writers_[subdomain]->addVertexData(vector, name, ncomps);
  } // ... addVertexData(...)

  template <class V>
  void addVertexData(const std::vector<std::vector<V>>& vectors, const std::string& name, const int ncomps = 1)
  {
    DUNE_THROW_IF(vectors.size() != glued_grid_.num_subdomains(),
                  Common::Exceptions::shapes_do_not_match,
                  "vectors.size() =  " << vectors.size()
                                       << "\n   glued_grid_.num_subdomains(): " << glued_grid_.num_subdomains());
    for (size_t ss = 0; ss < glued_grid_.num_subdomains(); ++ss)
      local_vtk_writers_[ss]->addVertexData(vectors[ss], name, ncomps);
  } // ... addVertexData(...)

  template <class VTKFunctionType>
  void addVertexData(const size_t subdomain, const std::shared_ptr<VTKFunctionType>& function)
  {
    DUNE_THROW_IF(subdomain > glued_grid_.num_subdomains(),
                  Common::Exceptions::index_out_of_range,
                  "subdomain = " << subdomain
                                 << "\n   glued_grid_.num_subdomains() = " << glued_grid_.num_subdomains());
    local_vtk_writers_[subdomain]->addVertexData(function);
  } // ... addVertexData(...)

  template <class VTKFunctionType>
  void addVertexData(const std::vector<std::shared_ptr<VTKFunctionType>>& functions)
  {
    DUNE_THROW_IF(functions.size() != glued_grid_.num_subdomains(),
                  Common::Exceptions::shapes_do_not_match,
                  "functions.size() =  " << functions.size()
                                         << "\n   glued_grid_.num_subdomains(): " << glued_grid_.num_subdomains());
    for (size_t ss = 0; ss < glued_grid_.num_subdomains(); ++ss)
      local_vtk_writers_[ss]->addVertexData(functions[ss]);
  } // ... addVertexData(...)

  void clear()
  {
    for (auto& local_vtk_writer : local_vtk_writers_)
      local_vtk_writer->clear();
  }

  void write(const std::string& name, VTK::OutputType type = VTK::ascii)
  {
    for (size_t ss = 0; ss < glued_grid_.num_subdomains(); ++ss)
      local_vtk_writers_[ss]->write_locally(name, type);
  }

private:
  template <Layers l, bool anything = true>
  struct extract_local_view
  {
    template <class G>
    auto operator()(G& g, int /*lv*/) -> decltype(g.leaf_view())
    {
      return g.leaf_view();
    }
  };

  template <bool anything>
  struct extract_local_view<Layers::level, anything>
  {
    template <class G>
    auto operator()(G& g, int lv) -> decltype(g.level_view(lv))
    {
      return g.level_view(lv);
    }
  };

  void prepare_local_vtk_writers()
  {
    for (auto&& macro_entity : elements(glued_grid_.macro_grid_view())) {
      const size_t subdomain = glued_grid_.subdomain(macro_entity);
      local_vtk_writers_.emplace_back(new LocalVTKWriter(
          extract_local_view<layer>()(glued_grid_.local_grid(macro_entity), local_levels_[subdomain]),
          subdomain,
          glued_grid_.num_subdomains()));
    }
  } // ... prepare_local_vtk_writers(...)

  const GluedGridType& glued_grid_;
  std::vector<int> local_levels_;
  std::vector<std::unique_ptr<LocalVTKWriter>> local_vtk_writers_;
}; // class GluedVTKWriter


#else // HAVE_DUNE_GRID_GLUE


template <class MacroGridImp, class LocalGridImp, Layers layer = Layers::level>
class GluedVTKWriter
{
  static_assert(AlwaysFalse<MacroGridImp>::value, "You are missing dune-grid-glue!");
};


template <class MacroGridImp, class LocalGridImp, Layers layer = Layers::level>
class Glued
{
  static_assert(AlwaysFalse<MacroGridImp>::value, "You are missing dune-grid-glue!");
};


#endif // HAVE_DUNE_GRID_GLUE


} // namespace Dune::XT::Grid::DD

#endif // DUNE_XT_GRID_DD_GLUED_HH
