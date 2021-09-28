// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   René Fritze     (2013 - 2016, 2018 - 2020)
//   Tobias Leibner  (2014 - 2017, 2020)

#ifndef DUNE_XT_GRID_SEARCH_HH
#define DUNE_XT_GRID_SEARCH_HH

#include <vector>

#include <boost/range/iterator_range.hpp>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/parallel/threadstorage.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/type_traits.hh>

namespace Dune::XT::Grid {


template <class GridLayerType, int codim = 0>
class EntitySearchBase
{
  static_assert(is_layer<GridLayerType>::value);

public:
  using EntityType = extract_entity_t<GridLayerType, codim>;
  using LocalCoordinateType = typename EntityType::Geometry::LocalCoordinate;
  using GlobalCoordinateType = typename EntityType::Geometry::GlobalCoordinate;
  using EntityVectorType = std::vector<std::unique_ptr<EntityType>>;
}; // class EntitySearchBase


template <int codim>
struct CheckInside
{
  template <class GeometryType, class GlobalCoordinateType>
  static bool check(const GeometryType& geometry, const GlobalCoordinateType& point)
  {
    if (codim == point.size()) {
      return Common::FloatCmp::eq(geometry.center(), point);
    } else {
      if (Common::FloatCmp::ne(
              geometry.global(geometry.local(point)),
              point,
              2. * Common::FloatCmp::DefaultEpsilon<typename GlobalCoordinateType::value_type>::value(),
              2. * Common::FloatCmp::DefaultEpsilon<typename GlobalCoordinateType::value_type>::value()))
        return false;
      const auto& refElement = reference_element(geometry);
      return refElement.checkInside(geometry.local(point));
    }
  }
};


template <>
struct CheckInside<0>
{
  template <class GeometryType, class GlobalCoordinateType>
  static bool check(const GeometryType& geometry, const GlobalCoordinateType& point)
  {
    const auto& refElement = reference_element(geometry);
    return refElement.checkInside(geometry.local(point));
  }
};


/** Provides a facility to search a given grid layer for entities with arbitrary codim
 * that contain a set of points. The search is "in-level", meaning no grid hierarchy
 * is used in the search. The search position iterator on the grid persists
 * between searches, reducing complexity of repeated searches on the grid.
 * \attention This makes it inherently not thread safe
 **/
template <class GridLayerType, int codim = 0>
class EntityInlevelSearch : public EntitySearchBase<GridLayerType, codim>
{
  using BaseType = EntitySearchBase<GridLayerType, codim>;
  using IteratorType = extract_iterator_t<GridLayerType, codim>;

public:
  using EntityVectorType = typename BaseType::EntityVectorType;

private:
  inline typename EntityVectorType::value_type check_add(const typename BaseType::EntityType& entity,
                                                         const typename BaseType::GlobalCoordinateType& point) const
  {
    const auto& geometry = entity.geometry();
    if (CheckInside<codim>::check(geometry, point)) {
      return std::make_unique<typename BaseType::EntityType>(entity);
    }
    return nullptr;
  }

public:
  EntityInlevelSearch(const GridLayerType& grid_layer)
    : grid_layer_(grid_layer)
    , it_last_(grid_layer_.template begin<codim>())
  {}

  /** \arg points iterable sequence of global coordinates to search for
   *  \return a vector of size points.size() of, potentially nullptr if no corresponding one was found,
   *          unique_ptr<Entity>
   **/
  template <class PointContainerType>
  EntityVectorType operator()(const PointContainerType& points)
  {
    const IteratorType begin = grid_layer_.template begin<codim>();
    const IteratorType end = grid_layer_.template end<codim>();
    EntityVectorType ret(points.size());
    typename EntityVectorType::size_type idx(0);
    for (const auto& point : points) {
      IteratorType it_current = *it_last_;
      bool it_reset = true;
      typename EntityVectorType::value_type tmp_ptr(nullptr);
      for (; it_current != end; ++it_current) {
        if ((tmp_ptr = check_add(*it_current, point))) {
          ret[idx++] = std::move(tmp_ptr);
          tmp_ptr = nullptr;
          it_reset = false;
          *it_last_ = it_current;
          break;
        }
      }
      if (!it_reset)
        continue;
      for (it_current = begin; it_current != *it_last_; ++it_current) {
        if ((tmp_ptr = check_add(*it_current, point))) {
          ret[idx++] = std::move(tmp_ptr);
          tmp_ptr = nullptr;
          it_reset = false;
          *it_last_ = it_current;
          break;
        }
      }
      if (!it_reset)
        continue;
      idx++;
    }
    return ret;
  } // ... operator()

private:
  const GridLayerType grid_layer_;
  Common::PerThreadValue<IteratorType> it_last_;
}; // class EntityInlevelSearch


/** Like EntityInlevelSearch, but works also for grids where there is no iterator for the given
 * codimension (e.g. UGGrid for codim != dim && codim != 0)
 * \see EntityInlevelSearch
 **/
template <class GridLayerType, int codim = 0>
class FallbackEntityInlevelSearch : public EntitySearchBase<GridLayerType>
{
  using BaseType = EntitySearchBase<GridLayerType, codim>;
  using IteratorType = typename extract_iterator<GridLayerType, 0>::type;

public:
  using EntityVectorType = typename BaseType::EntityVectorType;

private:
  inline typename EntityVectorType::value_type check_add(const typename BaseType::EntityType& entity,
                                                         const typename BaseType::GlobalCoordinateType& point) const
  {
    const auto& geometry = entity.geometry();
    if (CheckInside<codim>::check(geometry, point)) {
      return std::make_unique<typename BaseType::EntityType>(entity);
    }
    return nullptr;
  }

public:
  FallbackEntityInlevelSearch(const GridLayerType& grid_layer)
    : grid_layer_(grid_layer)
    , it_last_(grid_layer_.template begin<0>())
  {}

  /** \arg points iterable sequence of global coordinates to search for
   *  \return a vector of size points.size() of, potentially nullptr if no corresponding one was found,
   *          unique_ptr<Entity>
   **/
  template <class PointContainerType>
  EntityVectorType operator()(const PointContainerType& points)
  {
    const IteratorType begin = grid_layer_.template begin<0>();
    const IteratorType end = grid_layer_.template end<0>();
    EntityVectorType ret(points.size());
    typename EntityVectorType::size_type idx(0);
    for (const auto& point : points) {
      IteratorType it_current = it_last_;
      bool it_reset = true;
      typename EntityVectorType::value_type tmp_ptr(nullptr);
      for (; it_current != end; ++it_current) {
        const auto& entity = *it_current;
        for (unsigned int local_index = 0; local_index < entity.subEntities(codim); ++local_index) {
          const auto& sub_entity = entity.template subEntity<codim>(local_index);
          if ((tmp_ptr = check_add(sub_entity, point))) {
            ret[idx++] = std::move(tmp_ptr);
            tmp_ptr = nullptr;
            it_reset = false;
            it_last_ = it_current;
            break;
          }
        } // loop over subentities
        if (!it_reset)
          break;
      } // loop over codim 0 entities
      if (!it_reset)
        continue;
      for (it_current = begin; it_current != it_last_; ++it_current) {
        const auto& entity = *it_current;
        for (unsigned int local_index = 0; local_index < entity.subEntities(codim); ++local_index) {
          const auto& sub_entity = entity.template subEntity<codim>(local_index);
          if ((tmp_ptr = check_add(sub_entity, point))) {
            ret[idx++] = std::move(tmp_ptr);
            tmp_ptr = nullptr;
            it_reset = false;
            it_last_ = it_current;
            break;
          }
        } // loop over subentities
        if (!it_reset)
          break;
      } // loop over codim 0 entities
      if (!it_reset)
        continue;
      idx++;
    } // loop over points
    return ret;
  } // ... operator()

private:
  const GridLayerType grid_layer_;
  IteratorType it_last_;
}; // class FallbackEntityInlevelSearch


template <class GridLayerType>
class EntityHierarchicSearch : public EntitySearchBase<GridLayerType>
{
  using BaseType = EntitySearchBase<GridLayerType>;

  const GridLayerType grid_layer_;
  const int start_level_;

public:
  EntityHierarchicSearch(const GridLayerType& grid_layer)
    : grid_layer_(grid_layer)
    , start_level_(0)
  {}

  using EntityVectorType = typename BaseType::EntityVectorType;

  template <class PointContainerType>
  EntityVectorType operator()(const PointContainerType& points) const
  {
    auto level = std::min(grid_layer_.grid().maxLevel(), start_level_);
    auto range = elements(grid_layer_.grid().levelView(level));
    return process(points, range);
  }

private:
  template <class QuadpointContainerType, class RangeType>
  EntityVectorType process(const QuadpointContainerType& quad_points, const RangeType& range) const
  {
    EntityVectorType ret;

    for (const auto& my_ent : range) {
      const auto my_level = my_ent.level();
      const auto& geometry = my_ent.geometry();
      const auto& refElement = reference_element(geometry);
      for (const auto& point : quad_points) {
        if (refElement.checkInside(geometry.local(point))) {
          // if I cannot descend further add this entity even if it's not my view
          if (grid_layer_.grid().maxLevel() <= my_level || grid_layer_.contains(my_ent)) {
            ret.emplace_back(my_ent);
          } else {
            const auto h_end = my_ent.hend(my_level + 1);
            const auto h_begin = my_ent.hbegin(my_level + 1);
            const auto h_range = boost::make_iterator_range(h_begin, h_end);
            const auto kids = process(QuadpointContainerType(1, point), h_range);
            ret.insert(ret.end(), kids.begin(), kids.end());
          }
        }
      }
    }
    return ret;
  }
}; // class EntityHierarchicSearch


template <class GV, int codim = 0>
EntityInlevelSearch<GV, codim> make_entity_in_level_search(const GV& grid_view)
{
  return EntityInlevelSearch<GV, codim>(grid_view);
}


template <class GV>
EntityHierarchicSearch<GV> make_entity_hierarchic_search(const GV& grid_view)
{
  return EntityHierarchicSearch<GV>(grid_view);
}


} // namespace Dune::XT::Grid

#endif // DUNE_XT_GRID_SEARCH_HH
