// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2013 - 2015)
//   Rene Milk       (2013 - 2015)
//   Tobias Leibner  (2014 - 2015)

#ifndef DUNE_XT_GRID_SEARCH_HH
#define DUNE_XT_GRID_SEARCH_HH

#if HAVE_DUNE_GRID

#include <vector>

#include <boost/range/iterator_range.hpp>

#if HAVE_DUNE_GEOMETRY
#include <dune/geometry/referenceelements.hh>
#else
#include <dune/grid/common/genericreferenceelements.hh>
#endif
#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/ranges.hh>

#include <dune/xt/grid/entity.hh>

namespace Dune {
namespace XT {
namespace Grid {

template <class GridViewType>
class EntitySearchBase
{
  typedef typename GridViewType::Traits ViewTraits;
  static_assert(std::is_base_of<GridView<ViewTraits>, GridViewType>::value,
                "GridViewType has to be derived from GridView!");

public:
  typedef typename ViewTraits::template Codim<0>::Entity EntityType;
  typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;
  typedef typename EntityType::Geometry::GlobalCoordinate GlobalCoordinateType;
  typedef std::vector<std::unique_ptr<EntityType>> EntityVectorType;
}; // class EntitySearchBase

template <class GridViewType>
class EntityInlevelSearch : public EntitySearchBase<GridViewType>
{
  typedef EntitySearchBase<GridViewType> BaseType;

  typedef typename GridViewType::template Codim<0>::Iterator IteratorType;

public:
  typedef typename BaseType::EntityVectorType EntityVectorType;

private:
  inline typename EntityVectorType::value_type check_add(const typename BaseType::EntityType& entity,
                                                         const typename BaseType::GlobalCoordinateType& point) const
  {
    const auto& geometry   = entity.geometry();
    const auto& refElement = reference_element(geometry);
    if (refElement.checkInside(geometry.local(point))) {
      return Common::make_unique<typename BaseType::EntityType>(entity);
    }
    return nullptr;
  }

public:
  EntityInlevelSearch(const GridViewType& gridview)
    : gridview_(gridview)
    , it_last_(gridview_.template begin<0>())
  {
  }

  template <class PointContainerType>
  EntityVectorType operator()(const PointContainerType& points)
  {
    const IteratorType begin = gridview_.template begin<0>();
    const IteratorType end = gridview_.template end<0>();
    EntityVectorType ret(points.size());
    typename EntityVectorType::size_type idx(0);
    for (const auto& point : points) {
      IteratorType it_current = it_last_;
      bool it_reset = true;
      typename EntityVectorType::value_type tmp_ptr(nullptr);
      for (; it_current != end; ++it_current) {
        if ((tmp_ptr = check_add(*it_current, point))) {
          ret[idx++] = std::move(tmp_ptr);
          tmp_ptr    = nullptr;
          it_reset   = false;
          it_last_   = it_current;
          break;
        }
      }
      if (!it_reset)
        continue;
      for (it_current = begin; it_current != it_last_; ++it_current) {
        if ((tmp_ptr = check_add(*it_current, point))) {
          ret[idx++] = std::move(tmp_ptr);
          tmp_ptr    = nullptr;
          it_reset   = false;
          it_last_   = it_current;
          break;
        }
      }
    }
    return ret;
  } // ... operator()

private:
  const GridViewType gridview_;
  IteratorType it_last_;
}; // class EntityInlevelSearch

template <class GridViewType>
class EntityHierarchicSearch : public EntitySearchBase<GridViewType>
{
  typedef EntitySearchBase<GridViewType> BaseType;

  const GridViewType gridview_;
  const int start_level_;

public:
  EntityHierarchicSearch(const GridViewType& gridview)
    : gridview_(gridview)
    , start_level_(0)
  {
  }

  typedef typename BaseType::EntityVectorType EntityVectorType;

  template <class PointContainerType>
  EntityVectorType operator()(const PointContainerType& points) const
  {
    auto level = std::min(gridview_.grid().maxLevel(), start_level_);
    auto range = elements(gridview_.grid().levelView(level));
    return process(points, range);
  }

private:
  template <class QuadpointContainerType, class RangeType>
  EntityVectorType process(const QuadpointContainerType& quad_points, const RangeType& range) const
  {
    EntityVectorType ret;

    for (const auto& my_ent : range) {
      const auto my_level    = my_ent.level();
      const auto& geometry   = my_ent.geometry();
      const auto& refElement = reference_element(geometry);
      for (const auto& point : quad_points) {
        if (refElement.checkInside(geometry.local(point))) {
          // if I cannot descend further add this entity even if it's not my view
          if (gridview_.grid().maxLevel() <= my_level || gridview_.contains(my_ent)) {
            ret.emplace_back(my_ent);
          } else {
            const auto h_end   = my_ent.hend(my_level + 1);
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

template <class GV>
EntityInlevelSearch<GV> make_entity_in_level_search(const GV& grid_view)
{
  return EntityInlevelSearch<GV>(grid_view);
}

template <class GV>
EntityHierarchicSearch<GV> make_entity_hierarchic_search(const GV& grid_view)
{
  return EntityHierarchicSearch<GV>(grid_view);
}

} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // HAVE_DUNE_GRID

#endif // DUNE_XT_GRID_SEARCH_HH
