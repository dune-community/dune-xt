// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2018, 2020)
//   René Fritze     (2014 - 2020)
//   Tobias Leibner  (2015 - 2020)

#ifndef DUNE_XT_GRID_FILTERS_INTERSECTION_HH
#define DUNE_XT_GRID_FILTERS_INTERSECTION_HH

#include <functional>

#include <dune/grid/common/partitionset.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/timedlogging.hh>

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/print.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/print.hh>

#include "base.hh"

namespace Dune::XT::Grid::ApplyOn {


/**
 * \brief A filter which selects all intersections.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GV>
class AllIntersections : public IntersectionFilter<GV>
{
  using BaseType = IntersectionFilter<GV>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit AllIntersections() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new AllIntersections<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& /*intersection*/) const override final
  {
    return true;
  }
}; // class AllIntersections


/**
 * \brief A filter which selects all intersections only once.
 *
 * \note Use a proper mapper, the indixes are not unique if the view contains more than one codim 0 element type!
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GV>
class AllIntersectionsOnce : public IntersectionFilter<GV>
{
  using BaseType = IntersectionFilter<GV>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit AllIntersectionsOnce() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new AllIntersectionsOnce<GridViewType>();
  }

  bool contains(const GridViewType& grid_layer, const IntersectionType& intersection) const override final
  {
    if (!intersection.neighbor())
      return true;
    else {
      const auto inside_element = intersection.inside();
      const auto outside_element = intersection.outside();
      return grid_layer.indexSet().index(inside_element) < grid_layer.indexSet().index(outside_element);
    }
  }
}; // class AllIntersectionsOnce


/**
 * \brief A filter which selects no intersections.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GV>
class NoIntersections : public IntersectionFilter<GV>
{
  using BaseType = IntersectionFilter<GV>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit NoIntersections() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new NoIntersections<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& /*intersection*/) const override final
  {
    return false;
  }
}; // class NoIntersections


/**
 * \brief A filter which selects each inner intersection.
 *
 * \note To decide if this in an inner intersection,
\code
intersection.neighbor() && !intersection.boundary()
\endcode
 *       is used.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GV>
class InnerIntersections : public IntersectionFilter<GV>
{
  using BaseType = IntersectionFilter<GV>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit InnerIntersections() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new InnerIntersections<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return intersection.neighbor() && !intersection.boundary();
  }
}; // class InnerIntersections


/**
 * \brief A filter which selects each inner intersection only once.
 *
 * \note Use a proper mapper, the indixes are not unique if the view contains more than one codim 0 element type!
 *
 * \note To decide if this in an inner intersection,
\code
intersection.neighbor() && !intersection.boundary()
\endcode
 *       is used, and true is returned, if the index of the inside() element is smaller than the index of the outside()
 *       element.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GV>
class InnerIntersectionsOnce : public IntersectionFilter<GV>
{
  using BaseType = IntersectionFilter<GV>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit InnerIntersectionsOnce() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new InnerIntersectionsOnce<GridViewType>();
  }

  bool contains(const GridViewType& grid_layer, const IntersectionType& intersection) const override final
  {
    if (intersection.neighbor() && !intersection.boundary()) {
      const auto inside_element = intersection.inside();
      const auto outside_element = intersection.outside();
      const auto inside_index = grid_layer.indexSet().index(inside_element);
      const auto outside_index = grid_layer.indexSet().index(outside_element);
      return inside_index < outside_index;
    } else
      return false;
  }
}; // class InnerIntersectionsOnce

/**
 * \brief A filter which selects each inner intersection only once. Other than InnerIntersectionsOnce,
 * this also works if we only walk over a subset of the grid elements (if the correct map is provided).
 *
 * InnerIntersectionOnce selects each inner intersection only once. However, this only works if we always
 * visit the elements on both sides of each intersection. If we only walk over a subset of the grid elements,
 * we will visit some of the intersections only from one side. In that case, InnerIntersectionsOnceMap can
 * be used. However, as it is hard to tell which intersections are visited from both sides, the user has to
 * provide a map which contains for each entity index a set of outer element indices for which the
 * intersection should be skipped. To visit each intersection exactly once, one side of intersections which
 * are visited from both sides has to show up in the ignore map.
 * The map can be created from the set of entities by calling the fill_ignore_map method.
 *
 * \sa InnerIntersectionsOnce
 */
template <class GV>
class InnerIntersectionsOnceMap : public IntersectionFilter<GV>
{
  using BaseType = IntersectionFilter<GV>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  static void fill_ignore_map(std::map<size_t, std::set<size_t>>& outside_indices_to_ignore,
                              const std::set<extract_entity_t<GV>, EntityLess<GV>>& entity_range,
                              const GV& grid_view)
  {

    for (auto&& entity : entity_range) {
      const auto inside_index = grid_view.indexSet().index(entity);
      outside_indices_to_ignore[inside_index] = std::set<size_t>();
      for (auto&& intersection : intersections(grid_view, entity)) {
        if (intersection.neighbor() && entity_range.count(intersection.outside())) {
          const auto outside_index = grid_view.indexSet().index(intersection.outside());
          if (inside_index > outside_index)
            outside_indices_to_ignore[inside_index].insert(outside_index);
        }
      }
    }
  }

  InnerIntersectionsOnceMap(const std::map<size_t, std::set<size_t>>& outside_indices_to_ignore)
    : outside_indices_to_ignore_(outside_indices_to_ignore)
  {}

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new InnerIntersectionsOnceMap<GridViewType>(outside_indices_to_ignore_);
  }

  bool contains(const GridViewType& grid_layer, const IntersectionType& intersection) const override final
  {
    if (intersection.neighbor() && !intersection.boundary()) {
      const auto inside_element = intersection.inside();
      const auto outside_element = intersection.outside();
      const auto inside_index = grid_layer.indexSet().index(inside_element);
      const auto outside_index = grid_layer.indexSet().index(outside_element);
      if (outside_indices_to_ignore_.at(inside_index).count(outside_index))
        return false;
      else
        return true;
    } else
      return false;
  }

  const std::map<size_t, std::set<size_t>>& outside_indices_to_ignore_;
}; // class InnerIntersectionsOnceMap


/**
 * \brief Selects each inner intersection in given partition only once.
 *
 * \note Use a proper mapper, the indixes are not unique if the view contains more than one codim 0 element type!
 *
 * \sa InnerIntersectionsOnce
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GV, class PartitionSetType>
class PartitionSetInnerIntersectionsOnce : public IntersectionFilter<GV>
{
  using BaseType = IntersectionFilter<GV>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit PartitionSetInnerIntersectionsOnce() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new PartitionSetInnerIntersectionsOnce<GridViewType, PartitionSetType>();
  }

  bool contains(const GridViewType& grid_layer, const IntersectionType& intersection) const override final
  {
    if (intersection.neighbor() && !intersection.boundary()
        && PartitionSetType::contains(intersection.inside().partitionType())) {
      const auto inside_element = intersection.inside();
      const auto outside_element = intersection.outside();
      if (!PartitionSetType::contains(intersection.outside().partitionType()))
        return true;
      return grid_layer.indexSet().index(inside_element) < grid_layer.indexSet().index(outside_element);
    } else
      return false;
  }
}; // class PartitionSetInnerIntersectionsOnce


/**
 * \brief A filter which selects boundary intersections.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GV>
class BoundaryIntersections : public IntersectionFilter<GV>
{
  using BaseType = IntersectionFilter<GV>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit BoundaryIntersections() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new BoundaryIntersections<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return intersection.boundary();
  }
}; // class BoundaryIntersections


/**
 * \brief A filter which selects intersections on a non-periodic boundary.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GV>
class NonPeriodicBoundaryIntersections : public IntersectionFilter<GV>
{
  using BaseType = IntersectionFilter<GV>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit NonPeriodicBoundaryIntersections() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new NonPeriodicBoundaryIntersections<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return intersection.boundary() && !intersection.neighbor();
  }
}; // class BoundaryIntersections


/**
 * \brief A filter which selects intersections on a periodic boundary.
 *
 *        To decide if this a periodic intersection,
\code
intersection.neighbor() && intersection.boundary()
\endcode
 *        is used.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GV>
class PeriodicBoundaryIntersections : public IntersectionFilter<GV>
{
  using BaseType = IntersectionFilter<GV>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit PeriodicBoundaryIntersections() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new PeriodicBoundaryIntersections<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return intersection.neighbor() && intersection.boundary();
  }
}; // class PeriodicBoundaryIntersections


/**
 * \brief A filter which selects intersections on a periodic boundary only once.
 *
 * \note Use a proper mapper, the indixes are not unique if the view contains more than one codim 0 element type!
 *
 *        To decide if this in an periodic intersection,
\code
intersection.neighbor() && intersection.boundary()
\endcode
 *        is used, and true is returned, if the index of the inside() element is smaller than the index of the outside()
 *        element.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GV>
class PeriodicBoundaryIntersectionsOnce : public IntersectionFilter<GV>
{
  using BaseType = IntersectionFilter<GV>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit PeriodicBoundaryIntersectionsOnce() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new PeriodicBoundaryIntersectionsOnce<GridViewType>();
  }

  bool contains(const GridViewType& grid_layer, const IntersectionType& intersection) const override final
  {
    if (intersection.neighbor() && intersection.boundary()) {
      const auto inside_element = intersection.inside();
      const auto outside_element = intersection.outside();
      return grid_layer.indexSet().index(inside_element) < grid_layer.indexSet().index(outside_element);
    } else {
      return false;
    }
  }
}; // class PeriodicBoundaryIntersectionsOnce

/**
 * \brief A filter which selects each periodic boundary intersections only once. Other than
 * PeriodicBoundaryIntersectionsOnce, this also works if we only walk over a subset of the grid elements
 * (if the correct map is provided).
 *
 * \sa InnerIntersectionsOnceMap
 */
template <class GV>
class PeriodicBoundaryIntersectionsOnceMap : public IntersectionFilter<GV>
{
  using BaseType = IntersectionFilter<GV>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  static void fill_ignore_map(std::map<size_t, std::set<size_t>>& outside_indices_to_ignore,
                              const std::set<extract_entity_t<GV>, EntityLess<GV>>& entity_range,
                              const GV& grid_view)
  {
    InnerIntersectionsOnceMap<GV>::fill_ignore_map(outside_indices_to_ignore, entity_range, grid_view);
  }

  PeriodicBoundaryIntersectionsOnceMap(const std::map<size_t, std::set<size_t>>& outside_indices_to_ignore)
    : outside_indices_to_ignore_(outside_indices_to_ignore)
  {}

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new PeriodicBoundaryIntersectionsOnceMap<GridViewType>(outside_indices_to_ignore_);
  }

  bool contains(const GridViewType& grid_layer, const IntersectionType& intersection) const override final
  {
    if (intersection.neighbor() && intersection.boundary()) {
      const auto inside_element = intersection.inside();
      const auto outside_element = intersection.outside();
      const auto inside_index = grid_layer.indexSet().index(inside_element);
      const auto outside_index = grid_layer.indexSet().index(outside_element);
      if (outside_indices_to_ignore_.at(inside_index).count(outside_index))
        return false;
      else
        return true;
    } else {
      return false;
    }
  }

  const std::map<size_t, std::set<size_t>>& outside_indices_to_ignore_;
}; // class PeriodicBoundaryIntersectionsOnceMap

/**
 * \brief A filter which selects intersections based on a lambda expression.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GV>
class GenericFilteredIntersections : public IntersectionFilter<GV>
{
  using BaseType = IntersectionFilter<GV>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;
  using GenericFunctionType = std::function<bool(const GridViewType&, const IntersectionType&)>;

  explicit GenericFilteredIntersections(GenericFunctionType func)
    : filter_(func)
  {}

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new GenericFilteredIntersections<GridViewType>(filter_);
  }

  bool contains(const GridViewType& grid_layer, const IntersectionType& intersection) const override final
  {
    return filter_(grid_layer, intersection);
  }

private:
  const GenericFunctionType filter_;
}; // class BoundaryIntersections


/**
 * \brief A filter which selects intersections on a given part of the boundary.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GV>
class CustomBoundaryIntersections : public IntersectionFilter<GV>
{
  using ThisType = CustomBoundaryIntersections<GV>;
  using BaseType = IntersectionFilter<GV>;

public:
  using BaseType::logger;
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  /**
   * \attention Takes ownership of boundary_type, do not delete manually!
   */
  explicit CustomBoundaryIntersections(const BoundaryInfo<IntersectionType>& boundary_info,
                                       BoundaryType*&& boundary_type,
                                       const std::string& logging_prefix = "",
                                       const std::array<bool, 3>& logging_state = Common::default_logger_state())
    : BaseType(logging_prefix.empty() ? "xt.grid.customboundaryintersections" : logging_prefix, logging_state)
    , boundary_info_(boundary_info)
    , boundary_type_(boundary_type)
  {}

  explicit CustomBoundaryIntersections(const BoundaryInfo<IntersectionType>& boundary_info,
                                       const std::shared_ptr<BoundaryType>& boundary_type,
                                       const std::string& logging_prefix = "",
                                       const std::array<bool, 3>& logging_state = Common::default_logger_state())
    : BaseType(logging_prefix.empty() ? "xt.grid.customboundaryintersections" : logging_prefix, logging_state)
    , boundary_info_(boundary_info)
    , boundary_type_(boundary_type)
  {}

  CustomBoundaryIntersections(const ThisType& other)
    : BaseType(other)
    , boundary_info_(other.boundary_info_)
    , boundary_type_(other.boundary_type_)
  {}

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new CustomBoundaryIntersections<GridViewType>(*this);
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    LOG_(debug) << "contains(intersection=" << print(intersection)
                << "):\n  boundary_info_.type(intersection) = " << boundary_info_.type(intersection)
                << ", *boundary_type_ = " << *boundary_type_ << ", returning "
                << (boundary_info_.type(intersection) == *boundary_type_) << std::endl;
    return boundary_info_.type(intersection) == *boundary_type_;
  }

protected:
  const BoundaryInfo<IntersectionType>& boundary_info_;
  const std::shared_ptr<BoundaryType> boundary_type_;
}; // class CustomBoundaryIntersections


/**
 * \brief A filter which selects intersections on a given part of the physical boundary and the process boundary.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GV>
class CustomBoundaryAndProcessIntersections : public IntersectionFilter<GV>
{
  using BaseType = IntersectionFilter<GV>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  /**
   * \attention Takes ownership of boundary_type, do not delete manually!
   */
  explicit CustomBoundaryAndProcessIntersections(const BoundaryInfo<IntersectionType>& boundary_info,
                                                 BoundaryType*&& boundary_type)
    : boundary_info_(boundary_info)
    , boundary_type_(boundary_type)
  {}

  explicit CustomBoundaryAndProcessIntersections(const BoundaryInfo<IntersectionType>& boundary_info,
                                                 const std::shared_ptr<BoundaryType>& boundary_type)
    : boundary_info_(boundary_info)
    , boundary_type_(boundary_type)
  {}

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new CustomBoundaryAndProcessIntersections<GridViewType>(boundary_info_, boundary_type_);
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    const bool process_boundary = !intersection.neighbor() && !intersection.boundary();
    const bool physical_boundary = boundary_info_.type(intersection) == *boundary_type_;
    return process_boundary || physical_boundary;
  }

protected:
  const BoundaryInfo<IntersectionType>& boundary_info_;
  const std::shared_ptr<BoundaryType> boundary_type_;
}; // class CustomBoundaryAndProcessIntersections


/**
 * \brief A filter which selects intersections on a the process boundary.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GV>
class ProcessIntersections : public IntersectionFilter<GV>
{
  using BaseType = IntersectionFilter<GV>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit ProcessIntersections() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new ProcessIntersections<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return (!intersection.neighbor() && !intersection.boundary());
  }
}; // class ProcessIntersections


} // namespace Dune::XT::Grid::ApplyOn

#endif // DUNE_XT_GRID_FILTERS_INTERSECTION_HH
