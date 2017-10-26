// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   Rene Milk       (2014 - 2018)
//   Tobias Leibner  (2015 - 2017)

#ifndef DUNE_XT_GRID_WALKER_FILTERS_HH
#define DUNE_XT_GRID_WALKER_FILTERS_HH

#include <functional>

#include <dune/xt/common/memory.hh>

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace internal {


// forwards
template <class GL>
class CombinedEntityFilters;

template <class GL>
class CombinedIntersectionFilters;

template <class GL>
class NegatedEntityFilter;

template <class GL>
class NegatedIntersectionFilter;


} // namespace internal


/**
 *  \brief Interface for entity filters to restrict the range of entities a functor is applied on by the grid walker.
 *  \note  Filters are provided within the ApplyOn namespace.
 *  \note  See also the \sa GridWalker for typical use cases.
 */
template <class GL>
class EntityFilter
{
  static_assert(is_layer<GL>::value, "");

public:
  using GridLayerType = GL;
  using EntityType = extract_entity_t<GridLayerType>;

  virtual ~EntityFilter() = default;

  virtual EntityFilter<GridLayerType>* copy() const = 0;

  virtual bool contains(const GridLayerType& /*grid_layer*/, const EntityType& /*entity*/) const = 0;

  EntityFilter<GridLayerType>* operator!() const
  {
    return new internal::NegatedEntityFilter<GridLayerType>(*this);
  }

  EntityFilter<GridLayerType>* operator&&(const EntityFilter<GridLayerType>& other) const
  {
    return new internal::CombinedEntityFilters<GridLayerType>(
        *this, other, [](const auto& left, const auto& right) { return left && right; });
  }

  EntityFilter<GridLayerType>* operator&&(EntityFilter<GridLayerType>*&& other) const
  {
    return new internal::CombinedEntityFilters<GridLayerType>(
        *this, std::move(other), [](const auto& left, const auto& right) { return left && right; });
  }

  EntityFilter<GridLayerType>* operator||(const EntityFilter<GridLayerType>& other) const
  {
    return new internal::CombinedEntityFilters<GridLayerType>(
        *this, other, [](const auto& left, const auto& right) { return left || right; });
  }

  EntityFilter<GridLayerType>* operator||(EntityFilter<GridLayerType>*&& other) const
  {
    return new internal::CombinedEntityFilters<GridLayerType>(
        *this, std::move(other), [](const auto& left, const auto& right) { return left || right; });
  }
}; // class EntityFilter


/**
 *  \brief Interface for intersection filters to restrict the range of intersections a functor is applied on by the grid
 *         walker.
 *  \note  Filters are provided within the ApplyOn namespace.
 *  \note  See also the \sa GridWalker for typical use cases.
 */
template <class GL>
class IntersectionFilter
{
  static_assert(is_layer<GL>::value, "");

public:
  using GridLayerType = GL;
  using IntersectionType = extract_intersection_t<GridLayerType>;

  virtual ~IntersectionFilter<GridLayerType>() = default;

  virtual IntersectionFilter<GridLayerType>* copy() const = 0;

  virtual bool contains(const GridLayerType& /*grid_layer*/, const IntersectionType& /*intersection*/) const = 0;

  IntersectionFilter<GridLayerType>* operator!() const
  {
    return new internal::NegatedIntersectionFilter<GridLayerType>(*this);
  }

  IntersectionFilter<GridLayerType>* operator&&(const IntersectionFilter<GridLayerType>& other) const
  {
    return new internal::CombinedIntersectionFilters<GridLayerType>(
        *this, other, [](const auto& left, const auto& right) { return left && right; });
  }

  IntersectionFilter<GridLayerType>* operator&&(IntersectionFilter<GridLayerType>*&& other) const
  {
    return new internal::CombinedIntersectionFilters<GridLayerType>(
        *this, std::move(other), [](const auto& left, const auto& right) { return left && right; });
  }

  IntersectionFilter<GridLayerType>* operator||(const IntersectionFilter<GridLayerType>& other) const
  {
    return new internal::CombinedIntersectionFilters<GridLayerType>(
        *this, other, [](const auto& left, const auto& right) { return left || right; });
  }

  IntersectionFilter<GridLayerType>* operator||(IntersectionFilter<GridLayerType>*&& other) const
  {
    return new internal::CombinedIntersectionFilters<GridLayerType>(
        *this, std::move(other), [](const auto& left, const auto& right) { return left || right; });
  }
}; // class IntersectionFilter


namespace ApplyOn {


/**
 *  \brief A filter which selects all entities.
 */
template <class GL>
class AllEntities : public EntityFilter<GL>
{
  using BaseType = EntityFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;

  explicit AllEntities() = default;

  EntityFilter<GridLayerType>* copy() const override final
  {
    return new AllEntities<GridLayerType>();
  }

  bool contains(const GridLayerType& /*grid_layer*/, const EntityType& /*entity*/) const override final
  {
    return true;
  }
}; // class AllEntities


/**
 *  \brief A filter which selects entities which have a boundary intersection.
 */
template <class GL>
class BoundaryEntities : public EntityFilter<GL>
{
  using BaseType = EntityFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;

  explicit BoundaryEntities() = default;

  EntityFilter<GridLayerType>* copy() const override final
  {
    return new BoundaryEntities<GridLayerType>();
  }

  bool contains(const GridLayerType& /*grid_layer*/, const EntityType& entity) const override final
  {
    return entity.hasBoundaryIntersections();
  }
}; // class BoundaryEntities


/**
 *  \brief A filter which selects entities based on a lambda expression.
 */
template <class GL>
class LambdaFilteredEntities : public EntityFilter<GL>
{
  using BaseType = EntityFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using LambdaType = std::function<bool(const GridLayerType&, const EntityType&)>;

  explicit LambdaFilteredEntities(LambdaType lambda)
    : lambda_(lambda)
  {
  }

  EntityFilter<GridLayerType>* copy() const override final
  {
    return new LambdaFilteredEntities<GridLayerType>(lambda_);
  }

  bool contains(const GridLayerType& grid_layer, const EntityType& entity) const override final
  {
    return lambda_(grid_layer, entity);
  }

private:
  const LambdaType lambda_;
}; // class LambdaFilteredEntities


/**
 *  \brief A filter which selects entities in a compatible PartitionSet.
 */
template <class GL, class PartitionSetType>
class PartitionSetEntities : public EntityFilter<GL>
{
  using BaseType = EntityFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  typedef typename BaseType::EntityType EntityType;

  explicit PartitionSetEntities() = default;

  EntityFilter<GridLayerType>* copy() const override final
  {
    return new PartitionSetEntities<GridLayerType, PartitionSetType>();
  }

  bool contains(const GridLayerType& /*grid_layer*/, const EntityType& entity) const override final
  {
    return PartitionSetType::contains(entity.partitionType());
  }
}; // class PartitionSetEntities


/**
 *  \brief A filter which selects all intersections.
 */
template <class GL>
class AllIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  explicit AllIntersections() = default;

  IntersectionFilter<GridLayerType>* copy() const override final
  {
    return new AllIntersections<GridLayerType>();
  }

  bool contains(const GridLayerType& /*grid_layer*/, const IntersectionType& /*intersection*/) const override final
  {
    return true;
  }
}; // class AllIntersections


/**
 *  \brief A filter which selects each inner intersection.
 *
 *  \note  To decide if this in an inner intersection,
\code
intersection.neighbor() && !intersection.boundary()
\endcode
 *         is used.
 */
template <class GL>
class InnerIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  explicit InnerIntersections() = default;

  IntersectionFilter<GridLayerType>* copy() const override final
  {
    return new InnerIntersections<GridLayerType>();
  }

  bool contains(const GridLayerType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return intersection.neighbor() && !intersection.boundary();
  }
}; // class InnerIntersections


/**
 *  \brief A filter which selects each inner intersection only once.
 *
 *  \note  To decide if this in an inner intersection,
\code
intersection.neighbor() && !intersection.boundary()
\endcode
 *         is used, and true is returned, if the index of the inside() entity is smaller than the index of the outside()
 *         entity.
 */
template <class GL>
class InnerIntersectionsOnce : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  explicit InnerIntersectionsOnce() = default;

  IntersectionFilter<GridLayerType>* copy() const override final
  {
    return new InnerIntersectionsOnce<GridLayerType>();
  }

  bool contains(const GridLayerType& grid_layer, const IntersectionType& intersection) const override final
  {
    if (intersection.neighbor() && !intersection.boundary()) {
      const auto insideEntity = intersection.inside();
      const auto outsideNeighbor = intersection.outside();
      return grid_layer.indexSet().index(insideEntity) < grid_layer.indexSet().index(outsideNeighbor);
    } else
      return false;
  }
}; // class InnerIntersectionsOnce


/**
 *  \brief Selects each inner intersection in given partition only once.
 *  \see InnerIntersectionsOnce
 */
template <class GL, class PartitionSetType>
class PartitionSetInnerIntersectionsOnce : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  explicit PartitionSetInnerIntersectionsOnce() = default;

  IntersectionFilter<GridLayerType>* copy() const override final
  {
    return new PartitionSetInnerIntersectionsOnce<GridLayerType, PartitionSetType>();
  }

  bool apply_on(const GridLayerType& grid_layer, const IntersectionType& intersection) const override final
  {
    if (intersection.neighbor() && !intersection.boundary()
        && PartitionSetType::contains(intersection.inside().partitionType())) {
      const auto insideEntity = intersection.inside();
      const auto outsideNeighbor = intersection.outside();
      if (!PartitionSetType::contains(intersection.outside().partitionType()))
        return true;
      return grid_layer.indexSet().index(insideEntity) < grid_layer.indexSet().index(outsideNeighbor);
    } else
      return false;
  }
}; // class PartitionSetInnerIntersectionsOnce


/**
 *  \brief A filter which selects boundary intersections.
 *
 */
template <class GL>
class BoundaryIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  explicit BoundaryIntersections() = default;

  IntersectionFilter<GridLayerType>* copy() const override final
  {
    return new BoundaryIntersections<GridLayerType>();
  }

  bool contains(const GridLayerType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return intersection.boundary();
  }
}; // class BoundaryIntersections


/**
 *  \brief A filter which selects intersections on a non-periodic boundary.
 *
 */
template <class GL>
class NonPeriodicBoundaryIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  explicit NonPeriodicBoundaryIntersections() = default;

  IntersectionFilter<GridLayerType>* copy() const override final
  {
    return new NonPeriodicBoundaryIntersections<GridLayerType>();
  }

  bool contains(const GridLayerType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return intersection.boundary() && !intersection.neighbor();
  }
}; // class BoundaryIntersections


/**
 *  \brief A filter which selects intersections on a periodic boundary only once.
 *
 *         To decide if this a periodic intersection,
\code
intersection.neighbor() && intersection.boundary()
\endcode
 *         is used.
 */
template <class GL>
class PeriodicBoundaryIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  explicit PeriodicBoundaryIntersections() = default;

  IntersectionFilter<GridLayerType>* copy() const override final
  {
    return new PeriodicBoundaryIntersections<GridLayerType>();
  }

  bool contains(const GridLayerType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return intersection.neighbor() && intersection.boundary();
  }
}; // class PeriodicBoundaryIntersections


/**
 *  \brief A filter which selects intersections on a periodic boundary only once.
 *
 *         To decide if this in an periodic intersection,
\code
intersection.neighbor() && intersection.boundary()
\endcode
 *         is used, and true is returned, if the index of the inside() entity is smaller than the index of the outside()
 *         entity.
 */
template <class GL>
class PeriodicBoundaryIntersectionsOnce : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  explicit PeriodicBoundaryIntersectionsOnce() = default;

  IntersectionFilter<GridLayerType>* copy() const override final
  {
    return new PeriodicBoundaryIntersectionsOnce<GridLayerType>();
  }

  bool contains(const GridLayerType& grid_layer, const IntersectionType& intersection) const override final
  {
    if (intersection.neighbor() && intersection.boundary()) {
      const auto insideEntity = intersection.inside();
      const auto outsideNeighbor = intersection.outside();
      return grid_layer.indexSet().index(insideEntity) < grid_layer.indexSet().index(outsideNeighbor);
    } else {
      return false;
    }
  }
}; // class PeriodicBoundaryIntersectionsOnce


/**
 *  \brief A filter which selects intersections based on a lambda expression.
 *
 */
template <class GL>
class LambdaFilteredIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;
  using LambdaType = std::function<bool(const GridLayerType&, const IntersectionType&)>;

  explicit LambdaFilteredIntersections(LambdaType lambda)
    : lambda_(lambda)
  {
  }

  virtual IntersectionFilter<GridLayerType>* copy() const override final
  {
    return new LambdaFilteredIntersections<GridLayerType>(lambda_);
  }

  bool contains(const GridLayerType& grid_layer, const IntersectionType& intersection) const override final
  {
    return lambda_(grid_layer, intersection);
  }

private:
  const LambdaType lambda_;
}; // class BoundaryIntersections


/**
 *  \brief A filter which selects intersections on a given part of the boundary.
 *
 */
template <class GL>
class CustomBoundaryIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  /**
   * \attention Takes ownership of boundary_type, do not delete manually!
   */
  explicit CustomBoundaryIntersections(const BoundaryInfo<IntersectionType>& boundary_info,
                                       BoundaryType*&& boundary_type)
    : boundary_info_(boundary_info)
    , boundary_type_(boundary_type)
  {
  }

  explicit CustomBoundaryIntersections(const BoundaryInfo<IntersectionType>& boundary_info,
                                       const std::shared_ptr<BoundaryType>& boundary_type)
    : boundary_info_(boundary_info)
    , boundary_type_(boundary_type)
  {
  }

  IntersectionFilter<GridLayerType>* copy() const override final
  {
    return new CustomBoundaryIntersections<GridLayerType>(boundary_info_, boundary_type_);
  }

  bool contains(const GridLayerType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return boundary_info_.type(intersection) == *boundary_type_;
  }

protected:
  const BoundaryInfo<IntersectionType>& boundary_info_;
  const std::shared_ptr<BoundaryType> boundary_type_;
}; // class CustomBoundaryIntersections


/**
 *  \brief A filter which selects intersections on a given part of the physical boundary and the process boundary.
 *
 */
template <class GL>
class CustomBoundaryAndProcessIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  /**
   * \attention Takes ownership of boundary_type, do not delete manually!
   */
  explicit CustomBoundaryAndProcessIntersections(const BoundaryInfo<IntersectionType>& boundary_info,
                                                 BoundaryType*&& boundary_type)
    : boundary_info_(boundary_info)
    , boundary_type_(std::move(boundary_type))
  {
  }

  explicit CustomBoundaryAndProcessIntersections(const BoundaryInfo<IntersectionType>& boundary_info,
                                                 const std::shared_ptr<BoundaryType>& boundary_type)
    : boundary_info_(boundary_info)
    , boundary_type_(boundary_type)
  {
  }

  IntersectionFilter<GridLayerType>* copy() const override final
  {
    return new CustomBoundaryAndProcessIntersections<GridLayerType>(boundary_info_, boundary_type_);
  }

  bool contains(const GridLayerType& /*grid_layer*/, const IntersectionType& intersection) const override final
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
 *  \brief A filter which selects intersections on a the process boundary.
 *
 */
template <class GL>
class ProcessIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  explicit ProcessIntersections() = default;

  IntersectionFilter<GridLayerType>* copy() const override final
  {
    return new ProcessIntersections<GridLayerType>();
  }

  bool contains(const GridLayerType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return (!intersection.neighbor() && !intersection.boundary());
  }
}; // class ProcessIntersections


} // namespace ApplyOn
namespace internal {


template <class GL>
class CombinedIntersectionFilters : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  CombinedIntersectionFilters(const BaseType& left,
                              const BaseType& right,
                              std::function<bool(const bool&, const bool&)> combine_lambda)
    : left_(left.copy())
    , right_(right.copy())
    , combine_lambda_(combine_lambda)
  {
  }

  CombinedIntersectionFilters(const BaseType& left,
                              BaseType*&& right,
                              std::function<bool(const bool&, const bool&)> combine_lambda)
    : left_(left.copy())
    , right_(std::move(right))
    , combine_lambda_(combine_lambda)
  {
  }

  IntersectionFilter<GridLayerType>* copy() const override final
  {
    return new CombinedIntersectionFilters<GridLayerType>(*left_, *right_, combine_lambda_);
  }

  bool contains(const GridLayerType& grid_layer, const IntersectionType& intersection) const override final
  {
    return combine_lambda_(left_->contains(grid_layer, intersection), right_->contains(grid_layer, intersection));
  }

private:
  const std::unique_ptr<BaseType> left_;
  const std::unique_ptr<BaseType> right_;
  const std::function<bool(const bool&, const bool&)> combine_lambda_;
}; // class CombinedIntersectionFilters


template <class GL>
class CombinedEntityFilters : public EntityFilter<GL>
{
  using BaseType = EntityFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;

  CombinedEntityFilters(const BaseType& left,
                        const BaseType& right,
                        std::function<bool(const bool&, const bool&)> combine_lambda)
    : left_(left.copy())
    , right_(right.copy())
    , combine_lambda_(combine_lambda)
  {
  }

  CombinedEntityFilters(const BaseType& left,
                        BaseType*&& right,
                        std::function<bool(const bool&, const bool&)> combine_lambda)
    : left_(left.copy())
    , right_(std::move(right))
    , combine_lambda_(combine_lambda)
  {
  }

  EntityFilter<GridLayerType>* copy() const override final
  {
    return new CombinedEntityFilters<GridLayerType>(*left_, *right_, combine_lambda_);
  }

  bool contains(const GridLayerType& grid_layer, const EntityType& entity) const override final
  {
    return combine_lambda_(left_->contains(grid_layer, entity), right_->contains(grid_layer, entity));
  }

private:
  const std::unique_ptr<BaseType> left_;
  const std::unique_ptr<BaseType> right_;
  const std::function<bool(const bool&, const bool&)> combine_lambda_;
}; // class CombinedEntityFilters


template <class GL>
class NegatedIntersectionFilter : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  NegatedIntersectionFilter(const BaseType& filter)
    : filter_(filter.copy())
  {
  }

  NegatedIntersectionFilter(BaseType*&& filter)
    : filter_(std::move(filter))
  {
  }

  IntersectionFilter<GridLayerType>* copy() const override final
  {
    return new NegatedIntersectionFilter<GridLayerType>(*filter_);
  }

  bool contains(const GridLayerType& grid_layer, const IntersectionType& intersection) const override final
  {
    return !filter_->contains(grid_layer, intersection);
  }

private:
  const std::unique_ptr<BaseType> filter_;
}; // class NegatedIntersectionFilter


template <class GL>
class NegatedEntityFilter : public EntityFilter<GL>
{
  using BaseType = EntityFilter<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;

  NegatedEntityFilter(const BaseType& filter)
    : filter_(filter.copy())
  {
  }

  NegatedEntityFilter(BaseType*&& filter)
    : filter_(std::move(filter))
  {
  }

  EntityFilter<GridLayerType>* copy() const override final
  {
    return new NegatedEntityFilter<GridLayerType>(*filter_);
  }

  bool contains(const GridLayerType& grid_layer, const EntityType& entity) const override final
  {
    return !filter_->contains(grid_layer, entity);
  }

private:
  const std::unique_ptr<BaseType> filter_;
}; // class NegatedEntityFilter


} // namespace internal
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_WALKER_FILTERS_HH
