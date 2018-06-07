// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   Rene Milk       (2014 - 2018)
//   Tobias Leibner  (2015 - 2018)

#ifndef DUNE_XT_GRID_WALKER_APPLY_ON_HH
#define DUNE_XT_GRID_WALKER_APPLY_ON_HH

#include <functional>

#include <dune/grid/common/partitionset.hh>

#include <dune/xt/common/memory.hh>

#include <dune/xt/grid/boundaryinfo.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace ApplyOn {
namespace internal {


// forwards
template <class GL>
class CombinedIntersectionFilters;

template <class GL>
class NegatedIntersectionFilter;


} // namespace internal


/**
 *  \brief Interface for functors to select on which entities to apply.
 */
template <class GridLayerImp>
class WhichEntity
{
public:
  typedef GridLayerImp GridLayerType;
  using EntityType = extract_entity_t<GridLayerType>;

  virtual ~WhichEntity()
  {
  }

  virtual WhichEntity<GridLayerImp>* copy() const = 0; // required for the python bindings

  virtual bool apply_on(const GridLayerType& /*grid_layer*/, const EntityType& /*entity*/) const = 0;
}; // class WhichEntity


/**
 *  \brief Selects all entities.
 */
template <class GridLayerImp>
class AllEntities : public WhichEntity<GridLayerImp>
{
  typedef WhichEntity<GridLayerImp> BaseType;

public:
  using typename BaseType::GridLayerType;
  typedef typename BaseType::EntityType EntityType;

  virtual WhichEntity<GridLayerImp>* copy() const override final
  {
    return new AllEntities<GridLayerImp>();
  }

  virtual bool apply_on(const GridLayerType& /*grid_layer*/, const EntityType& /*entity*/) const override final
  {
    return true;
  }
}; // class AllEntities

template <class GridLayerImp>
const AllEntities<GridLayerImp>* all_entities()
{
  static const AllEntities<GridLayerImp> all_entities_singleton;
  return &all_entities_singleton;
}


/**
 *  \brief Selects entities which have a boundary intersection.
 */
template <class GridLayerImp>
class BoundaryEntities : public WhichEntity<GridLayerImp>
{
  typedef WhichEntity<GridLayerImp> BaseType;

public:
  using typename BaseType::GridLayerType;
  typedef typename BaseType::EntityType EntityType;

  virtual WhichEntity<GridLayerImp>* copy() const override final
  {
    return new BoundaryEntities<GridLayerImp>();
  }

  virtual bool apply_on(const GridLayerType& /*grid_layer*/, const EntityType& entity) const override final
  {
    return entity.hasBoundaryIntersections();
  }
}; // class BoundaryEntities


/**
 *  \brief Interface for functors to select on which intersections to apply.
 *
 *  \note  When deriving from this class, one can use the \sa internal::WhichIntersectionBase to provide a generic
 *         implementation of WhichIntersection::copy in two circumstances: (i) if one requires a BoundaryInfo, derive
 *         from internal::WhichIntersectionBase<..., true>, \sa DirichletIntersections; (ii) if the derived class does
 *         not have any members, derive from internal::WhichIntersectionBase<...>, \sa InnerIntersections. In all other
 *         circumstances, one has to manually implement copy, \sa FilteredIntersections.
 */
template <class GridLayerImp>
class WhichIntersection
{
public:
  typedef GridLayerImp GridLayerType;
  using IntersectionType = extract_intersection_t<GridLayerType>;

  virtual ~WhichIntersection<GridLayerImp>()
  {
  }

  virtual WhichIntersection<GridLayerImp>* copy() const = 0; // required for redirect lambdas, i.e., in python bindings

  virtual bool apply_on(const GridLayerType& /*grid_layer*/, const IntersectionType& /*intersection*/) const = 0;

  WhichIntersection<GridLayerType>* operator!() const
  {
    return new internal::NegatedIntersectionFilter<GridLayerType>(*this);
  }

  WhichIntersection<GridLayerType>* operator&&(const WhichIntersection<GridLayerType>& other) const
  {
    return new internal::CombinedIntersectionFilters<GridLayerType>(
        *this, other, [](const auto& left, const auto& right) { return left && right; });
  }

  WhichIntersection<GridLayerType>* operator&&(WhichIntersection<GridLayerType>*&& other) const
  {
    return new internal::CombinedIntersectionFilters<GridLayerType>(
        *this, std::move(other), [](const auto& left, const auto& right) { return left && right; });
  }

  WhichIntersection<GridLayerType>* operator||(const WhichIntersection<GridLayerType>& other) const
  {
    return new internal::CombinedIntersectionFilters<GridLayerType>(
        *this, other, [](const auto& left, const auto& right) { return left || right; });
  }

  WhichIntersection<GridLayerType>* operator||(WhichIntersection<GridLayerType>*&& other) const
  {
    return new internal::CombinedIntersectionFilters<GridLayerType>(
        *this, std::move(other), [](const auto& left, const auto& right) { return left || right; });
  }
}; // class WhichIntersection


namespace internal {


template <class GL>
class CombinedIntersectionFilters : public WhichIntersection<GL>
{
  using BaseType = WhichIntersection<GL>;

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

  WhichIntersection<GridLayerType>* copy() const override final
  {
    return new CombinedIntersectionFilters<GL>(*left_, *right_, combine_lambda_);
  }

  bool apply_on(const GridLayerType& grid_layer, const IntersectionType& intersection) const override final
  {
    return combine_lambda_(left_->apply_on(grid_layer, intersection), right_->apply_on(grid_layer, intersection));
  }

private:
  const std::unique_ptr<BaseType> left_;
  const std::unique_ptr<BaseType> right_;
  const std::function<bool(const bool&, const bool&)> combine_lambda_;
}; // class CombinedIntersectionFilters


template <class GL>
class NegatedIntersectionFilter : public WhichIntersection<GL>
{
  using BaseType = WhichIntersection<GL>;

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

  WhichIntersection<GridLayerType>* copy() const override final
  {
    return new NegatedIntersectionFilter<GL>(*filter_);
  }

  bool apply_on(const GridLayerType& grid_layer, const IntersectionType& intersection) const override final
  {
    return !filter_->apply_on(grid_layer, intersection);
  }

private:
  const std::unique_ptr<BaseType> filter_;
}; // class NegatedIntersectionFilter


template <class GV, class Imp, bool ctor_with_boundary_info = false>
class WhichIntersectionBase : public WhichIntersection<GV>
{
public:
  virtual WhichIntersection<GV>* copy() const override final
  {
    return new Imp();
  }
};

template <class GV, class Imp>
class WhichIntersectionBase<GV, Imp, true> : public WhichIntersection<GV>
{
public:
  using typename WhichIntersection<GV>::IntersectionType;

  WhichIntersectionBase(const BoundaryInfo<IntersectionType>& boundary_info)
    : boundary_info_(boundary_info)
  {
  }

  virtual WhichIntersection<GV>* copy() const override final
  {
    return new Imp(boundary_info_.access());
  }

protected:
  Common::ConstStorageProvider<BoundaryInfo<IntersectionType>> boundary_info_;
};


} // namespace internal


/**
 *  \brief Selects all intersections.
 */
template <class GridLayerImp>
class AllIntersections : public internal::WhichIntersectionBase<GridLayerImp, AllIntersections<GridLayerImp>>
{
  typedef WhichIntersection<GridLayerImp> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  virtual bool apply_on(const GridLayerType& /*grid_layer*/,
                        const IntersectionType& /*intersection*/) const override final
  {
    return true;
  }
}; // class AllIntersections

template <class GridLayerImp>
const AllIntersections<GridLayerImp>* all_intersections()
{
  static const AllIntersections<GridLayerImp> all_intersections_singleton;
  return &all_intersections_singleton;
}


/**
 *  \brief Selects no intersections.
 */
template <class GridLayerImp>
class NoIntersections : public internal::WhichIntersectionBase<GridLayerImp, NoIntersections<GridLayerImp>>
{
  typedef WhichIntersection<GridLayerImp> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  virtual bool apply_on(const GridLayerType& /*grid_layer*/,
                        const IntersectionType& /*intersection*/) const override final
  {
    return false;
  }
}; // class NoIntersections


/**
 *  \brief Selects each inner intersection.
 *
 *  To decide if this in an inner intersection,
\code
intersection.neighbor() && !intersection.boundary()
\endcode
 *  is used.
 */
template <class GridLayerImp>
class InnerIntersections : public internal::WhichIntersectionBase<GridLayerImp, InnerIntersections<GridLayerImp>>
{
  typedef WhichIntersection<GridLayerImp> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  virtual bool apply_on(const GridLayerType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return intersection.neighbor() && !intersection.boundary();
  }
}; // class InnerIntersections


/**
 *  \brief Selects each inner intersection (in given partition) only once.
 */
template <class GridLayerImp, class PartitionSetType = Dune::Partitions::All>
class InnerIntersectionsPrimally
    : public internal::WhichIntersectionBase<GridLayerImp, InnerIntersectionsPrimally<GridLayerImp, PartitionSetType>>
{
  typedef internal::WhichIntersectionBase<GridLayerImp, InnerIntersectionsPrimally> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  virtual bool apply_on(const GridLayerType& grid_layer, const IntersectionType& intersection) const override final
  {
    if (intersection.neighbor() && !intersection.boundary()
        && PartitionSetType::contains(intersection.inside().partitionType())) {
      const auto insideEntity = intersection.inside();
      const auto outsideNeighbor = intersection.outside();
      // if outside entity is not in Partition, we need to use this intersection even if index is greater
      if (!PartitionSetType::contains(intersection.outside().partitionType()))
        return true;
      return grid_layer.indexSet().index(insideEntity) < grid_layer.indexSet().index(outsideNeighbor);
    } else
      return false;
  }
}; // class InnerIntersectionsPrimally


template <class GridLayerImp>
class BoundaryIntersections : public internal::WhichIntersectionBase<GridLayerImp, BoundaryIntersections<GridLayerImp>>
{
  typedef internal::WhichIntersectionBase<GridLayerImp, BoundaryIntersections> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  virtual bool apply_on(const GridLayerType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return intersection.boundary();
  }
}; // class BoundaryIntersections


template <class GridLayerImp, class PartitionSetType = Dune::Partitions::All>
class NonPeriodicBoundaryIntersections
    : public internal::WhichIntersectionBase<GridLayerImp,
                                             NonPeriodicBoundaryIntersections<GridLayerImp, PartitionSetType>>
{
  typedef internal::WhichIntersectionBase<GridLayerImp,
                                          NonPeriodicBoundaryIntersections<GridLayerImp, PartitionSetType>>
      BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  virtual bool apply_on(const GridLayerType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return intersection.boundary() && !intersection.neighbor()
           && PartitionSetType::contains(intersection.inside().partitionType());
  }
}; // class NonPeriodicBoundaryIntersections


/**
 *  \brief Selects each periodic intersection.
 *
 *  To decide if this in an periodic intersection,
\code
intersection.neighbor() && intersection.boundary()
\endcode
 *  is used.
 */
template <class GridLayerImp, class PartitionSetType = Dune::Partitions::All>
class PeriodicIntersections
    : public internal::WhichIntersectionBase<GridLayerImp, PeriodicIntersections<GridLayerImp, PartitionSetType>>
{
  typedef internal::WhichIntersectionBase<GridLayerImp, PeriodicIntersections> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  virtual bool apply_on(const GridLayerType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return intersection.neighbor() && intersection.boundary()
           && PartitionSetType::contains(intersection.inside().partitionType());
  }
}; // class PeriodicIntersections

template <class GridLayerImp, class PartitionSetType = Dune::Partitions::All>
class PeriodicIntersectionsPrimally
    : public internal::WhichIntersectionBase<GridLayerImp,
                                             PeriodicIntersectionsPrimally<GridLayerImp, PartitionSetType>>
{
  typedef internal::WhichIntersectionBase<GridLayerImp, PeriodicIntersectionsPrimally> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  virtual bool apply_on(const GridLayerType& grid_layer, const IntersectionType& intersection) const override final
  {
    if (intersection.neighbor() && intersection.boundary()
        && PartitionSetType::contains(intersection.inside().partitionType())) {
      const auto insideEntity = intersection.inside();
      const auto outsideNeighbor = intersection.outside();
      if (!PartitionSetType::contains(intersection.outside().partitionType()))
        return true;
      return grid_layer.indexSet().index(insideEntity) < grid_layer.indexSet().index(outsideNeighbor);
    } else {
      return false;
    }
  }
}; // class PeriodicIntersectionsPrimally


template <class GridLayerImp>
class FilteredIntersections : public WhichIntersection<GridLayerImp>
{
  typedef WhichIntersection<GridLayerImp> BaseType;
  typedef FilteredIntersections<GridLayerImp> ThisType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;
  typedef std::function<bool(const GridLayerType&, const IntersectionType&)> FilterType;

  FilteredIntersections(FilterType filter)
    : filter_(filter)
  {
  }

  virtual WhichIntersection<GridLayerImp>* copy() const override final
  {
    return new FilteredIntersections<GridLayerImp>(filter_);
  }

  virtual bool apply_on(const GridLayerType& grid_layer, const IntersectionType& intersection) const override final
  {
    return filter_(grid_layer, intersection);
  }

private:
  const FilterType filter_;
}; // class BoundaryIntersections


template <class GridLayerImp>
class CustomBoundaryIntersections : public WhichIntersection<GridLayerImp>
{
  typedef WhichIntersection<GridLayerImp> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  /**
   * \attention Takes ownership of boundary_type, do not delete manually!
   */
  explicit CustomBoundaryIntersections(const BoundaryInfo<IntersectionType>& boundary_info,
                                       BoundaryType*&& boundary_type)
    : boundary_info_(boundary_info)
    , boundary_type_(std::move(boundary_type))
  {
  }

  explicit CustomBoundaryIntersections(const BoundaryInfo<IntersectionType>& boundary_info,
                                       const std::shared_ptr<BoundaryType>& boundary_type)
    : boundary_info_(boundary_info)
    , boundary_type_(boundary_type)
  {
  }

  WhichIntersection<GridLayerType>* copy() const override final
  {
    return new CustomBoundaryIntersections<GridLayerType>(boundary_info_, boundary_type_);
  }

  bool apply_on(const GridLayerType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return boundary_info_.type(intersection) == *boundary_type_;
  }

protected:
  const BoundaryInfo<IntersectionType>& boundary_info_;
  const std::shared_ptr<BoundaryType> boundary_type_;
}; // class CustomBoundaryIntersections


template <class GridLayerImp>
class DirichletIntersections
    : public internal::WhichIntersectionBase<GridLayerImp, DirichletIntersections<GridLayerImp>, true>
{
  typedef internal::WhichIntersectionBase<GridLayerImp, DirichletIntersections<GridLayerImp>, true> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  explicit DirichletIntersections(const BoundaryInfo<IntersectionType>& boundary_info)
    : BaseType(boundary_info)
  {
  }

  bool apply_on(const GridLayerType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    static constexpr const XT::Grid::DirichletBoundary dirichlet{};
    return boundary_info_.access().type(intersection) == dirichlet;
  }

protected:
  using BaseType::boundary_info_;
}; // class DirichletIntersections

template <class GridLayerImp>
class DirichletAndProcessIntersections
    : public internal::WhichIntersectionBase<GridLayerImp, DirichletAndProcessIntersections<GridLayerImp>, true>
{
  typedef internal::WhichIntersectionBase<GridLayerImp, DirichletAndProcessIntersections<GridLayerImp>, true> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  explicit DirichletAndProcessIntersections(const BoundaryInfo<IntersectionType>& boundary_info)
    : BaseType(boundary_info)
  {
  }

  virtual bool apply_on(const GridLayerType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    static constexpr const XT::Grid::DirichletBoundary dirichlet{};
    const bool process_boundary = (!intersection.neighbor() && !intersection.boundary());
    return (boundary_info_.access().type(intersection) == dirichlet) || process_boundary;
  }

protected:
  using BaseType::boundary_info_;
}; // class DirichletAndProcessIntersections

template <class GridLayerImp>
class ProcessIntersections : public internal::WhichIntersectionBase<GridLayerImp, ProcessIntersections<GridLayerImp>>
{
  typedef internal::WhichIntersectionBase<GridLayerImp, ProcessIntersections<GridLayerImp>> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  virtual bool apply_on(const GridLayerType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return (!intersection.neighbor() && !intersection.boundary());
  }

}; // class ProcessIntersections


template <class GridLayerImp>
class NeumannIntersections
    : public internal::WhichIntersectionBase<GridLayerImp, NeumannIntersections<GridLayerImp>, true>
{
  typedef internal::WhichIntersectionBase<GridLayerImp, NeumannIntersections<GridLayerImp>, true> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  explicit NeumannIntersections(const BoundaryInfo<IntersectionType>& boundary_info)
    : BaseType(boundary_info)
  {
  }

  bool apply_on(const GridLayerType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return boundary_info_.access().type(intersection) == NeumannBoundary();
  }

protected:
  using BaseType::boundary_info_;
}; // class NeumannIntersections

template <class GridLayerImp>
class ReflectingIntersections
    : public internal::WhichIntersectionBase<GridLayerImp, ReflectingIntersections<GridLayerImp>, true>
{
  typedef internal::WhichIntersectionBase<GridLayerImp, ReflectingIntersections<GridLayerImp>, true> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  explicit ReflectingIntersections(const BoundaryInfo<IntersectionType>& boundary_info)
    : BaseType(boundary_info)
  {
  }

  virtual bool apply_on(const GridLayerType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return boundary_info_.access().type(intersection) == ReflectingBoundary();
  }

protected:
  using BaseType::boundary_info_;
}; // class ReflectingIntersections

/**
 *  \brief Selects entities in the compatible PartitionSet.
 */
template <class GridLayerImp, class PartitionSetType>
class PartitionSetEntities : public WhichEntity<GridLayerImp>
{
  typedef WhichEntity<GridLayerImp> BaseType;

public:
  using typename BaseType::GridLayerType;
  typedef typename BaseType::EntityType EntityType;

  virtual WhichEntity<GridLayerImp>* copy() const override final
  {
    return new PartitionSetEntities<GridLayerImp, PartitionSetType>();
  }

  virtual bool apply_on(const GridLayerType& /*grid_layer*/, const EntityType& entity) const override final
  {
    return PartitionSetType::contains(entity.partitionType());
  }
}; // class PartitionSetEntities
} // namespace ApplyOn
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_WALKER_APPLY_ON_HH
