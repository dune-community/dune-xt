// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014 - 2018)
//   Tobias Leibner  (2014 - 2015, 2017)

#ifndef DUNE_XT_GRID_WALKER_HH
#define DUNE_XT_GRID_WALKER_HH

#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

#if HAVE_TBB
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#include <tbb/tbb_stddef.h>
#endif

#include <dune/common/version.hh>

#include <dune/grid/common/rangegenerators.hh>

#if HAVE_TBB
#include <dune/xt/grid/parallel/partitioning/ranged.hh>
#endif
#include <dune/xt/common/parallel/threadmanager.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/unused.hh>
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/grid/functors/lambda.hh>
#include <dune/xt/grid/type_traits.hh>

#include "walker/filters.hh"
#include "walker/wrapper.hh"

namespace Dune {
namespace XT {
namespace Grid {


template <class GL>
class Walker : public EntityAndIntersectionFunctor<GL>
{
  using BaseType = EntityAndIntersectionFunctor<GL>;
  using ThisType = Walker<GL>;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;

  explicit Walker(GridLayerType grd_lr)
    : grid_layer_(grd_lr)
  {
  }

  Walker(const Walker& other) = delete;
  Walker(Walker&& source) = default;

  const GridLayerType& grid_layer() const
  {
    return grid_layer_;
  }

  GridLayerType& grid_layer()
  {
    return grid_layer_;
  }

  /**
   * \name These methods can be used to append an \sa EntityFunctor.
   * \{
   */

  ThisType& append(EntityFunctor<GL>& functor, EntityFilter<GL>*&& filter = new ApplyOn::AllEntities<GL>())
  {
    entity_functor_wrappers_.emplace_back(new internal::EntityFunctorWrapper<GL>(functor, std::move(filter)));
    return *this;
  }

  ThisType& append(EntityFunctor<GL>& functor, std::function<bool(const GL&, const EntityType&)> entity_filter)
  {
    entity_functor_wrappers_.emplace_back(
        new internal::EntityFunctorWrapper<GL>(functor, new ApplyOn::LambdaFilteredEntities<GL>(entity_filter)));
    return *this;
  }

  ThisType& append(EntityFunctor<GL>*&& functor, EntityFilter<GL>*&& filter)
  {
    entity_functor_wrappers_.emplace_back(
        new internal::EntityFunctorWrapper<GL>(std::move(functor), std::move(filter)));
    return *this;
  }

  /**
   * \}
   * \name These methods can be used to append an entity lambda expression.
   * \{
   */

  ThisType& append(std::function<void(const EntityType&)> apply_lambda,
                   EntityFilter<GL>*&& filter = new ApplyOn::AllEntities<GL>())
  {
    return this->append(new EntityLambdaFunctor<GL>(apply_lambda), std::move(filter));
  }

  ThisType& append(std::function<void(const EntityType&)> apply_lambda,
                   std::function<bool(const GL&, const EntityType&)> filter)
  {
    return this->append(new EntityLambdaFunctor<GL>(apply_lambda), new ApplyOn::LambdaFilteredEntities<GL>(filter));
  }

  ThisType& append(std::function<void(const EntityType&)> apply_lambda,
                   std::function<void()> prepare_lambda,
                   std::function<void()> finalize_lambda,
                   EntityFilter<GL>*&& filter = new ApplyOn::AllEntities<GL>())
  {
    return this->append(new EntityLambdaFunctor<GL>(apply_lambda, prepare_lambda, finalize_lambda), std::move(filter));
  }

  ThisType& append(std::function<void(const EntityType&)> apply_lambda,
                   std::function<void()> prepare_lambda,
                   std::function<void()> finalize_lambda,
                   std::function<bool(const GL&, const EntityType&)> filter)
  {
    return this->append(new EntityLambdaFunctor<GL>(apply_lambda, prepare_lambda, finalize_lambda),
                        new ApplyOn::LambdaFilteredEntities<GL>(filter));
  }

  /**
   * \}
   * \name These methods can be used to append an \sa IntersectionFunctor.
   * \{
   */

  ThisType& append(IntersectionFunctor<GL>& functor,
                   IntersectionFilter<GL>*&& filter = new ApplyOn::AllIntersections<GL>())
  {
    intersection_functor_wrappers_.emplace_back(
        new internal::IntersectionFunctorWrapper<GL>(functor, std::move(filter)));
    return *this;
  }

  ThisType& append(IntersectionFunctor<GL>& functor, std::function<bool(const GL&, const IntersectionType&)> filter)
  {
    intersection_functor_wrappers_.emplace_back(
        new internal::IntersectionFunctorWrapper<GL>(functor, new ApplyOn::LambdaFilteredIntersections<GL>(filter)));
    return *this;
  }

  ThisType& append(IntersectionFunctor<GL>*&& functor, IntersectionFilter<GL>*&& filter)
  {
    intersection_functor_wrappers_.emplace_back(
        new internal::IntersectionFunctorWrapper<GL>(std::move(functor), std::move(filter)));
    return *this;
  }

  /**
   * \}
   * \name These methods can be used to append an intersection lambda expression.
   * \{
   */

  ThisType& append(std::function<void(const IntersectionType&, const EntityType&, const EntityType&)> apply_lambda,
                   IntersectionFilter<GL>*&& filter = new ApplyOn::AllIntersections<GL>())
  {
    return this->append(new IntersectionLambdaFunctor<GL>(apply_lambda), std::move(filter));
  }

  ThisType& append(std::function<void(const IntersectionType&, const EntityType&, const EntityType&)> apply_lambda,
                   std::function<bool(const GL&, const IntersectionType&)> filter)
  {
    return this->append(new IntersectionLambdaFunctor<GL>(apply_lambda),
                        new ApplyOn::LambdaFilteredIntersections<GL>(filter));
  }

  ThisType& append(std::function<void(const IntersectionType&, const EntityType&, const EntityType&)> apply_lambda,
                   std::function<void()> prepare_lambda,
                   std::function<void()> finalize_lambda,
                   IntersectionFilter<GL>*&& filter = new ApplyOn::AllIntersections<GL>())
  {
    return this->append(new IntersectionLambdaFunctor<GL>(apply_lambda, prepare_lambda, finalize_lambda),
                        std::move(filter));
  }

  ThisType& append(std::function<void(const IntersectionType&, const EntityType&, const EntityType&)> apply_lambda,
                   std::function<void()> prepare_lambda,
                   std::function<void()> finalize_lambda,
                   std::function<bool(const GL&, const IntersectionType&)> filter)
  {
    return this->append(new IntersectionLambdaFunctor<GL>(apply_lambda, prepare_lambda, finalize_lambda),
                        new ApplyOn::LambdaFilteredIntersections<GL>(filter));
  }

  /**
   * \}
   * \name These methods can be used to append an \sa EntityAndIntersectionFunctor.
   * \{
   */

  ThisType& append(EntityAndIntersectionFunctor<GL>& functor,
                   EntityFilter<GL>*&& entity_filter = new ApplyOn::AllEntities<GL>(),
                   IntersectionFilter<GL>*&& intersection_filter = new ApplyOn::AllIntersections<GL>())
  {
    if (&functor == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    entity_and_intersection_functor_wrappers_.emplace_back(new internal::EntityAndIntersectionFunctorWrapper<GL>(
        functor, std::move(entity_filter), std::move(intersection_filter)));
    return *this;
  }

  ThisType& append(EntityAndIntersectionFunctor<GL>& functor,
                   IntersectionFilter<GL>*&& intersection_filter,
                   EntityFilter<GL>*&& entity_filter = new ApplyOn::AllEntities<GL>())
  {
    if (&functor == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    entity_and_intersection_functor_wrappers_.emplace_back(new internal::EntityAndIntersectionFunctorWrapper<GL>(
        functor, std::move(entity_filter), std::move(intersection_filter)));
    return *this;
  }

  ThisType& append(EntityAndIntersectionFunctor<GL>& functor,
                   std::function<bool(const GL&, const EntityType&)> entity_filter,
                   std::function<bool(const GL&, const IntersectionType&)> intersection_filter)
  {
    if (&functor == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    entity_and_intersection_functor_wrappers_.emplace_back(new internal::EntityAndIntersectionFunctorWrapper<GL>(
        functor,
        new ApplyOn::LambdaFilteredEntities<GL>(entity_filter),
        new ApplyOn::LambdaFilteredIntersections<GL>(intersection_filter)));
    return *this;
  }

  ThisType& append(EntityAndIntersectionFunctor<GL>*&& functor,
                   EntityFilter<GL>*&& entity_filter,
                   IntersectionFilter<GL>*&& intersection_filter)
  {
    if (functor == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    entity_and_intersection_functor_wrappers_.emplace_back(new internal::EntityAndIntersectionFunctorWrapper<GL>(
        std::move(functor), std::move(entity_filter), std::move(intersection_filter)));
    return *this;
  }

  /**
   * \}
   * \name These methods can be used to append entity and intersection lambda expressions.
   * \{
   */

  ThisType&
  append(std::function<void(const EntityType&)> entity_apply_on,
         std::function<void(const IntersectionType&, const EntityType&, const EntityType&)> intersection_apply_on,
         EntityFilter<GL>*&& entity_filter = new ApplyOn::AllEntities<GL>(),
         IntersectionFilter<GL>*&& intersection_filter = new ApplyOn::AllIntersections<GL>())
  {
    return this->append(new EntityAndIntersectionLambdaFunctor<GL>(entity_apply_on, intersection_apply_on),
                        std::move(entity_filter),
                        std::move(intersection_filter));
  }

  ThisType&
  append(std::function<void(const EntityType&)> entity_apply_on,
         std::function<void(const IntersectionType&, const EntityType&, const EntityType&)> intersection_apply_on,
         IntersectionFilter<GL>*&& intersection_filter,
         EntityFilter<GL>*&& entity_filter = new ApplyOn::AllEntities<GL>())
  {
    return this->append(new EntityAndIntersectionLambdaFunctor<GL>(entity_apply_on, intersection_apply_on),
                        std::move(entity_filter),
                        std::move(intersection_filter));
  }

  ThisType&
  append(std::function<void(const EntityType&)> entity_apply_on,
         std::function<void(const IntersectionType&, const EntityType&, const EntityType&)> intersection_apply_on,
         std::function<bool(const GL&, const EntityType&)> entity_filter,
         std::function<bool(const GL&, const IntersectionType&)> intersection_filter)
  {
    return this->append(new EntityAndIntersectionLambdaFunctor<GL>(entity_apply_on, intersection_apply_on),
                        new ApplyOn::LambdaFilteredEntities<GL>(entity_filter),
                        new ApplyOn::LambdaFilteredIntersections<GL>(intersection_filter));
  }

  ThisType&
  append(std::function<void(const EntityType&)> entity_apply_on,
         std::function<void(const IntersectionType&, const EntityType&, const EntityType&)> intersection_apply_on,
         std::function<void()> prepare_lambda,
         std::function<void()> finalize_lambda,
         EntityFilter<GL>*&& entity_filter = new ApplyOn::AllEntities<GL>(),
         IntersectionFilter<GL>*&& intersection_filter = new ApplyOn::AllIntersections<GL>())
  {
    return this->append(new EntityAndIntersectionLambdaFunctor<GL>(
                            entity_apply_on, intersection_apply_on, prepare_lambda, finalize_lambda),
                        std::move(entity_filter),
                        std::move(intersection_filter));
  }

  ThisType&
  append(std::function<void(const EntityType&)> entity_apply_on,
         std::function<void(const IntersectionType&, const EntityType&, const EntityType&)> intersection_apply_on,
         std::function<void()> prepare_lambda,
         std::function<void()> finalize_lambda,
         IntersectionFilter<GL>*&& intersection_filter,
         EntityFilter<GL>*&& entity_filter = new ApplyOn::AllEntities<GL>())
  {
    return this->append(new EntityAndIntersectionLambdaFunctor<GL>(
                            entity_apply_on, intersection_apply_on, prepare_lambda, finalize_lambda),
                        std::move(entity_filter),
                        std::move(intersection_filter));
  }

  ThisType&
  append(std::function<void(const EntityType&)> entity_apply_on,
         std::function<void(const IntersectionType&, const EntityType&, const EntityType&)> intersection_apply_on,
         std::function<void()> prepare_lambda,
         std::function<void()> finalize_lambda,
         std::function<bool(const GL&, const EntityType&)> entity_filter,
         std::function<bool(const GL&, const IntersectionType&)> intersection_filter)
  {
    return this->append(new EntityAndIntersectionLambdaFunctor<GL>(
                            entity_apply_on, intersection_apply_on, prepare_lambda, finalize_lambda),
                        new ApplyOn::LambdaFilteredEntities<GL>(entity_filter),
                        new ApplyOn::LambdaFilteredIntersections<GL>(intersection_filter));
  }

  /**
   * \}
   * \name These methods can be used to append another Walker.
   * \{
   */

  /**
   * \note The other_walker will be applied on the intersection of the given entity_filter (intersection_filter) and the
   *       filters of its EntityFunctors (IntersectionFunctors).
   * \sa   WalkerWrapper
   */
  ThisType& append(Walker<GL>& other_walker,
                   EntityFilter<GL>*&& entity_filter = new ApplyOn::AllEntities<GL>(),
                   IntersectionFilter<GL>*&& intersection_filter = new ApplyOn::AllIntersections<GL>())
  {
    if (&other_walker == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    entity_and_intersection_functor_wrappers_.emplace_back(
        new internal::WalkerWrapper<GL>(other_walker, std::move(entity_filter), std::move(intersection_filter)));
    return *this;
  }

  /**
   * \note The other_walker will be applied on the intersection of the given entity_filter (intersection_filter) and the
   *       filters of its EntityFunctors (IntersectionFunctors).
   * \sa   WalkerWrapper
   */
  ThisType& append(Walker<GL>& other_walker,
                   IntersectionFilter<GL>*&& intersection_filter,
                   EntityFilter<GL>*&& entity_filter = new ApplyOn::AllEntities<GL>())
  {
    if (&other_walker == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    entity_and_intersection_functor_wrappers_.emplace_back(
        new internal::WalkerWrapper<GL>(other_walker, std::move(entity_filter), std::move(intersection_filter)));
    return *this;
  }

  /**
   * \note The other_walker will be applied on the intersection of the given entity_filter (intersection_filter) and the
   *       filters of its EntityFunctors (IntersectionFunctors).
   * \sa   WalkerWrapper
   */
  ThisType& append(Walker<GL>& other_walker,
                   std::function<bool(const GL&, const EntityType&)> entity_filter,
                   std::function<bool(const GL&, const IntersectionType&)> intersection_filter)
  {
    if (&other_walker == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    entity_and_intersection_functor_wrappers_.emplace_back(
        new internal::WalkerWrapper<GL>(other_walker,
                                        new ApplyOn::LambdaFilteredEntities<GL>(entity_filter),
                                        new ApplyOn::LambdaFilteredIntersections<GL>(intersection_filter)));
    return *this;
  }

  /**
   * \}
   * \name These methods are required by EntityAndIntersectionFunctor.
   * \{
   */

  virtual void prepare() override
  {
    for (auto& wraper : entity_functor_wrappers_)
      wraper->functor().prepare();
    for (auto& wraper : intersection_functor_wrappers_)
      wraper->functor().prepare();
    for (auto& wraper : entity_and_intersection_functor_wrappers_)
      wraper->functor().prepare();
  } // ... prepare()

  virtual void apply_local(const EntityType& entity) override
  {
    for (auto& wraper : entity_functor_wrappers_) {
      if (wraper->filter().contains(grid_layer_, entity))
        wraper->functor().apply_local(entity);
    }
    for (auto& wraper : entity_and_intersection_functor_wrappers_) {
      if (wraper->entity_filter().contains(grid_layer_, entity))
        wraper->functor().apply_local(entity);
    }
  } // ... apply_local(...)

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& inside_entity,
                           const EntityType& outside_entity) override
  {
    for (auto& wraper : intersection_functor_wrappers_) {
      if (wraper->filter().contains(grid_layer_, intersection))
        wraper->functor().apply_local(intersection, inside_entity, outside_entity);
    }
    for (auto& wraper : entity_and_intersection_functor_wrappers_) {
      if (wraper->intersection_filter().contains(grid_layer_, intersection))
        wraper->functor().apply_local(intersection, inside_entity, outside_entity);
    }
  } // ... apply_local(...)

  virtual void finalize() override
  {
    for (auto& wraper : entity_functor_wrappers_)
      wraper->functor().finalize();
    for (auto& wraper : intersection_functor_wrappers_)
      wraper->functor().finalize();
    for (auto& wraper : entity_and_intersection_functor_wrappers_)
      wraper->functor().finalize();
  } // ... finalize()

  /**
   * \}
   */

  void walk(const bool use_tbb = false, const bool clear_functors = true)
  {
#if HAVE_TBB
    if (use_tbb) {
      const auto num_partitions =
          DXTC_CONFIG_GET("threading.partition_factor", 1u) * XT::Common::threadManager().current_threads();
      RangedPartitioning<GridLayerType, 0> partitioning(grid_layer_, num_partitions);
      this->walk(partitioning, clear_functors);
      return;
    }
#else
    DUNE_UNUSED const auto no_warning_for_use_tbb = use_tbb;
#endif
    // prepare functors
    prepare();

    // only do something, if we have to
    if ((entity_functor_wrappers_.size() + intersection_functor_wrappers_.size()
         + entity_and_intersection_functor_wrappers_.size())
        > 0) {
      walk_range(elements(grid_layer_));
    } // only do something, if we have to

    // finalize functors
    finalize();

    if (clear_functors)
      clear();
  } // ... walk(...)

  void clear()
  {
    entity_functor_wrappers_.clear();
    intersection_functor_wrappers_.clear();
    entity_and_intersection_functor_wrappers_.clear();
  }

#if HAVE_TBB
protected:
  template <class PartioningType, class WalkerType>
  struct Body
  {
    Body(WalkerType& walker, PartioningType& partitioning)
      : walker_(walker)
      , partitioning_(partitioning)
    {
    }

    Body(Body& other, tbb::split /*split*/)
      : walker_(other.walker_)
      , partitioning_(other.partitioning_)
    {
    }

    void operator()(const tbb::blocked_range<std::size_t>& range) const
    {
      // for all partitions in tbb-range
      for (std::size_t p = range.begin(); p != range.end(); ++p) {
        auto partition = partitioning_.partition(p);
        walker_.walk_range(partition);
      }
    }

    void join(Body& /*other*/)
    {
    }

    WalkerType& walker_;
    const PartioningType& partitioning_;
  }; // struct Body

public:
  template <class PartioningType>
  void walk(PartioningType& partitioning, const bool clear_functors = true)
  {
    // prepare functors
    prepare();

    // only do something, if we have to
    if ((entity_functor_wrappers_.size() + intersection_functor_wrappers_.size()
         + entity_and_intersection_functor_wrappers_.size())
        > 0) {
      tbb::blocked_range<std::size_t> range(0, partitioning.partitions());
      Body<PartioningType, ThisType> body(*this, partitioning);
      tbb::parallel_reduce(range, body);
    }

    // finalize functors
    finalize();

    if (clear_functors)
      clear();
  } // ... tbb_walk(...)
#else
public:
  template <class PartioningType>
  void walk(PartioningType& partitioning)
  {
    // prepare functors
    prepare();

    // only do something, if we have to
    if ((codim0_functors_.size() + codim1_functors_.size()) > 0) {
      // no actual SMP walk, use range as is
      walk_range(partitioning.everything());
    }

    // finalize functors
    finalize();
    clear();
  } // ... tbb_walk(...)
#endif // HAVE_TBB

protected:
  template <class EntityRange>
  void walk_range(const EntityRange& entity_range)
  {
#ifdef __INTEL_COMPILER
    const auto it_end = entity_range.end();
    for (auto it = entity_range.begin(); it != it_end; ++it) {
      const EntityType& entity = *it;
#else
    for (const EntityType& entity : entity_range) {
#endif
      // apply entity functors
      apply_local(entity);

      // only walk the intersections, if there are codim1 functors present
      if ((intersection_functor_wrappers_.size() + entity_and_intersection_functor_wrappers_.size()) > 0) {
        // Do not use intersections(...) here, since that does not work for a SubdomainGridPart which is based on
        // alugrid and then wrapped as a grid view (see also https://github.com/dune-community/dune-xt-grid/issues/26)
        const auto intersection_it_end = grid_layer_.iend(entity);
        for (auto intersection_it = grid_layer_.ibegin(entity); intersection_it != intersection_it_end;
             ++intersection_it) {
          const auto& intersection = *intersection_it;
          if (intersection.neighbor()) {
            const auto neighbor = intersection.outside();
            apply_local(intersection, entity, neighbor);
          } else
            apply_local(intersection, entity, entity);
        } // walk the intersections
      } // only walk the intersections, if there are codim1 functors present
    } // .. walk entities
  } // ... walk_range(...)

  friend class internal::WalkerWrapper<GridLayerType>;

  GridLayerType grid_layer_;
  std::list<std::shared_ptr<internal::EntityFunctorWrapper<GridLayerType>>> entity_functor_wrappers_;
  std::list<std::shared_ptr<internal::IntersectionFunctorWrapper<GridLayerType>>> intersection_functor_wrappers_;
  std::list<std::shared_ptr<internal::EntityAndIntersectionFunctorWrapper<GridLayerType>>>
      entity_and_intersection_functor_wrappers_;
}; // class Walker


template <class GL>
Walker<GL> make_walker(GL grid_layer)
{
  return Walker<GL>(grid_layer);
}


} // namespace Grid
} // namespace XT
} // namespace Dune


#include "walker.lib.hh"


#endif // DUNE_XT_GRID_WALKER_HH
