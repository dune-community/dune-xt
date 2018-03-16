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
#include <dune/common/deprecated.hh>

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
#include <dune/xt/grid/type_traits.hh>

#include <dune/xt/grid/walker/apply-on.hh>
#include <dune/xt/grid/walker/functors.hh>
#include <dune/xt/grid/walker/wrapper.hh>

namespace Dune {
namespace XT {
namespace Grid {


template <class GridLayerImp>
class Walker : public Functor::Codim0And1<GridLayerImp>
{
  static_assert(is_layer<GridLayerImp>::value, "");
  typedef Walker<GridLayerImp> ThisType;

public:
  typedef GridLayerImp GridLayerType;
  using EntityType = extract_entity_t<GridLayerType>;
  using IntersectionType = extract_intersection_t<GridLayerType>;

  explicit Walker(GridLayerType grd_lr)
    : grid_layer_(grd_lr)
  {
  }

  /// \sa https://github.com/dune-community/dune-gdt/issues/89
  Walker(const Walker& other) = delete; // <- b.c. of the functors: type-erasue = no copy!

  Walker(Walker&& source) = default;

  const GridLayerType& grid_layer() const
  {
    return grid_layer_;
  }

  GridLayerType& grid_layer()
  {
    return grid_layer_;
  }

  ThisType& append(std::function<void(const EntityType&)> lambda,
                   const ApplyOn::WhichEntity<GridLayerType>* where = new ApplyOn::AllEntities<GridLayerType>())
  {
    codim0_functors_.emplace_back(new internal::Codim0LambdaWrapper<GridLayerType>(lambda, where));
    return *this;
  }

  ThisType&
  append(std::function<void(const IntersectionType&, const EntityType&, const EntityType&)> lambda,
         const ApplyOn::WhichIntersection<GridLayerType>* where = new ApplyOn::AllIntersections<GridLayerType>())
  {
    codim1_functors_.emplace_back(new internal::Codim1LambdaWrapper<GridLayerType>(lambda, where));
    return *this;
  }

  ThisType& append(Functor::Codim0<GridLayerType>& functor,
                   const ApplyOn::WhichEntity<GridLayerType>* where = new ApplyOn::AllEntities<GridLayerType>())
  {
    codim0_functors_.emplace_back(
        new internal::Codim0FunctorWrapper<GridLayerType, Functor::Codim0<GridLayerType>>(functor, where));
    return *this;
  }

  ThisType&
  append(Functor::Codim1<GridLayerType>& functor,
         const ApplyOn::WhichIntersection<GridLayerType>* where = new ApplyOn::AllIntersections<GridLayerType>())
  {
    codim1_functors_.emplace_back(
        new internal::Codim1FunctorWrapper<GridLayerType, Functor::Codim1<GridLayerType>>(functor, where));
    return *this;
  }

  ThisType&
  append(Functor::Codim0And1<GridLayerType>& functor,
         const ApplyOn::WhichEntity<GridLayerType>* which_entities = new ApplyOn::AllEntities<GridLayerType>(),
         const ApplyOn::WhichIntersection<GridLayerType>* which_intersections =
             new ApplyOn::AllIntersections<GridLayerType>())
  {
    codim0_functors_.emplace_back(
        new internal::Codim0FunctorWrapper<GridLayerType, Functor::Codim0And1<GridLayerType>>(functor, which_entities));
    codim1_functors_.emplace_back(new internal::Codim1FunctorWrapper<GridLayerType, Functor::Codim0And1<GridLayerType>>(
        functor, which_intersections));
    return *this;
  }

  ThisType&
  append(Functor::Codim0And1<GridLayerType>& functor,
         const ApplyOn::WhichIntersection<GridLayerType>* which_intersections,
         const ApplyOn::WhichEntity<GridLayerType>* which_entities = new ApplyOn::AllEntities<GridLayerType>())
  {
    codim0_functors_.emplace_back(
        new internal::Codim0FunctorWrapper<GridLayerType, Functor::Codim0And1<GridLayerType>>(functor, which_entities));
    codim1_functors_.emplace_back(new internal::Codim1FunctorWrapper<GridLayerType, Functor::Codim0And1<GridLayerType>>(
        functor, which_intersections));
    return *this;
  }

  ThisType&
  append(ThisType& other,
         const ApplyOn::WhichEntity<GridLayerType>* which_entities = new ApplyOn::AllEntities<GridLayerType>(),
         const ApplyOn::WhichIntersection<GridLayerType>* which_intersections =
             new ApplyOn::AllIntersections<GridLayerType>())
  {
    if (&other == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    codim0_functors_.emplace_back(new internal::WalkerWrapper<GridLayerType, ThisType>(other, which_entities));
    codim1_functors_.emplace_back(new internal::WalkerWrapper<GridLayerType, ThisType>(other, which_intersections));
    return *this;
  } // ... append(...)

  ThisType&
  append(ThisType& other,
         const ApplyOn::WhichIntersection<GridLayerType>* which_intersections,
         const ApplyOn::WhichEntity<GridLayerType>* which_entities = new ApplyOn::AllEntities<GridLayerType>())
  {
    if (&other == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    codim0_functors_.emplace_back(new internal::WalkerWrapper<GridLayerType, ThisType>(other, which_entities));
    codim1_functors_.emplace_back(new internal::WalkerWrapper<GridLayerType, ThisType>(other, which_intersections));
    return *this;
  } // ... append(...)

  void clear()
  {
    codim0_functors_.clear();
    codim1_functors_.clear();
  } // ... clear()

  virtual void prepare()
  {
    for (auto& functor : codim0_functors_)
      functor->prepare();
    for (auto& functor : codim1_functors_)
      functor->prepare();
  } // ... prepare()

  bool apply_on(const EntityType& entity) const
  {
    for (const auto& functor : codim0_functors_)
      if (functor->apply_on(grid_layer_, entity))
        return true;
    return false;
  } // ... apply_on(...)

  bool apply_on(const IntersectionType& intersection) const
  {
    for (const auto& functor : codim1_functors_)
      if (functor->apply_on(grid_layer_, intersection))
        return true;
    return false;
  } // ... apply_on(...)

  virtual void apply_local(const EntityType& entity)
  {
    for (auto& functor : codim0_functors_)
      if (functor->apply_on(grid_layer_, entity))
        functor->apply_local(entity);
  } // ... apply_local(...)

  virtual void
  apply_local(const IntersectionType& intersection, const EntityType& inside_entity, const EntityType& outside_entity)
  {
    for (auto& functor : codim1_functors_)
      if (functor->apply_on(grid_layer_, intersection))
        functor->apply_local(intersection, inside_entity, outside_entity);
  } // ... apply_local(...)

  virtual void finalize()
  {
    for (auto& functor : codim0_functors_)
      functor->finalize();
    for (auto& functor : codim1_functors_)
      functor->finalize();
  } // ... finalize()

  void walk(const bool use_tbb = false)
  {
#if HAVE_TBB
    if (use_tbb) {
      const auto num_partitions =
          DXTC_CONFIG_GET("threading.partition_factor", 1u) * XT::Common::threadManager().current_threads();
      RangedPartitioning<GridLayerType, 0> partitioning(grid_layer_, num_partitions);
      this->walk(partitioning);
      return;
    }
#else
    DUNE_UNUSED const auto no_warning_for_use_tbb = use_tbb;
#endif
    // prepare functors
    prepare();

    // only do something, if we have to
    if ((codim0_functors_.size() + codim1_functors_.size()) > 0) {
      walk_range(elements(grid_layer_));
    } // only do something, if we have to

    // finalize functors
    finalize();
    clear();
  } // ... walk(...)

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
  void walk(PartioningType& partitioning)
  {
    // prepare functors
    prepare();

    // only do something, if we have to
    if ((codim0_functors_.size() + codim1_functors_.size()) > 0) {
      tbb::blocked_range<std::size_t> range(0, partitioning.partitions());
      Body<PartioningType, ThisType> body(*this, partitioning);
      tbb::parallel_reduce(range, body);
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
      // apply codim0 functors
      apply_local(entity);

      // only walk the intersections, if there are codim1 functors present
      if (codim1_functors_.size() > 0) {
        // walk the intersections, do not use intersections(...) here, since that does not work for a SubdomainGridView
        // which is based on alugrid and then wrapped as a grid view (see also
        // https://github.com/dune-community/dune-xt-grid/issues/26)
        const auto intersection_it_end = grid_layer_.iend(entity);
        for (auto intersection_it = grid_layer_.ibegin(entity); intersection_it != intersection_it_end;
             ++intersection_it) {
          const auto& intersection = *intersection_it;

          // apply codim1 functors
          if (intersection.neighbor()) {
            const auto neighbor = intersection.outside();
            apply_local(intersection, entity, neighbor);
          } else
            apply_local(intersection, entity, entity);

        } // walk the intersections
      } // only walk the intersections, if there are codim1 functors present
    }
  } // ... walk_range(...)

  GridLayerType grid_layer_;
  std::vector<std::unique_ptr<internal::Codim0Object<GridLayerType>>> codim0_functors_;
  std::vector<std::unique_ptr<internal::Codim1Object<GridLayerType>>> codim1_functors_;
}; // class Walker


} // namespace Grid
} // namespace XT
} // namespace Dune


#include "walker.lib.hh"


#endif // DUNE_XT_GRID_WALKER_HH
