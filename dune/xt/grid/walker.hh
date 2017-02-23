// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014 - 2016)
//   Tobias Leibner  (2014 - 2015, 2017)

#ifndef DUNE_XT_GRID_WALKER_HH
#define DUNE_XT_GRID_WALKER_HH

#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

#include <dune/common/version.hh>
#include <dune/common/deprecated.hh>

#if HAVE_TBB
#include <dune/xt/grid/parallel/partitioning/ranged.hh>
#include <dune/xt/common/parallel/threadmanager.hh>

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#include <tbb/tbb_stddef.h>
#endif // HAVE_TBB

#include <dune/xt/common/parallel/threadmanager.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/unused.hh>
#include <dune/xt/common/timedlogging.hh>

#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/rangegenerators.hh>

#include "walker/apply-on.hh"
#include "walker/functors.hh"
#include "walker/wrapper.hh"

namespace Dune {
namespace XT {
namespace Grid {
namespace internal {


template <class GPV, bool is_grd_vw = is_grid_view<GPV>::value>
struct GridPartViewHolder
{
  typedef GPV type;

  GridPartViewHolder(GPV grid_view)
    : grid_view_(grid_view)
  {
  }

  GridPartViewHolder(const GridPartViewHolder& other) = default;
  GridPartViewHolder(GridPartViewHolder&& source) = default;

  const type& real_grid_view() const
  {
    return grid_view_;
  }

  const GPV grid_view_;
};

template <class GPV>
struct GridPartViewHolder<GPV, false>
{
  typedef typename GPV::GridViewType type;

  GridPartViewHolder(const GridPartViewHolder& other) = default;
  GridPartViewHolder(GridPartViewHolder&& source) = default;

  GridPartViewHolder(GPV grid_part)
    : grid_view_(grid_part)
    , grid_part_(grid_view_.gridView())
  {
  }

  const type& real_grid_view() const
  {
    return grid_part_;
  }

  const GPV grid_view_;
  const type grid_part_;
};


} // namespace internal


template <class GridViewImp>
class Walker : internal::GridPartViewHolder<GridViewImp>, public Functor::Codim0And1<GridViewImp>
{
  typedef Walker<GridViewImp> ThisType;

public:
  typedef GridViewImp GridViewType;
  typedef typename XT::Grid::Entity<GridViewType>::Type EntityType;
  typedef typename XT::Grid::Intersection<GridViewType>::Type IntersectionType;

  explicit Walker(GridViewType grd_vw)
    : internal::GridPartViewHolder<GridViewImp>(grd_vw)
  {
  }

  /// \sa https://github.com/dune-community/dune-gdt/issues/89
  Walker(const Walker& other) = delete; // <- b.c. of the functors: type-erasue = no copy!

  Walker(Walker&& source) = default;

  const GridViewType& grid_view() const
  {
    return this->grid_view_;
  }

  ThisType& append(std::function<void(const EntityType&)> lambda,
                   const ApplyOn::WhichEntity<GridViewType>* where = new ApplyOn::AllEntities<GridViewType>())
  {
    codim0_functors_.emplace_back(new internal::Codim0LambdaWrapper<GridViewType>(lambda, where));
    return *this;
  }

  ThisType&
  append(std::function<void(const IntersectionType&, const EntityType&, const EntityType&)> lambda,
         const ApplyOn::WhichIntersection<GridViewType>* where = new ApplyOn::AllIntersections<GridViewType>())
  {
    codim1_functors_.emplace_back(new internal::Codim1LambdaWrapper<GridViewType>(lambda, where));
    return *this;
  }

  ThisType& append(Functor::Codim0<GridViewType>& functor,
                   const ApplyOn::WhichEntity<GridViewType>* where = new ApplyOn::AllEntities<GridViewType>())
  {
    codim0_functors_.emplace_back(
        new internal::Codim0FunctorWrapper<GridViewType, Functor::Codim0<GridViewType>>(functor, where));
    return *this;
  }

  ThisType&
  append(Functor::Codim1<GridViewType>& functor,
         const ApplyOn::WhichIntersection<GridViewType>* where = new ApplyOn::AllIntersections<GridViewType>())
  {
    codim1_functors_.emplace_back(
        new internal::Codim1FunctorWrapper<GridViewType, Functor::Codim1<GridViewType>>(functor, where));
    return *this;
  }

  ThisType& append(Functor::Codim0And1<GridViewType>& functor,
                   const ApplyOn::WhichEntity<GridViewType>* which_entities = new ApplyOn::AllEntities<GridViewType>(),
                   const ApplyOn::WhichIntersection<GridViewType>* which_intersections =
                       new ApplyOn::AllIntersections<GridViewType>())
  {
    codim0_functors_.emplace_back(
        new internal::Codim0FunctorWrapper<GridViewType, Functor::Codim0And1<GridViewType>>(functor, which_entities));
    codim1_functors_.emplace_back(new internal::Codim1FunctorWrapper<GridViewType, Functor::Codim0And1<GridViewType>>(
        functor, which_intersections));
    return *this;
  }

  ThisType& append(Functor::Codim0And1<GridViewType>& functor,
                   const ApplyOn::WhichIntersection<GridViewType>* which_intersections,
                   const ApplyOn::WhichEntity<GridViewType>* which_entities = new ApplyOn::AllEntities<GridViewType>())
  {
    codim0_functors_.emplace_back(
        new internal::Codim0FunctorWrapper<GridViewType, Functor::Codim0And1<GridViewType>>(functor, which_entities));
    codim1_functors_.emplace_back(new internal::Codim1FunctorWrapper<GridViewType, Functor::Codim0And1<GridViewType>>(
        functor, which_intersections));
    return *this;
  }

  ThisType& append(ThisType& other,
                   const ApplyOn::WhichEntity<GridViewType>* which_entities = new ApplyOn::AllEntities<GridViewType>(),
                   const ApplyOn::WhichIntersection<GridViewType>* which_intersections =
                       new ApplyOn::AllIntersections<GridViewType>())
  {
    if (&other == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    codim0_functors_.emplace_back(new internal::WalkerWrapper<GridViewType, ThisType>(other, which_entities));
    codim1_functors_.emplace_back(new internal::WalkerWrapper<GridViewType, ThisType>(other, which_intersections));
    return *this;
  } // ... append(...)

  ThisType& append(ThisType& other,
                   const ApplyOn::WhichIntersection<GridViewType>* which_intersections,
                   const ApplyOn::WhichEntity<GridViewType>* which_entities = new ApplyOn::AllEntities<GridViewType>())
  {
    if (&other == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    codim0_functors_.emplace_back(new internal::WalkerWrapper<GridViewType, ThisType>(other, which_entities));
    codim1_functors_.emplace_back(new internal::WalkerWrapper<GridViewType, ThisType>(other, which_intersections));
    return *this;
  } // ... append(...)

  template <class... Args>
  DUNE_DEPRECATED_MSG("Use append() instead (since 11.01.2017)!")
  ThisType& add(Args&&... args)
  {
    return append(std::forward<Args>(args)...);
  }

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
      if (functor->apply_on(this->grid_view_, entity))
        return true;
    return false;
  } // ... apply_on(...)

  bool apply_on(const IntersectionType& intersection) const
  {
    for (const auto& functor : codim1_functors_)
      if (functor->apply_on(this->grid_view_, intersection))
        return true;
    return false;
  } // ... apply_on(...)

  virtual void apply_local(const EntityType& entity)
  {
    for (auto& functor : codim0_functors_)
      if (functor->apply_on(this->grid_view_, entity))
        functor->apply_local(entity);
  } // ... apply_local(...)

  virtual void
  apply_local(const IntersectionType& intersection, const EntityType& inside_entity, const EntityType& outside_entity)
  {
    for (auto& functor : codim1_functors_)
      if (functor->apply_on(this->grid_view_, intersection))
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
      RangedPartitioning<typename internal::GridPartViewHolder<GridViewType>::type, 0> partitioning(
          this->real_grid_view(), num_partitions);
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
      walk_range(elements(this->grid_view_));
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
        // walk the intersections
        for (auto&& intersection : intersections(this->grid_view_, entity)) {

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

  std::vector<std::unique_ptr<internal::Codim0Object<GridViewType>>> codim0_functors_;
  std::vector<std::unique_ptr<internal::Codim1Object<GridViewType>>> codim1_functors_;
}; // class Walker


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_WALKER_HH
