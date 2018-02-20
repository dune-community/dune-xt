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
class Walker : public ElementAndIntersectionFunctor<GL>
{
  using BaseType = ElementAndIntersectionFunctor<GL>;
  using ThisType = Walker<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;
  using typename BaseType::IntersectionType;

  explicit Walker(GridViewType grd_vw)
    : grid_view_(grd_vw)
  {
  }

  Walker(const Walker& other) = delete;
  Walker(Walker&& source) = default;

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  GridViewType& grid_view()
  {
    return grid_view_;
  }

  /**
   * \name These methods can be used to append an \sa ElementFunctor.
   * \{
   */

  ThisType& append(ElementFunctor<GL>& functor, ElementFilter<GL>*&& filter = new ApplyOn::AllElements<GL>())
  {
    element_functor_wrappers_.emplace_back(new internal::ElementFunctorWrapper<GL>(functor, std::move(filter)));
    return *this;
  }

  ThisType& append(ElementFunctor<GL>& functor, std::function<bool(const GL&, const ElementType&)> element_filter)
  {
    element_functor_wrappers_.emplace_back(
        new internal::ElementFunctorWrapper<GL>(functor, new ApplyOn::LambdaFilteredElements<GL>(element_filter)));
    return *this;
  }

  ThisType& append(ElementFunctor<GL>*&& functor, ElementFilter<GL>*&& filter)
  {
    element_functor_wrappers_.emplace_back(
        new internal::ElementFunctorWrapper<GL>(std::move(functor), std::move(filter)));
    return *this;
  }

  /**
   * \}
   * \name These methods can be used to append an element lambda expression.
   * \{
   */

  ThisType& append(std::function<void(const ElementType&)> apply_lambda,
                   ElementFilter<GL>*&& filter = new ApplyOn::AllElements<GL>())
  {
    return this->append(new ElementLambdaFunctor<GL>(apply_lambda), std::move(filter));
  }

  ThisType& append(std::function<void(const ElementType&)> apply_lambda,
                   std::function<bool(const GL&, const ElementType&)> filter)
  {
    return this->append(new ElementLambdaFunctor<GL>(apply_lambda), new ApplyOn::LambdaFilteredElements<GL>(filter));
  }

  ThisType& append(std::function<void(const ElementType&)> apply_lambda,
                   std::function<void()> prepare_lambda,
                   std::function<void()> finalize_lambda,
                   ElementFilter<GL>*&& filter = new ApplyOn::AllElements<GL>())
  {
    return this->append(new ElementLambdaFunctor<GL>(apply_lambda, prepare_lambda, finalize_lambda), std::move(filter));
  }

  ThisType& append(std::function<void(const ElementType&)> apply_lambda,
                   std::function<void()> prepare_lambda,
                   std::function<void()> finalize_lambda,
                   std::function<bool(const GL&, const ElementType&)> filter)
  {
    return this->append(new ElementLambdaFunctor<GL>(apply_lambda, prepare_lambda, finalize_lambda),
                        new ApplyOn::LambdaFilteredElements<GL>(filter));
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

  ThisType& append(std::function<void(const IntersectionType&, const ElementType&, const ElementType&)> apply_lambda,
                   IntersectionFilter<GL>*&& filter = new ApplyOn::AllIntersections<GL>())
  {
    return this->append(new IntersectionLambdaFunctor<GL>(apply_lambda), std::move(filter));
  }

  ThisType& append(std::function<void(const IntersectionType&, const ElementType&, const ElementType&)> apply_lambda,
                   std::function<bool(const GL&, const IntersectionType&)> filter)
  {
    return this->append(new IntersectionLambdaFunctor<GL>(apply_lambda),
                        new ApplyOn::LambdaFilteredIntersections<GL>(filter));
  }

  ThisType& append(std::function<void(const IntersectionType&, const ElementType&, const ElementType&)> apply_lambda,
                   std::function<void()> prepare_lambda,
                   std::function<void()> finalize_lambda,
                   IntersectionFilter<GL>*&& filter = new ApplyOn::AllIntersections<GL>())
  {
    return this->append(new IntersectionLambdaFunctor<GL>(apply_lambda, prepare_lambda, finalize_lambda),
                        std::move(filter));
  }

  ThisType& append(std::function<void(const IntersectionType&, const ElementType&, const ElementType&)> apply_lambda,
                   std::function<void()> prepare_lambda,
                   std::function<void()> finalize_lambda,
                   std::function<bool(const GL&, const IntersectionType&)> filter)
  {
    return this->append(new IntersectionLambdaFunctor<GL>(apply_lambda, prepare_lambda, finalize_lambda),
                        new ApplyOn::LambdaFilteredIntersections<GL>(filter));
  }

  /**
   * \}
   * \name These methods can be used to append an \sa ElementAndIntersectionFunctor.
   * \{
   */

  ThisType& append(ElementAndIntersectionFunctor<GL>& functor,
                   ElementFilter<GL>*&& element_filter = new ApplyOn::AllElements<GL>(),
                   IntersectionFilter<GL>*&& intersection_filter = new ApplyOn::AllIntersections<GL>())
  {
    if (&functor == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    element_and_intersection_functor_wrappers_.emplace_back(new internal::ElementAndIntersectionFunctorWrapper<GL>(
        functor, std::move(element_filter), std::move(intersection_filter)));
    return *this;
  }

  ThisType& append(ElementAndIntersectionFunctor<GL>& functor,
                   IntersectionFilter<GL>*&& intersection_filter,
                   ElementFilter<GL>*&& element_filter = new ApplyOn::AllElements<GL>())
  {
    if (&functor == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    element_and_intersection_functor_wrappers_.emplace_back(new internal::ElementAndIntersectionFunctorWrapper<GL>(
        functor, std::move(element_filter), std::move(intersection_filter)));
    return *this;
  }

  ThisType& append(ElementAndIntersectionFunctor<GL>& functor,
                   std::function<bool(const GL&, const ElementType&)> element_filter,
                   std::function<bool(const GL&, const IntersectionType&)> intersection_filter)
  {
    if (&functor == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    element_and_intersection_functor_wrappers_.emplace_back(new internal::ElementAndIntersectionFunctorWrapper<GL>(
        functor,
        new ApplyOn::LambdaFilteredElements<GL>(element_filter),
        new ApplyOn::LambdaFilteredIntersections<GL>(intersection_filter)));
    return *this;
  }

  ThisType& append(ElementAndIntersectionFunctor<GL>*&& functor,
                   ElementFilter<GL>*&& element_filter,
                   IntersectionFilter<GL>*&& intersection_filter)
  {
    if (functor == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    element_and_intersection_functor_wrappers_.emplace_back(new internal::ElementAndIntersectionFunctorWrapper<GL>(
        std::move(functor), std::move(element_filter), std::move(intersection_filter)));
    return *this;
  }

  /**
   * \}
   * \name These methods can be used to append element and intersection lambda expressions.
   * \{
   */

  ThisType&
  append(std::function<void(const ElementType&)> element_apply_on,
         std::function<void(const IntersectionType&, const ElementType&, const ElementType&)> intersection_apply_on,
         ElementFilter<GL>*&& element_filter = new ApplyOn::AllElements<GL>(),
         IntersectionFilter<GL>*&& intersection_filter = new ApplyOn::AllIntersections<GL>())
  {
    return this->append(new ElementAndIntersectionLambdaFunctor<GL>(element_apply_on, intersection_apply_on),
                        std::move(element_filter),
                        std::move(intersection_filter));
  }

  ThisType&
  append(std::function<void(const ElementType&)> element_apply_on,
         std::function<void(const IntersectionType&, const ElementType&, const ElementType&)> intersection_apply_on,
         IntersectionFilter<GL>*&& intersection_filter,
         ElementFilter<GL>*&& element_filter = new ApplyOn::AllElements<GL>())
  {
    return this->append(new ElementAndIntersectionLambdaFunctor<GL>(element_apply_on, intersection_apply_on),
                        std::move(element_filter),
                        std::move(intersection_filter));
  }

  ThisType&
  append(std::function<void(const ElementType&)> element_apply_on,
         std::function<void(const IntersectionType&, const ElementType&, const ElementType&)> intersection_apply_on,
         std::function<bool(const GL&, const ElementType&)> element_filter,
         std::function<bool(const GL&, const IntersectionType&)> intersection_filter)
  {
    return this->append(new ElementAndIntersectionLambdaFunctor<GL>(element_apply_on, intersection_apply_on),
                        new ApplyOn::LambdaFilteredElements<GL>(element_filter),
                        new ApplyOn::LambdaFilteredIntersections<GL>(intersection_filter));
  }

  ThisType&
  append(std::function<void(const ElementType&)> element_apply_on,
         std::function<void(const IntersectionType&, const ElementType&, const ElementType&)> intersection_apply_on,
         std::function<void()> prepare_lambda,
         std::function<void()> finalize_lambda,
         ElementFilter<GL>*&& element_filter = new ApplyOn::AllElements<GL>(),
         IntersectionFilter<GL>*&& intersection_filter = new ApplyOn::AllIntersections<GL>())
  {
    return this->append(new ElementAndIntersectionLambdaFunctor<GL>(
                            element_apply_on, intersection_apply_on, prepare_lambda, finalize_lambda),
                        std::move(element_filter),
                        std::move(intersection_filter));
  }

  ThisType&
  append(std::function<void(const ElementType&)> element_apply_on,
         std::function<void(const IntersectionType&, const ElementType&, const ElementType&)> intersection_apply_on,
         std::function<void()> prepare_lambda,
         std::function<void()> finalize_lambda,
         IntersectionFilter<GL>*&& intersection_filter,
         ElementFilter<GL>*&& element_filter = new ApplyOn::AllElements<GL>())
  {
    return this->append(new ElementAndIntersectionLambdaFunctor<GL>(
                            element_apply_on, intersection_apply_on, prepare_lambda, finalize_lambda),
                        std::move(element_filter),
                        std::move(intersection_filter));
  }

  ThisType&
  append(std::function<void(const ElementType&)> element_apply_on,
         std::function<void(const IntersectionType&, const ElementType&, const ElementType&)> intersection_apply_on,
         std::function<void()> prepare_lambda,
         std::function<void()> finalize_lambda,
         std::function<bool(const GL&, const ElementType&)> element_filter,
         std::function<bool(const GL&, const IntersectionType&)> intersection_filter)
  {
    return this->append(new ElementAndIntersectionLambdaFunctor<GL>(
                            element_apply_on, intersection_apply_on, prepare_lambda, finalize_lambda),
                        new ApplyOn::LambdaFilteredElements<GL>(element_filter),
                        new ApplyOn::LambdaFilteredIntersections<GL>(intersection_filter));
  }

  /**
   * \}
   * \name These methods can be used to append another Walker.
   * \{
   */

  /**
   * \note The other_walker will be applied on the intersection of the given element_filter (intersection_filter) and
   *       the filters of its ElementFunctors (IntersectionFunctors).
   * \sa   WalkerWrapper
   */
  ThisType& append(Walker<GL>& other_walker,
                   ElementFilter<GL>*&& element_filter = new ApplyOn::AllElements<GL>(),
                   IntersectionFilter<GL>*&& intersection_filter = new ApplyOn::AllIntersections<GL>())
  {
    if (&other_walker == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    element_and_intersection_functor_wrappers_.emplace_back(
        new internal::WalkerWrapper<GL>(other_walker, std::move(element_filter), std::move(intersection_filter)));
    return *this;
  }

  /**
   * \note The other_walker will be applied on the intersection of the given element_filter (intersection_filter) and
   *       the filters of its ElementFunctors (IntersectionFunctors).
   * \sa   WalkerWrapper
   */
  ThisType& append(Walker<GL>& other_walker,
                   IntersectionFilter<GL>*&& intersection_filter,
                   ElementFilter<GL>*&& element_filter = new ApplyOn::AllElements<GL>())
  {
    if (&other_walker == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    element_and_intersection_functor_wrappers_.emplace_back(
        new internal::WalkerWrapper<GL>(other_walker, std::move(element_filter), std::move(intersection_filter)));
    return *this;
  }

  /**
   * \note The other_walker will be applied on the intersection of the given element_filter (intersection_filter) and
   *       the filters of its ElementFunctors (IntersectionFunctors).
   * \sa   WalkerWrapper
   */
  ThisType& append(Walker<GL>& other_walker,
                   std::function<bool(const GL&, const ElementType&)> element_filter,
                   std::function<bool(const GL&, const IntersectionType&)> intersection_filter)
  {
    if (&other_walker == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    element_and_intersection_functor_wrappers_.emplace_back(
        new internal::WalkerWrapper<GL>(other_walker,
                                        new ApplyOn::LambdaFilteredElements<GL>(element_filter),
                                        new ApplyOn::LambdaFilteredIntersections<GL>(intersection_filter)));
    return *this;
  }

  /**
   * \}
   * \name These methods are required by ElementAndIntersectionFunctor.
   * \{
   */

  virtual void prepare() override
  {
    for (auto& wraper : element_functor_wrappers_)
      wraper->functor().prepare();
    for (auto& wraper : intersection_functor_wrappers_)
      wraper->functor().prepare();
    for (auto& wraper : element_and_intersection_functor_wrappers_)
      wraper->functor().prepare();
  } // ... prepare()

  virtual void apply_local(const ElementType& element) override
  {
    for (auto& wraper : element_functor_wrappers_) {
      if (wraper->filter().contains(grid_view_, element))
        wraper->functor().apply_local(element);
    }
    for (auto& wraper : element_and_intersection_functor_wrappers_) {
      if (wraper->element_filter().contains(grid_view_, element))
        wraper->functor().apply_local(element);
    }
  } // ... apply_local(...)

  virtual void apply_local(const IntersectionType& intersection,
                           const ElementType& inside_element,
                           const ElementType& outside_element) override
  {
    for (auto& wraper : intersection_functor_wrappers_) {
      if (wraper->filter().contains(grid_view_, intersection))
        wraper->functor().apply_local(intersection, inside_element, outside_element);
    }
    for (auto& wraper : element_and_intersection_functor_wrappers_) {
      if (wraper->intersection_filter().contains(grid_view_, intersection))
        wraper->functor().apply_local(intersection, inside_element, outside_element);
    }
  } // ... apply_local(...)

  virtual void finalize() override
  {
    for (auto& wraper : element_functor_wrappers_)
      wraper->functor().finalize();
    for (auto& wraper : intersection_functor_wrappers_)
      wraper->functor().finalize();
    for (auto& wraper : element_and_intersection_functor_wrappers_)
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
      RangedPartitioning<GridViewType, 0> partitioning(grid_view_, num_partitions);
      this->walk(partitioning, clear_functors);
      return;
    }
#else
    DUNE_UNUSED const auto no_warning_for_use_tbb = use_tbb;
#endif
    // prepare functors
    prepare();

    // only do something, if we have to
    if ((element_functor_wrappers_.size() + intersection_functor_wrappers_.size()
         + element_and_intersection_functor_wrappers_.size())
        > 0) {
      walk_range(elements(grid_view_));
    } // only do something, if we have to

    // finalize functors
    finalize();

    if (clear_functors)
      clear();
  } // ... walk(...)

  void clear()
  {
    element_functor_wrappers_.clear();
    intersection_functor_wrappers_.clear();
    element_and_intersection_functor_wrappers_.clear();
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
    if ((element_functor_wrappers_.size() + intersection_functor_wrappers_.size()
         + element_and_intersection_functor_wrappers_.size())
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
  template <class ElementRange>
  void walk_range(const ElementRange& element_range)
  {
#ifdef __INTEL_COMPILER
    const auto it_end = element_range.end();
    for (auto it = element_range.begin(); it != it_end; ++it) {
      const ElementType& element = *it;
#else
    for (const ElementType& element : element_range) {
#endif
      // apply element functors
      apply_local(element);

      // only walk the intersections, if there are codim1 functors present
      if ((intersection_functor_wrappers_.size() + element_and_intersection_functor_wrappers_.size()) > 0) {
        // Do not use intersections(...) here, since that does not work for a SubdomainGridPart which is based on
        // alugrid and then wrapped as a grid view (see also https://github.com/dune-community/dune-xt-grid/issues/26)
        const auto intersection_it_end = grid_view_.iend(element);
        for (auto intersection_it = grid_view_.ibegin(element); intersection_it != intersection_it_end;
             ++intersection_it) {
          const auto& intersection = *intersection_it;
          if (intersection.neighbor()) {
            const auto neighbor = intersection.outside();
            apply_local(intersection, element, neighbor);
          } else
            apply_local(intersection, element, element);
        } // walk the intersections
      } // only walk the intersections, if there are codim1 functors present
    } // .. walk elements
  } // ... walk_range(...)

  friend class internal::WalkerWrapper<GridViewType>;

  GridViewType grid_view_;
  std::list<std::shared_ptr<internal::ElementFunctorWrapper<GridViewType>>> element_functor_wrappers_;
  std::list<std::shared_ptr<internal::IntersectionFunctorWrapper<GridViewType>>> intersection_functor_wrappers_;
  std::list<std::shared_ptr<internal::ElementAndIntersectionFunctorWrapper<GridViewType>>>
      element_and_intersection_functor_wrappers_;
}; // class Walker


template <class GL>
Walker<GL> make_walker(GL grid_view)
{
  return Walker<GL>(grid_view);
}


} // namespace Grid
} // namespace XT
} // namespace Dune


#include "walker.lib.hh"


#endif // DUNE_XT_GRID_WALKER_HH
