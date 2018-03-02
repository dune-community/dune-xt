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
#include <dune/xt/common/parallel/threadstorage.hh>
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

static const std::function<void()> dxt_void_noop = [](...) {};

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

private:
  using IntersectionElementElementFunction =
      std::function<void(const IntersectionType&, const ElementType&, const ElementType&)>;
  using ElementFunction = std::function<void(const ElementType&)>;
  using ViewElementFunction = std::function<bool(const GL&, const ElementType&)>;
  using ViewIntersectionFunction = std::function<bool(const GL&, const IntersectionType&)>;
  using VoidFunction = std::function<void()>;

  template <typename WrapperType, class... Args>
  void emplace_all(Common::PerThreadValue<std::list<std::shared_ptr<WrapperType>>>& thread_storage, Args&&... args)
  {
    for (auto&& local_list : thread_storage) {
      local_list->emplace_back(new WrapperType(std::forward<Args>(args)...));
    }
  }

public:
  explicit Walker(GridViewType grd_vw)
    : grid_view_(grd_vw)
  {
  }

  Walker(const ThisType& other)
    : BaseType()
    , grid_view_(other.grid_view_)
  {
    // Since all Common::PerThreadValue are created with the same size given by the global singleton threadManage(),
    // we just assume they are of the same size!
    // Copy the element functors ...
    auto zip_emplace = [](auto& target_thread_ctr, const auto& source_thread_ctr, auto& gen_function) {
      auto target_wrappers_it = target_thread_ctr.begin();
      auto source_wrappers_it = source_thread_ctr.begin();
      for (; target_wrappers_it != target_thread_ctr.end() && source_wrappers_it != source_thread_ctr.end();
           ++target_wrappers_it, ++source_wrappers_it) {
        auto& target__wrappers = *(*target_wrappers_it);
        for (auto&& source_wrapper : **source_wrappers_it) {
          target__wrappers.emplace_back(gen_function(*source_wrapper));
        }
      }
    };
    auto el_wrapper = [](auto&& wrapper) {
      return new typename std::decay<decltype(wrapper)>::type(wrapper.functor(), wrapper.filter());
    };
    zip_emplace(element_functor_wrappers_, other.element_functor_wrappers_, el_wrapper);
    zip_emplace(intersection_functor_wrappers_, other.intersection_functor_wrappers_, el_wrapper);

    auto elint_wrapper = [](auto&& wrapper) {
      return new typename std::decay<decltype(wrapper)>::type(
          wrapper.functor(), wrapper.element_filter(), wrapper.intersection_filter());
    };
    zip_emplace(
        element_and_intersection_functor_wrappers_, other.element_and_intersection_functor_wrappers_, elint_wrapper);
  } // Walker(...)

  Walker(ThisType&& source) = default;

  virtual ~Walker() = default;

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

  ThisType& append(ElementFunctor<GL>& functor, const ElementFilter<GL>& filter = ApplyOn::AllElements<GL>())
  {
    emplace_all(element_functor_wrappers_, functor, filter);
    return *this;
  }

  ThisType& append(ElementFunctor<GL>*&& functor, const ElementFilter<GL>& filter = ApplyOn::AllElements<GL>())
  {
    emplace_all(element_functor_wrappers_, *functor, filter);
    delete functor;
    return *this;
  }

  ThisType& append(ElementFunctor<GL>& functor, ViewElementFunction element_filter)
  {
    emplace_all(element_functor_wrappers_, functor, ApplyOn::LambdaFilteredElements<GL>(element_filter));
    return *this;
  }

  ThisType& append(ElementFunctor<GL>*&& functor, ViewElementFunction element_filter)
  {
    emplace_all(element_functor_wrappers_, *functor, ApplyOn::LambdaFilteredElements<GL>(element_filter));
    delete functor;
    return *this;
  }

  /**
   * \}
   * \name These methods can be used to append an element lambda expression.
   * \{
   */

  ThisType& append(ElementFunction apply_lambda,
                   const ElementFilter<GL>& filter = ApplyOn::AllElements<GL>(),
                   VoidFunction prepare_lambda = dxt_void_noop,
                   VoidFunction finalize_lambda = dxt_void_noop)
  {
    return this->append(new ElementLambdaFunctor<GL>(apply_lambda, prepare_lambda, finalize_lambda), filter);
  }

  ThisType& append(ElementFunction apply_lambda,
                   ViewElementFunction filter,
                   VoidFunction prepare_lambda = dxt_void_noop,
                   VoidFunction finalize_lambda = dxt_void_noop)
  {
    return this->append(new ElementLambdaFunctor<GL>(apply_lambda, prepare_lambda, finalize_lambda),
                        ApplyOn::LambdaFilteredElements<GL>(filter));
  }

  /**
   * \}
   * \name These methods can be used to append an \sa IntersectionFunctor.
   * \{
   */

  ThisType& append(IntersectionFunctor<GL>& functor,
                   const IntersectionFilter<GL>& filter = ApplyOn::AllIntersections<GL>())
  {
    emplace_all(intersection_functor_wrappers_, functor, filter);
    return *this;
  }

  ThisType& append(IntersectionFunctor<GL>*&& functor,
                   const IntersectionFilter<GL>& filter = ApplyOn::AllIntersections<GL>())
  {
    emplace_all(intersection_functor_wrappers_, *functor, filter);
    delete functor;
    return *this;
  }

  ThisType& append(IntersectionFunctor<GL>& functor, ViewIntersectionFunction filter)
  {
    emplace_all(intersection_functor_wrappers_, functor, ApplyOn::LambdaFilteredIntersections<GL>(filter));
    return *this;
  }

  ThisType& append(IntersectionFunctor<GL>*&& functor, ViewIntersectionFunction filter)
  {
    emplace_all(intersection_functor_wrappers_, *functor, ApplyOn::LambdaFilteredIntersections<GL>(filter));
    delete functor;
    return *this;
  }

  /**
   * \}
   * \name These methods can be used to append an intersection lambda expression.
   * \{
   */

  ThisType& append(IntersectionElementElementFunction apply_lambda,
                   const IntersectionFilter<GL>& filter = ApplyOn::AllIntersections<GL>(),
                   VoidFunction prepare_lambda = dxt_void_noop,
                   VoidFunction finalize_lambda = dxt_void_noop)
  {
    return this->append(new IntersectionLambdaFunctor<GL>(apply_lambda, prepare_lambda, finalize_lambda), filter);
  }

  ThisType& append(IntersectionElementElementFunction apply_lambda,
                   ViewIntersectionFunction filter,
                   VoidFunction prepare_lambda = dxt_void_noop,
                   VoidFunction finalize_lambda = dxt_void_noop)
  {
    return this->append(new IntersectionLambdaFunctor<GL>(apply_lambda, prepare_lambda, finalize_lambda),
                        ApplyOn::LambdaFilteredIntersections<GL>(filter));
  }

  /**
   * \}
   * \name These methods can be used to append an \sa ElementAndIntersectionFunctor.
   * \{
   */

  ThisType& append(ElementAndIntersectionFunctor<GL>& functor,
                   const IntersectionFilter<GL>& intersection_filter = ApplyOn::AllIntersections<GL>(),
                   const ElementFilter<GL>& element_filter = ApplyOn::AllElements<GL>())
  {
    if (&functor == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    emplace_all(element_and_intersection_functor_wrappers_, functor, element_filter, intersection_filter);
    return *this;
  }

  ThisType& append(ElementAndIntersectionFunctor<GL>*&& functor,
                   const IntersectionFilter<GL>& intersection_filter = ApplyOn::AllIntersections<GL>(),
                   const ElementFilter<GL>& element_filter = ApplyOn::AllElements<GL>())
  {
    if (functor == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    emplace_all(element_and_intersection_functor_wrappers_, *functor, element_filter, intersection_filter);
    delete functor;
    return *this;
  }

  ThisType& append(ElementAndIntersectionFunctor<GL>& functor,
                   ViewElementFunction element_filter,
                   ViewIntersectionFunction intersection_filter)
  {
    if (&functor == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    emplace_all(element_and_intersection_functor_wrappers_,
                functor,
                ApplyOn::LambdaFilteredElements<GL>(element_filter),
                ApplyOn::LambdaFilteredIntersections<GL>(intersection_filter));
    return *this;
  }

  ThisType& append(ElementAndIntersectionFunctor<GL>*&& functor,
                   ViewElementFunction element_filter,
                   ViewIntersectionFunction intersection_filter)
  {
    if (functor == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    emplace_all(element_and_intersection_functor_wrappers_,
                *functor,
                ApplyOn::LambdaFilteredElements<GL>(element_filter),
                ApplyOn::LambdaFilteredIntersections<GL>(intersection_filter));
    delete functor;
    return *this;
  }

  /**
   * \}
   * \name These methods can be used to append element and intersection lambda expressions.
   * \{
   */

  ThisType& append(ElementFunction element_apply_on,
                   IntersectionElementElementFunction intersection_apply_on,
                   ViewElementFunction element_filter,
                   ViewIntersectionFunction intersection_filter)
  {
    return this->append(new ElementAndIntersectionLambdaFunctor<GL>(element_apply_on, intersection_apply_on),
                        ApplyOn::LambdaFilteredElements<GL>(element_filter),
                        ApplyOn::LambdaFilteredIntersections<GL>(intersection_filter));
  }

  ThisType& append(ElementFunction element_apply_on,
                   IntersectionElementElementFunction intersection_apply_on,
                   const ElementFilter<GL>& element_filter = ApplyOn::AllElements<GL>(),
                   const IntersectionFilter<GL>& intersection_filter = ApplyOn::AllIntersections<GL>(),
                   VoidFunction prepare_lambda = dxt_void_noop,
                   VoidFunction finalize_lambda = dxt_void_noop)
  {
    return this->append(new ElementAndIntersectionLambdaFunctor<GL>(
                            element_apply_on, intersection_apply_on, prepare_lambda, finalize_lambda),
                        element_filter,
                        intersection_filter);
  }

  ThisType& append(ElementFunction element_apply_on,
                   IntersectionElementElementFunction intersection_apply_on,
                   ViewElementFunction element_filter,
                   ViewIntersectionFunction intersection_filter,
                   VoidFunction prepare_lambda = dxt_void_noop,
                   VoidFunction finalize_lambda = dxt_void_noop)
  {
    return this->append(new ElementAndIntersectionLambdaFunctor<GL>(
                            element_apply_on, intersection_apply_on, prepare_lambda, finalize_lambda),
                        ApplyOn::LambdaFilteredElements<GL>(element_filter),
                        ApplyOn::LambdaFilteredIntersections<GL>(intersection_filter));
  }

  /**
   * \}
   * \name These methods are required by ElementAndIntersectionFunctor.
   * \{
   */

  virtual void prepare() override
  {
    auto prep = [](auto& pt) {
      for (auto&& list_ptr : pt) {
        for (auto&& wrapper : *list_ptr)
          wrapper->functor().prepare();
      }
    };
    prep(element_functor_wrappers_);
    prep(intersection_functor_wrappers_);
    prep(element_and_intersection_functor_wrappers_);
  } // ... prepare()

  virtual void apply_local(const ElementType& element) override
  {
    for (auto&& wrapper : *element_functor_wrappers_) {
      if (wrapper->filter().contains(grid_view_, element))
        wrapper->functor().apply_local(element);
    }
    for (auto&& wrapper : *element_and_intersection_functor_wrappers_) {
      if (wrapper->element_filter().contains(grid_view_, element))
        wrapper->functor().apply_local(element);
    }
  } // ... apply_local(...)

  virtual void apply_local(const IntersectionType& intersection,
                           const ElementType& inside_element,
                           const ElementType& outside_element) override
  {
    for (auto&& wrapper : *intersection_functor_wrappers_) {
      if (wrapper->filter().contains(grid_view_, intersection))
        wrapper->functor().apply_local(intersection, inside_element, outside_element);
    }
    for (auto&& wrapper : *element_and_intersection_functor_wrappers_) {
      if (wrapper->intersection_filter().contains(grid_view_, intersection))
        wrapper->functor().apply_local(intersection, inside_element, outside_element);
    }
  } // ... apply_local(...)

  virtual void finalize() override
  {
    auto fin = [](auto& pt) {
      for (auto&& list_ptr : pt) {
        for (auto&& wrapper : *list_ptr)
          wrapper->functor().finalize();
      }
    };
    fin(element_and_intersection_functor_wrappers_);
    fin(element_functor_wrappers_);
    fin(intersection_functor_wrappers_);
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
    if ((element_functor_wrappers_->size() + intersection_functor_wrappers_->size()
         + element_and_intersection_functor_wrappers_->size())
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
    auto clr = [](auto& pt) {
      for (auto&& list_ptr : pt) {
        list_ptr->clear();
      }
    };
    clr(element_functor_wrappers_);
    clr(intersection_functor_wrappers_);
    clr(element_and_intersection_functor_wrappers_);
  }

  BaseType* copy() override
  {
    return new ThisType(*this);
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
    if ((element_functor_wrappers_->size() + intersection_functor_wrappers_->size()
         + element_and_intersection_functor_wrappers_->size())
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

private:
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
      if ((intersection_functor_wrappers_->size() + element_and_intersection_functor_wrappers_->size()) > 0) {
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

  GridViewType grid_view_;
  Common::PerThreadValue<std::list<std::shared_ptr<internal::ElementFunctorWrapper<GridViewType>>>>
      element_functor_wrappers_;
  Common::PerThreadValue<std::list<std::shared_ptr<internal::IntersectionFunctorWrapper<GridViewType>>>>
      intersection_functor_wrappers_;
  Common::PerThreadValue<std::list<std::shared_ptr<internal::ElementAndIntersectionFunctorWrapper<GridViewType>>>>
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
