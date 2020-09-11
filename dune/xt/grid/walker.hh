// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2019)
//   René Fritze     (2014 - 2019)
//   Tobias Leibner  (2014 - 2015, 2017 - 2020)

#ifndef DUNE_XT_GRID_WALKER_HH
#define DUNE_XT_GRID_WALKER_HH

#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

#if HAVE_TBB
#  include <tbb/blocked_range.h>
#  include <tbb/parallel_reduce.h>
#endif

#include <dune/common/version.hh>

#include <dune/grid/common/rangegenerators.hh>
#if HAVE_TBB
#  include <dune/xt/grid/parallel/partitioning/ranged.hh>
#endif

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/parallel/threadmanager.hh>
#include <dune/xt/common/parallel/threadstorage.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/grid/filters.hh>
#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/grid/functors/generic.hh>
#include <dune/xt/grid/type_traits.hh>

namespace Dune {
namespace XT {
namespace Grid {


// forward
template <class Gl>
class Walker;


namespace internal {


/**
 * \brief To be used within the Walker as internal storage type.
 *
 * \note Most likely you do not want to use this class directly, but instead append() an ElementFunctor to a Walker.
 *
 * \sa ElementFunctor
 * \sa Walker
 */
template <class GL>
class ElementFunctorWrapper
{
  static_assert(is_layer<GL>::value, "");
  using ThisType = ElementFunctorWrapper;

public:
  using FunctorType = ElementFunctor<GL>;
  using FilterType = ElementFilter<GL>;

  ElementFunctorWrapper(FunctorType& functr, const FilterType& filtr)
    : functor_(functr.copy())
    , filter_(filtr.copy())
  {}

  ElementFunctorWrapper(const ThisType& other)
    : functor_(other.functor_->copy())
    , filter_(other.filter_->copy())
  {}

  const FilterType& filter() const
  {
    return *filter_;
  }

  FunctorType& functor()
  {
    return *functor_;
  }

private:
  std::unique_ptr<FunctorType> functor_;
  const std::unique_ptr<const FilterType> filter_;
}; // class ElementFunctorWrapper


/**
 * \brief To be used within the \sa Walker as internal storage type.
 * \note  Most likely you do not want to use this class directly, but instead append() an \sa IntersectionFunctor to a
 *        Walker.
 */
template <class GL>
class IntersectionFunctorWrapper
{
  static_assert(is_layer<GL>::value, "");
  using ThisType = IntersectionFunctorWrapper;

public:
  using FunctorType = IntersectionFunctor<GL>;
  using FilterType = IntersectionFilter<GL>;

  IntersectionFunctorWrapper(FunctorType& functr, const FilterType& filtr)
    : functor_(functr.copy())
    , filter_(filtr.copy())
  {}

  IntersectionFunctorWrapper(const ThisType& other)
    : functor_(other.functor_->copy())
    , filter_(other.filter_->copy())
  {}

  const FilterType& filter() const
  {
    return *filter_;
  }

  FunctorType& functor()
  {
    return *functor_;
  }

private:
  std::unique_ptr<FunctorType> functor_;
  const std::unique_ptr<const FilterType> filter_;
}; // class IntersectionFunctorWrapper


/**
 * \brief To be used within the \sa Walker as internal storage type.
 * \note  Most likely you do not want to use this class directly, but instead append() an \sa
 *        ElementAndIntersectionFunctor to a Walker.
 */
template <class GL>
class ElementAndIntersectionFunctorWrapper
{
  static_assert(is_layer<GL>::value, "");
  using ThisType = ElementAndIntersectionFunctorWrapper;

public:
  using FunctorType = ElementAndIntersectionFunctor<GL>;
  using ElementFilterType = ElementFilter<GL>;
  using IntersectionFilterType = IntersectionFilter<GL>;

  ElementAndIntersectionFunctorWrapper(FunctorType& functr,
                                       const ElementFilterType& element_filtr,
                                       const IntersectionFilterType& intersection_filtr)
    : functor_(functr.copy())
    , element_filter_(element_filtr.copy())
    , intersection_filter_(intersection_filtr.copy())
  {}

  ElementAndIntersectionFunctorWrapper(const ThisType& other)
    : functor_(other.functor_->copy())
    , element_filter_(other.element_filter_->copy())
    , intersection_filter_(other.intersection_filter_->copy())
  {}

  const ElementFilterType& element_filter() const
  {
    return *element_filter_;
  }

  const IntersectionFilterType& intersection_filter() const
  {
    return *intersection_filter_;
  }

  FunctorType& functor()
  {
    return *functor_;
  }

private:
  std::unique_ptr<FunctorType> functor_;
  const std::unique_ptr<const ElementFilterType> element_filter_;
  const std::unique_ptr<const IntersectionFilterType> intersection_filter_;
}; // class ElementAndIntersectionFunctorWrapper


} // namespace internal


template <class GV>
class Walker : public ElementAndIntersectionFunctor<GV>
{
protected:
  using BaseType = ElementAndIntersectionFunctor<GV>;

private:
  using ThisType = Walker;

public:
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

private:
  using GenericElementApplyFunctionType = std::function<void(const ElementType&)>;
  using GenericIntersectionApplyFunctionType =
      std::function<void(const IntersectionType&, const ElementType&, const ElementType&)>;
  using GenericElementFilterFunctionType = std::function<bool(const GV&, const ElementType&)>;
  using GenericIntersectionFilterFunctionType = std::function<bool(const GV&, const IntersectionType&)>;
  using VoidFunctionType = std::function<void()>;

  template <typename WrapperType, class... Args>
  void emplace_all(std::list<WrapperType>& storage,
                   Common::PerThreadValue<std::list<WrapperType>>& thread_storage,
                   Args&&... args)
  {
    storage.emplace_back(std::forward<Args>(args)...);
    thread_storage = Common::PerThreadValue<std::list<WrapperType>>(storage);
  }

  void reinitialize_thread_storage()
  {
    element_functor_wrappers_ = Common::PerThreadValue<std::list<internal::ElementFunctorWrapper<GridViewType>>>(
        stored_element_functor_wrappers_);
    intersection_functor_wrappers_ =
        Common::PerThreadValue<std::list<internal::IntersectionFunctorWrapper<GridViewType>>>(
            stored_intersection_functor_wrappers_);
    element_and_intersection_functor_wrappers_ =
        Common::PerThreadValue<std::list<internal::ElementAndIntersectionFunctorWrapper<GridViewType>>>(
            stored_element_and_intersection_functor_wrappers_);
  }

public:
  explicit Walker(GridViewType grd_vw)
    : grid_view_(grd_vw)
  {}

  Walker(const ThisType& other)
    : BaseType(other)
    , grid_view_(other.grid_view_)
    , stored_element_functor_wrappers_(other.stored_element_functor_wrappers_)
    , stored_intersection_functor_wrappers_(other.stored_intersection_functor_wrappers_)
    , stored_element_and_intersection_functor_wrappers_(other.stored_element_and_intersection_functor_wrappers_)
  {
    reinitialize_thread_storage();
  }

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

  ThisType& append(ElementFunctor<GV>& functor, const ElementFilter<GV>& filter = ApplyOn::AllElements<GV>())
  {
    emplace_all(stored_element_functor_wrappers_, element_functor_wrappers_, functor, filter);
    return *this;
  }

  ThisType& append(ElementFunctor<GV>*&& functor, const ElementFilter<GV>& filter = ApplyOn::AllElements<GV>())
  {
    emplace_all(stored_element_functor_wrappers_, element_functor_wrappers_, *functor, filter);
    delete functor;
    return *this;
  }

  ThisType& append(ElementFunctor<GV>& functor, GenericElementFilterFunctionType filter)
  {
    emplace_all(stored_element_functor_wrappers_,
                element_functor_wrappers_,
                functor,
                ApplyOn::GenericFilteredElements<GV>(filter));
    return *this;
  }

  ThisType& append(ElementFunctor<GV>*&& functor, GenericElementFilterFunctionType filter)
  {
    emplace_all(stored_element_functor_wrappers_,
                element_functor_wrappers_,
                *functor,
                ApplyOn::GenericFilteredElements<GV>(filter));
    delete functor;
    return *this;
  }

  /**
   * \}
   * \name These methods can be used to append an element lambda expression.
   * \{
   */

  ThisType& append(VoidFunctionType prepare_func,
                   GenericElementApplyFunctionType apply_func,
                   VoidFunctionType finalize_func,
                   const ElementFilter<GV>& filter = ApplyOn::AllElements<GV>())
  {
    return this->append(new GenericElementFunctor<GV>(prepare_func, apply_func, finalize_func), filter);
  }

  ThisType& append(VoidFunctionType prepare_func,
                   GenericElementApplyFunctionType apply_func,
                   VoidFunctionType finalize_func,
                   GenericElementFilterFunctionType filter)
  {
    return this->append(new GenericElementFunctor<GV>(prepare_func, apply_func, finalize_func),
                        ApplyOn::GenericFilteredElements<GV>(filter));
  }

  /**
   * \}
   * \name These methods can be used to append an \sa IntersectionFunctor.
   * \{
   */

  ThisType& append(IntersectionFunctor<GV>& functor,
                   const IntersectionFilter<GV>& filter = ApplyOn::AllIntersections<GV>())
  {
    emplace_all(stored_intersection_functor_wrappers_, intersection_functor_wrappers_, functor, filter);
    return *this;
  }

  ThisType& append(IntersectionFunctor<GV>*&& functor,
                   const IntersectionFilter<GV>& filter = ApplyOn::AllIntersections<GV>())
  {
    emplace_all(stored_intersection_functor_wrappers_, intersection_functor_wrappers_, *functor, filter);
    delete functor;
    return *this;
  }

  ThisType& append(IntersectionFunctor<GV>& functor, GenericIntersectionFilterFunctionType filter)
  {
    emplace_all(stored_intersection_functor_wrappers_,
                intersection_functor_wrappers_,
                functor,
                ApplyOn::GenericFilteredIntersections<GV>(filter));
    return *this;
  }

  ThisType& append(IntersectionFunctor<GV>*&& functor, GenericIntersectionFilterFunctionType filter)
  {
    emplace_all(stored_intersection_functor_wrappers_,
                intersection_functor_wrappers_,
                *functor,
                ApplyOn::GenericFilteredIntersections<GV>(filter));
    delete functor;
    return *this;
  }

  /**
   * \}
   * \name These methods can be used to append an intersection lambda expression.
   * \{
   */

  ThisType& append(VoidFunctionType prepare_func,
                   GenericIntersectionApplyFunctionType apply_func,
                   VoidFunctionType finalize_func,
                   const IntersectionFilter<GV>& filter = ApplyOn::AllIntersections<GV>())
  {
    return this->append(new GenericIntersectionFunctor<GV>(prepare_func, apply_func, finalize_func), filter);
  }

  ThisType& append(VoidFunctionType prepare_func,
                   GenericIntersectionApplyFunctionType apply_func,
                   VoidFunctionType finalize_func,
                   GenericIntersectionFilterFunctionType filter)
  {
    return this->append(new GenericIntersectionFunctor<GV>(prepare_func, apply_func, finalize_func),
                        ApplyOn::GenericFilteredIntersections<GV>(filter));
  }

  /**
   * \}
   * \name These methods can be used to append an \sa ElementAndIntersectionFunctor.
   * \{
   */

  ThisType& append(ElementAndIntersectionFunctor<GV>& functor,
                   const IntersectionFilter<GV>& intersection_filter = ApplyOn::AllIntersections<GV>(),
                   const ElementFilter<GV>& element_filter = ApplyOn::AllElements<GV>())
  {
    if (&functor == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    emplace_all(stored_element_and_intersection_functor_wrappers_,
                element_and_intersection_functor_wrappers_,
                functor,
                element_filter,
                intersection_filter);
    return *this;
  }

  ThisType& append(ElementAndIntersectionFunctor<GV>*&& functor,
                   const IntersectionFilter<GV>& intersection_filter = ApplyOn::AllIntersections<GV>(),
                   const ElementFilter<GV>& element_filter = ApplyOn::AllElements<GV>())
  {
    if (functor == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    emplace_all(stored_element_and_intersection_functor_wrappers_,
                element_and_intersection_functor_wrappers_,
                *functor,
                element_filter,
                intersection_filter);
    delete functor;
    return *this;
  }

  ThisType& append(ElementAndIntersectionFunctor<GV>& functor,
                   GenericElementFilterFunctionType element_filter,
                   GenericIntersectionFilterFunctionType intersection_filter)
  {
    if (&functor == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    emplace_all(stored_element_and_intersection_functor_wrappers_,
                element_and_intersection_functor_wrappers_,
                functor,
                ApplyOn::GenericFilteredElements<GV>(element_filter),
                ApplyOn::GenericFilteredIntersections<GV>(intersection_filter));
    return *this;
  }

  ThisType& append(ElementAndIntersectionFunctor<GV>*&& functor,
                   GenericElementFilterFunctionType element_filter,
                   GenericIntersectionFilterFunctionType intersection_filter)
  {
    if (functor == this)
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong, "Do not append a Walker to itself!");
    emplace_all(stored_element_and_intersection_functor_wrappers_,
                element_and_intersection_functor_wrappers_,
                *functor,
                ApplyOn::GenericFilteredElements<GV>(element_filter),
                ApplyOn::GenericFilteredIntersections<GV>(intersection_filter));
    delete functor;
    return *this;
  }

  /**
   * \}
   * \name These methods can be used to append element and intersection lambda expressions.
   * \{
   */

  ThisType& append(VoidFunctionType prepare_func,
                   GenericElementApplyFunctionType element_apply_on,
                   GenericIntersectionApplyFunctionType intersection_apply_on,
                   VoidFunctionType finalize_func,
                   const ElementFilter<GV>& element_filter = ApplyOn::AllElements<GV>(),
                   const IntersectionFilter<GV>& intersection_filter = ApplyOn::AllIntersections<GV>())
  {
    return this->append(new GenericElementAndIntersectionFunctor<GV>(
                            prepare_func, element_apply_on, intersection_apply_on, finalize_func),
                        element_filter,
                        intersection_filter);
  }

  ThisType& append(VoidFunctionType prepare_func,
                   GenericElementApplyFunctionType element_apply_on,
                   GenericIntersectionApplyFunctionType intersection_apply_on,
                   VoidFunctionType finalize_func,
                   GenericElementFilterFunctionType element_filter,
                   GenericIntersectionFilterFunctionType intersection_filter)
  {
    return this->append(new GenericElementAndIntersectionFunctor<GV>(
                            prepare_func, element_apply_on, intersection_apply_on, finalize_func),
                        ApplyOn::GenericFilteredElements<GV>(element_filter),
                        ApplyOn::GenericFilteredIntersections<GV>(intersection_filter));
  }

  /**
   * \}
   * \name These methods are required by ElementAndIntersectionFunctor.
   * \{
   */

  virtual void prepare() override
  {
    auto prep = [](auto& wrapper_list) {
      for (auto&& wrapper : wrapper_list)
        wrapper.functor().prepare();
    };
    prep(stored_element_functor_wrappers_);
    prep(stored_intersection_functor_wrappers_);
    prep(stored_element_and_intersection_functor_wrappers_);
    // prepare is called in a single thread, so make sure all functors in the thread are also prepared
    reinitialize_thread_storage();
  } // ... prepare()

  // prepare calling thread
  void prepare_thread()
  {
    auto prep = [](auto& wrapper_list) {
      for (auto&& wrapper : wrapper_list)
        wrapper.functor().prepare();
    };
    prep(*element_functor_wrappers_);
    prep(*intersection_functor_wrappers_);
    prep(*element_and_intersection_functor_wrappers_);
  }

  virtual void apply_local(const ElementType& element) override
  {
    apply_local(element, *element_functor_wrappers_, *element_and_intersection_functor_wrappers_);
  } // ... apply_local(...)

  void apply_local(const ElementType& element,
                   std::list<internal::ElementFunctorWrapper<GridViewType>>& element_functor_wrappers,
                   std::list<internal::ElementAndIntersectionFunctorWrapper<GridViewType>>&
                       element_and_intersection_functor_wrappers)
  {
    for (auto&& wrapper : element_functor_wrappers) {
      if (wrapper.filter().contains(grid_view_, element))
        wrapper.functor().apply_local(element);
    }
    for (auto&& wrapper : element_and_intersection_functor_wrappers) {
      if (wrapper.element_filter().contains(grid_view_, element))
        wrapper.functor().apply_local(element);
    }
  } // ... apply_local(...)

  virtual void apply_local(const IntersectionType& intersection,
                           const ElementType& inside_element,
                           const ElementType& outside_element) override
  {
    apply_local(intersection,
                inside_element,
                outside_element,
                *intersection_functor_wrappers_,
                *element_and_intersection_functor_wrappers_);
  } // ... apply_local(...)

  virtual void apply_local(const IntersectionType& intersection,
                           const ElementType& inside_element,
                           const ElementType& outside_element,
                           std::list<internal::IntersectionFunctorWrapper<GridViewType>>& intersection_functor_wrappers,
                           std::list<internal::ElementAndIntersectionFunctorWrapper<GridViewType>>&
                               element_and_intersection_functor_wrappers)
  {
    for (auto&& wrapper : intersection_functor_wrappers) {
      if (wrapper.filter().contains(grid_view_, intersection))
        wrapper.functor().apply_local(intersection, inside_element, outside_element);
    }
    for (auto&& wrapper : element_and_intersection_functor_wrappers) {
      if (wrapper.intersection_filter().contains(grid_view_, intersection))
        wrapper.functor().apply_local(intersection, inside_element, outside_element);
    }
  } // ... apply_local(...)

  // finalize all threads
  virtual void finalize() override
  {
    auto fin = [](auto& per_thread_value) {
      for (auto&& list : per_thread_value) {
        for (auto&& wrapper : list)
          wrapper.functor().finalize();
      }
    };
    fin(element_functor_wrappers_);
    fin(intersection_functor_wrappers_);
    fin(element_and_intersection_functor_wrappers_);
  } // ... finalize()

  /**
   * \}
   */

  void walk([[maybe_unused]] const bool use_tbb = false, const bool clear_functors = true)
  {
#if HAVE_TBB
    if (use_tbb) {
      const auto num_partitions =
          DXTC_CONFIG_GET("threading.partition_factor", 1u) * XT::Common::threadManager().current_threads();
      RangedPartitioning<GridViewType, 0> partitioning(grid_view_, num_partitions);
      this->walk(partitioning, clear_functors);
      return;
    }
#endif
    // prepare functors
    prepare();

    // only do something, if we have to
    if ((element_functor_wrappers_->size() + intersection_functor_wrappers_->size()
         + element_and_intersection_functor_wrappers_->size())
        > 0) {
      walk_range(elements(grid_view_));
    }

    // finalize functors
    finalize();

    if (clear_functors)
      clear();
  } // ... walk(...)

  void clear()
  {
    stored_element_functor_wrappers_.clear();
    stored_intersection_functor_wrappers_.clear();
    stored_element_and_intersection_functor_wrappers_.clear();
    reinitialize_thread_storage();
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
    Body(WalkerType& walker, const PartioningType& partitioning)
      : walker_(walker)
      , partitioning_(partitioning)
    {}

    Body(Body& other, tbb::split /*split*/)
      : walker_(other.walker_)
      , partitioning_(other.partitioning_)
    {}

    void operator()(const tbb::blocked_range<std::size_t>& range) const
    {
      // for all partitions in tbb-range
      for (std::size_t p = range.begin(); p != range.end(); ++p) {
        auto partition = partitioning_.partition(p);
        walker_.walk_range(partition);
      }
    }

    void join(Body& /*other*/) {}

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
    if ((element_functor_wrappers_->size() + intersection_functor_wrappers_->size()
         + element_and_intersection_functor_wrappers_->size())
        > 0) {
      // no actual SMP walk, use range as is
      walk_range(partitioning.everything());
    }

    // finalize functors
    finalize();
    clear();
  } // ... tbb_walk(...)
#endif // HAVE_TBB

  template <class ElementRange>
  void walk_range(const ElementRange& element_range)
  {
    auto& element_functor_wrappers = *element_functor_wrappers_;
    auto& intersection_functor_wrappers = *intersection_functor_wrappers_;
    auto& element_and_intersection_functor_wrappers = *element_and_intersection_functor_wrappers_;
#ifdef __INTEL_COMPILER
    const auto it_end = element_range.end();
    for (auto it = element_range.begin(); it != it_end; ++it) {
      const ElementType& element = *it;
#else
    for (const ElementType& element : element_range) {
#endif
      // apply element functors
      apply_local(element, element_functor_wrappers, element_and_intersection_functor_wrappers);

      // only walk the intersections, if there are codim1 functors present
      if ((intersection_functor_wrappers.size() + element_and_intersection_functor_wrappers.size()) > 0) {
        for (auto&& intersection : intersections(grid_view_, element)) {
          if (intersection.neighbor()) {
            const auto neighbor = intersection.outside();
            apply_local(intersection,
                        element,
                        neighbor,
                        intersection_functor_wrappers,
                        element_and_intersection_functor_wrappers);
          } else
            apply_local(intersection,
                        element,
                        element,
                        intersection_functor_wrappers,
                        element_and_intersection_functor_wrappers);
        } // walk the intersections
      } // only walk the intersections, if there are codim1 functors present
    } // .. walk elements
  } // ... walk_range(...)

protected:
  GridViewType grid_view_;
  // We want each thread to have its own copy of each functor. However, as we do not know in advance how many different
  // threads we will have (even if DXTC_CONFIG["threading.max_count"] is set, there may be only max_threads at a time,
  // but we cannot guarantee that there will only be max_threads over the whole runtime of our program), we store the
  // wrappers in a list that is shared by all threads and copy that list to the thread once a new thread tries to access
  // the functors.
  std::list<internal::ElementFunctorWrapper<GridViewType>> stored_element_functor_wrappers_;
  std::list<internal::IntersectionFunctorWrapper<GridViewType>> stored_intersection_functor_wrappers_;
  std::list<internal::ElementAndIntersectionFunctorWrapper<GridViewType>>
      stored_element_and_intersection_functor_wrappers_;
  Common::PerThreadValue<std::list<internal::ElementFunctorWrapper<GridViewType>>> element_functor_wrappers_;
  Common::PerThreadValue<std::list<internal::IntersectionFunctorWrapper<GridViewType>>> intersection_functor_wrappers_;
  Common::PerThreadValue<std::list<internal::ElementAndIntersectionFunctorWrapper<GridViewType>>>
      element_and_intersection_functor_wrappers_;
}; // class Walker

template <class GV>
Walker<GV> make_walker(GV grid_view)
{
  return Walker<GV>(grid_view);
}


} // namespace Grid
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_GRID_WALKER_HH
