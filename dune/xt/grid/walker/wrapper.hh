// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   Rene Milk       (2014 - 2016, 2018)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_GRID_WALKER_WRAPPER_HH
#define DUNE_XT_GRID_WALKER_WRAPPER_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/grid/type_traits.hh>

#include "filters.hh"

namespace Dune {
namespace XT {
namespace Grid {

// forward
template <class Gl>
class Walker;


namespace internal {


/**
 * \brief To be used within the \sa Walker as internal storage type.
 * \note  Most likely you do not want to use this class directly, but instead append() an \sa ElementFunctor to a
 *        Walker.
 */
template <class GL>
class ElementFunctorWrapper
{
  static_assert(is_layer<GL>::value, "");
  using ThisType = ElementFunctorWrapper<GL>;

public:
  using FunctorType = ElementFunctor<GL>;
  using FilterType = ElementFilter<GL>;

  ElementFunctorWrapper(FunctorType& functr, const FilterType& filtr)
    : functor_(functr.copy())
    , filter_(filtr.copy())
  {
  }

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
  using ThisType = IntersectionFunctorWrapper<GL>;

public:
  using FunctorType = IntersectionFunctor<GL>;
  using FilterType = IntersectionFilter<GL>;

  IntersectionFunctorWrapper(FunctorType& functr, const FilterType& filtr)
    : functor_(functr.copy())
    , filter_(filtr.copy())
  {
  }

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
  using ThisType = ElementAndIntersectionFunctorWrapper<GL>;

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
  {
  }

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
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_WALKER_WRAPPER_HH
