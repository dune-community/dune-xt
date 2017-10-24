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
 * \note  Most likely you do not want to use this class directly, but instead append() an \sa EntityFunctor to a
 *        Walker.
 */
template <class GL>
class EntityFunctorWrapper
{
  static_assert(is_layer<GL>::value, "");
  using ThisType = EntityFunctorWrapper<GL>;

public:
  using FunctorType = EntityFunctor<GL>;
  using FilterType = EntityFilter<GL>;

  /**
   * \attention Takes ownership of filtr_ptr, do not delete manually!
   */
  EntityFunctorWrapper(FunctorType& functr, const FilterType*&& filtr_ptr)
    : functor_(functr)
    , filter_(std::move(filtr_ptr))
  {
  }

  /**
   * \attention Takes ownership of functr_ptr, do not delete manually!
   * \attention Takes ownership of filtr_ptr, do not delete manually!
   */
  EntityFunctorWrapper(FunctorType*&& functr_ptr, const FilterType*&& filtr_ptr)
    : functor_(std::move(functr_ptr))
    , filter_(std::move(filtr_ptr))
  {
  }

  bool operator==(const ThisType& other) const
  {
    return this == &other;
  }

  bool operator!=(const ThisType& other) const
  {
    return this != other;
  }

  virtual const FilterType& filter() const
  {
    return filter_.access();
  }

  virtual FunctorType& functor()
  {
    return functor_.access();
  }

private:
  Common::StorageProvider<FunctorType> functor_;
  Common::ConstStorageProvider<FilterType> filter_;
}; // class EntityFunctorWrapper


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

  /**
   * \attention Takes ownership of filtr_ptr, do not delete manually!
   */
  IntersectionFunctorWrapper(FunctorType& functr, const FilterType*&& filtr_ptr)
    : functor_(functr)
    , filter_(std::move(filtr_ptr))
  {
  }

  /**
   * \attention Takes ownership of functr_ptr, do not delete manually!
   * \attention Takes ownership of filtr_ptr, do not delete manually!
   */
  IntersectionFunctorWrapper(FunctorType*&& functr_ptr, const FilterType*&& filtr_ptr)
    : functor_(std::move(functr_ptr))
    , filter_(std::move(filtr_ptr))
  {
  }

  bool operator==(const ThisType& other) const
  {
    return this == &other;
  }

  bool operator!=(const ThisType& other) const
  {
    return this != other;
  }

  virtual const FilterType& filter() const
  {
    return filter_.access();
  }

  virtual FunctorType& functor()
  {
    return functor_.access();
  }

private:
  Common::StorageProvider<FunctorType> functor_;
  Common::ConstStorageProvider<FilterType> filter_;
}; // class IntersectionFunctorWrapper


/**
 * \brief To be used within the \sa Walker as internal storage type.
 * \note  Most likely you do not want to use this class directly, but instead append() an \sa
 *        EntityAndIntersectionFunctor to a Walker.
 */
template <class GL>
class EntityAndIntersectionFunctorWrapper
{
  static_assert(is_layer<GL>::value, "");
  using ThisType = EntityAndIntersectionFunctorWrapper<GL>;

public:
  using FunctorType = EntityAndIntersectionFunctor<GL>;
  using EntityFilterType = EntityFilter<GL>;
  using IntersectionFilterType = IntersectionFilter<GL>;

  /**
   * \attention Takes ownership of entity_filtr_ptr, do not delete manually!
   * \attention Takes ownership of intersection_filtr_ptr, do not delete manually!
   */
  EntityAndIntersectionFunctorWrapper(FunctorType& functr,
                                      const EntityFilterType*&& entity_filtr_ptr,
                                      const IntersectionFilterType*&& intersection_filtr_ptr)
    : functor_(functr)
    , entity_filter_(std::move(entity_filtr_ptr))
    , intersection_filter_(std::move(intersection_filtr_ptr))
  {
  }

  /**
   * \attention Takes ownership of functr_ptr, do not delete manually!
   * \attention Takes ownership of entity_filtr_ptr, do not delete manually!
   * \attention Takes ownership of intersection_filtr_ptr, do not delete manually!
   */
  EntityAndIntersectionFunctorWrapper(FunctorType*&& functr_ptr,
                                      const EntityFilterType*&& entity_filtr_ptr,
                                      const IntersectionFilterType*&& intersection_filtr_ptr)
    : functor_(std::move(functr_ptr))
    , entity_filter_(std::move(entity_filtr_ptr))
    , intersection_filter_(std::move(intersection_filtr_ptr))
  {
  }

  bool operator==(const ThisType& other) const
  {
    return this == &other;
  }

  bool operator!=(const ThisType& other) const
  {
    return this != other;
  }

  virtual const EntityFilterType& entity_filter() const
  {
    return entity_filter_.access();
  }

  virtual const IntersectionFilterType& intersection_filter() const
  {
    return intersection_filter_.access();
  }

  virtual FunctorType& functor()
  {
    return functor_.access();
  }

private:
  Common::StorageProvider<FunctorType> functor_;
  Common::ConstStorageProvider<EntityFilterType> entity_filter_;
  Common::ConstStorageProvider<IntersectionFilterType> intersection_filter_;
}; // class EntityAndIntersectionFunctorWrapper


/**
 * \brief To be used within the \sa Walker as internal storage type.
 * \note  Most likely you do not want to use this class directly, but instead append() a Walker to a Walker.
 */
template <class GL>
class WalkerWrapper : public EntityAndIntersectionFunctorWrapper<GL>
{
  using BaseType = EntityAndIntersectionFunctorWrapper<GL>;
  using ThisType = WalkerWrapper<GL>;

public:
  using typename BaseType::EntityFilterType;
  using typename BaseType::IntersectionFilterType;

  /**
   * \attention Takes ownership of entity_filtr_ptr, do not delete manually!
   * \attention Takes ownership of intersection_filtr_ptr, do not delete manually!
   */
  WalkerWrapper(Walker<GL>& walkr,
                const EntityFilterType*&& entity_filtr_ptr,
                const IntersectionFilterType*&& intersection_filtr_ptr)
    : BaseType(walkr,
               new ApplyOn::LambdaFilteredEntities<GL>([&](const auto& grid_layer, const auto& entity) {
                 if (restriction_entity_filter_.access().contains(grid_layer, entity)) {
                   for (const auto& wrapper : walkr.entity_functor_wrappers_) {
                     if (wrapper->filter().contains(grid_layer, entity))
                       return true;
                   }
                   for (const auto& wrapper : walkr.entity_and_intersection_functor_wrappers_) {
                     if (wrapper->entity_filter().contains(grid_layer, entity))
                       return true;
                   }
                   return false;
                 } else
                   return false;
               }),
               new ApplyOn::LambdaFilteredIntersections<GL>([&](const auto& grid_layer, const auto& intersection) {
                 if (restriction_intersection_filter_.access().contains(grid_layer, intersection)) {
                   for (const auto& wrapper : walkr.intersection_functor_wrappers_) {
                     if (wrapper->filter().contains(grid_layer, intersection))
                       return true;
                   }
                   for (const auto& wrapper : walkr.entity_and_intersection_functor_wrappers_) {
                     if (wrapper->intersection_filter().contains(grid_layer, intersection))
                       return true;
                   }
                   return false;
                 } else
                   return false;
               }))
    , restriction_entity_filter_(std::move(entity_filtr_ptr))
    , restriction_intersection_filter_(std::move(intersection_filtr_ptr))
  {
  }

private:
  Common::ConstStorageProvider<EntityFilterType> restriction_entity_filter_;
  Common::ConstStorageProvider<IntersectionFilterType> restriction_intersection_filter_;
}; // class WalkerWrapper


} // namespace internal
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_WALKER_WRAPPER_HH
