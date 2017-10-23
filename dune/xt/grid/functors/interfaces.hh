// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   Rene Milk       (2014 - 2018)

#ifndef DUNE_XT_GRID_FUNCTORS_INTERFACES_HH
#define DUNE_XT_GRID_FUNCTORS_INTERFACES_HH

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/grid/type_traits.hh>

namespace Dune {
namespace XT {
namespace Grid {


/**
 * \brief Interface for functors which are applied to (codim 0) entities of a grid layer by the \sa GridWalker.
 * \sa    EntityReturnFunctor
 */
template <class GL>
class EntityFunctor
{
  static_assert(is_layer<GL>::value, "");

public:
  using GridLayerType = GL;
  using EntityType = extract_entity_t<GridLayerType>;

  virtual ~EntityFunctor() = default;

  virtual void prepare()
  {
  }

  virtual void apply_local(const EntityType& entity) = 0;

  virtual void finalize()
  {
  }
}; // class EntityFunctor


/**
 * \brief Interface for functors which are like \sa EntityFunctor but return results.
 *
 *        Think of computing the L^2 norm of a function f:
\code
const auto local_l2_norm = functor.compute_locally(entity);
\endcode
          would compute and return \int_entity f^2 \dx. On the other hand, apply_local(entity) would be implemented
          to call compute_locally and store the result, i.e. an implementation might look like:
\code
void apply_local(const EntityType& entity) override final
{
  this->result_ += functor.compute_locally(entity);
}
\endcode
          In that case, result_ would be a private member which is set to zero in prepare(). An implementation of
          finalize() might then look like:
\code
void finalize() override final
{
  this->result_ = std::sqrt(this->result);
}
\endcode
          The purpose of result() is then simply to return result_;
 */
template <class GL, class R>
class EntityReturnFunctor : public EntityFunctor<GL>
{
  using BaseType = EntityFunctor<GL>;

public:
  using ResultType = R;
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;

  virtual ~EntityReturnFunctor() = default;

  /**
   * \note Calling this method should not alter the state of this functor, i.e. calling
\code
auto result = func.compute_locally(entity);
func.apply_local(entity);
\endcode
   *       should give the same as:
\code
func.apply_local(entity);
\endcode
   */
  virtual ResultType compute_locally(const EntityType& entity) = 0;

  /**
   * \brief Can be called to access the final result of the functor.
   */
  virtual ResultType result() const = 0;
}; // class EntityReturnFunctor


/**
 * \brief Interface for functors which are applied to (codim 1) intersection of a grid layer by the \sa GridWalker.
 * \sa    EntityFunctor
 * \sa    EntityReturnFunctor
 */
template <class GL>
class IntersectionFunctor
{
  static_assert(is_layer<GL>::value, "");

public:
  using GridLayerType = GL;
  using EntityType = extract_entity_t<GridLayerType>;
  using IntersectionType = extract_intersection_t<GridLayerType>;

  virtual ~IntersectionFunctor() = default;

  virtual void prepare()
  {
  }

  /**
   * \note The meaning of outside_intersection depends on the circumstances. In general, the result of
   *       intersection.outside() is given, but this might differ on periodic or boundary intersections.
   */
  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& inside_entity,
                           const EntityType& outside_entity) = 0;

  virtual void finalize()
  {
  }
}; // class IntersectionFunctor


/**
 * \brief Interface for functors which are like \sa IntersectionFunctor but return results.
 * \sa    EntityFunctor
 * \sa    EntityReturnFunctor
 * \sa    IntersectionFunctor
 */
template <class GL, class R>
class IntersectionReturnFunctor : public IntersectionFunctor<GL>
{
  using BaseType = IntersectionFunctor<GL>;

public:
  using ResultType = R;
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;

  virtual ~IntersectionReturnFunctor() = default;

  /**
   * \note Calling this method should not alter the state of this functor, i.e. calling
\code
auto result = func.compute_locally(intersection, inside_entity, outside_entity);
func.apply_local(intersection, inside_entity, outside_entity);
\endcode
   *       should give the same as:
\code
func.apply_local(intersection, inside_entity, outside_entity);
\endcode
   */
  virtual ResultType compute_locally(const IntersectionType& intersection,
                                     const EntityType& inside_entity,
                                     const EntityType& outside_entity) = 0;

  virtual ResultType result() const = 0;
}; // class IntersectionReturnFunctor


/**
 * \brief Interface for functors which are applied to entities and intersections of a grid layer by the \sa GridWalker.
 * \sa    EntityFunctor
 * \sa    IntersectionFunctor
 * \sa    EntityReturnFunctor
 * \sa    IntersectionReturnFunctor
 */
template <class GL>
class EntityAndIntersectionFunctor
{
  static_assert(is_layer<GL>::value, "");

public:
  using GridLayerType = GL;
  using EntityType = extract_entity_t<GridLayerType>;
  using IntersectionType = extract_intersection_t<GridLayerType>;

  virtual ~EntityAndIntersectionFunctor() = default;

  virtual void prepare()
  {
  }

  virtual void apply_local(const EntityType& entity) = 0;

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& inside_entity,
                           const EntityType& outside_entity) = 0;

  virtual void finalize()
  {
  }
}; // class EntityAndIntersectionFunctor


/**
 * \brief Interface for functors which are like \sa EntityAndIntersectionFunctor but return results.
 * \sa    EntityFunctor
 * \sa    EntityReturnFunctor
 * \sa    IntersectionFunctor
 * \sa    IntersectionReturnFunctor
 * \sa    EntityAndIntersectionFunctor
 */
template <class GL, class R>
class EntityAndIntersectionReturnFunctor : public EntityAndIntersectionFunctor<GL>
{
  using BaseType = EntityAndIntersectionFunctor<GL>;

public:
  using ResultType = R;
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;

  virtual ~EntityAndIntersectionReturnFunctor() = default;

  /**
   * \note Calling this method should not alter the state of this functor, i.e. calling
\code
auto result = func.compute_locally(entity);
func.apply_local(entity);
\endcode
   *       should give the same as:
\code
func.apply_local(entity);
\endcode
   */
  virtual ResultType compute_locally(const EntityType& entity) = 0;

  /**
   * \note Calling this method should not alter the state of this functor, i.e. calling
\code
auto result = func.compute_locally(intersection, inside_entity, outside_entity);
func.apply_local(intersection, inside_entity, outside_entity);
\endcode
   *       should give the same as:
\code
func.apply_local(intersection, inside_entity, outside_entity);
\endcode
   */
  virtual ResultType compute_locally(const IntersectionType& intersection,
                                     const EntityType& inside_entity,
                                     const EntityType& outside_entity) = 0;

  virtual ResultType result() const = 0;
}; // class EntityAndIntersectionReturnFunctor


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_FUNCTORS_INTERFACES_HH
