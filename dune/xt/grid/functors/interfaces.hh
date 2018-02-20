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
 * \brief Interface for functors which are applied to elements (codim 0 entities) of a grid layer by the \sa GridWalker.
 * \sa    ElementReturnFunctor
 */
template <class GL>
class ElementFunctor
{
  static_assert(is_layer<GL>::value, "");

public:
  using GridViewType = GL;
  using ElementType = extract_entity_t<GridViewType>;

  virtual ~ElementFunctor() = default;

  virtual void prepare()
  {
  }

  virtual void apply_local(const ElementType& element) = 0;

  virtual void finalize()
  {
  }
}; // class ElementFunctor


/**
 * \brief Interface for functors which are like \sa ElementFunctor but return results.
 *
 *        Think of computing the L^2 norm of a function f:
\code
const auto local_l2_norm = functor.compute_locally(entity);
\endcode
          would compute and return \int_entity f^2 \dx. On the other hand, apply_local(entity) would be implemented
          to call compute_locally and store the result, i.e. an implementation might look like:
\code
void apply_local(const ElementType& element) override final
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
class ElementReturnFunctor : public ElementFunctor<GL>
{
  using BaseType = ElementFunctor<GL>;

public:
  using ResultType = R;
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;

  virtual ~ElementReturnFunctor() = default;

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
  virtual ResultType compute_locally(const ElementType& element) = 0;

  /**
   * \brief Can be called to access the final result of the functor.
   */
  virtual ResultType result() const = 0;
}; // class ElementReturnFunctor


/**
 * \brief Interface for functors which are applied to (codim 1) intersection of a grid layer by the \sa GridWalker.
 * \sa    ElementFunctor
 * \sa    ElementReturnFunctor
 */
template <class GL>
class IntersectionFunctor
{
  static_assert(is_layer<GL>::value, "");

public:
  using GridViewType = GL;
  using ElementType = extract_entity_t<GridViewType>;
  using IntersectionType = extract_intersection_t<GridViewType>;

  virtual ~IntersectionFunctor() = default;

  virtual void prepare()
  {
  }

  /**
   * \note The meaning of outside_intersection depends on the circumstances. In general, the result of
   *       intersection.outside() is given, but this might differ on periodic or boundary intersections.
   */
  virtual void apply_local(const IntersectionType& intersection,
                           const ElementType& inside_element,
                           const ElementType& outside_element) = 0;

  virtual void finalize()
  {
  }
}; // class IntersectionFunctor


/**
 * \brief Interface for functors which are like \sa IntersectionFunctor but return results.
 * \sa    ElementFunctor
 * \sa    ElementReturnFunctor
 * \sa    IntersectionFunctor
 */
template <class GL, class R>
class IntersectionReturnFunctor : public IntersectionFunctor<GL>
{
  using BaseType = IntersectionFunctor<GL>;

public:
  using ResultType = R;
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;
  using typename BaseType::IntersectionType;

  virtual ~IntersectionReturnFunctor() = default;

  /**
   * \note Calling this method should not alter the state of this functor, i.e. calling
\code
auto result = func.compute_locally(intersection, inside_element, outside_element);
func.apply_local(intersection, inside_element, outside_element);
\endcode
   *       should give the same as:
\code
func.apply_local(intersection, inside_element, outside_element);
\endcode
   */
  virtual ResultType compute_locally(const IntersectionType& intersection,
                                     const ElementType& inside_element,
                                     const ElementType& outside_element) = 0;

  virtual ResultType result() const = 0;
}; // class IntersectionReturnFunctor


/**
 * \brief Interface for functors which are applied to entities and intersections of a grid layer by the \sa GridWalker.
 * \sa    ElementFunctor
 * \sa    IntersectionFunctor
 * \sa    ElementReturnFunctor
 * \sa    IntersectionReturnFunctor
 */
template <class GL>
class ElementAndIntersectionFunctor
{
  static_assert(is_layer<GL>::value, "");

public:
  using GridViewType = GL;
  using ElementType = extract_entity_t<GridViewType>;
  using IntersectionType = extract_intersection_t<GridViewType>;

  virtual ~ElementAndIntersectionFunctor() = default;

  virtual void prepare()
  {
  }

  virtual void apply_local(const ElementType& element) = 0;

  virtual void apply_local(const IntersectionType& intersection,
                           const ElementType& inside_element,
                           const ElementType& outside_element) = 0;

  virtual void finalize()
  {
  }
}; // class ElementAndIntersectionFunctor


/**
 * \brief Interface for functors which are like \sa ElementAndIntersectionFunctor but return results.
 * \sa    ElementFunctor
 * \sa    ElementReturnFunctor
 * \sa    IntersectionFunctor
 * \sa    IntersectionReturnFunctor
 * \sa    ElementAndIntersectionFunctor
 */
template <class GL, class R>
class ElementAndIntersectionReturnFunctor : public ElementAndIntersectionFunctor<GL>
{
  using BaseType = ElementAndIntersectionFunctor<GL>;

public:
  using ResultType = R;
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;
  using typename BaseType::IntersectionType;

  virtual ~ElementAndIntersectionReturnFunctor() = default;

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
  virtual ResultType compute_locally(const ElementType& element) = 0;

  /**
   * \note Calling this method should not alter the state of this functor, i.e. calling
\code
auto result = func.compute_locally(intersection, inside_element, outside_element);
func.apply_local(intersection, inside_element, outside_element);
\endcode
   *       should give the same as:
\code
func.apply_local(intersection, inside_element, outside_element);
\endcode
   */
  virtual ResultType compute_locally(const IntersectionType& intersection,
                                     const ElementType& inside_element,
                                     const ElementType& outside_element) = 0;

  virtual ResultType result() const = 0;
}; // class ElementAndIntersectionReturnFunctor


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_FUNCTORS_INTERFACES_HH
