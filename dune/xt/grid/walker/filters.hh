// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   Rene Milk       (2014 - 2018)
//   Tobias Leibner  (2015 - 2017)

#ifndef DUNE_XT_GRID_WALKER_FILTERS_HH
#define DUNE_XT_GRID_WALKER_FILTERS_HH

#include <functional>

#include <dune/xt/common/memory.hh>

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace internal {


// forwards
template <class GL>
class CombinedElementFilter;

template <class GL>
class CombinedIntersectionFilter;

template <class GL>
class NegatedElementFilter;

template <class GL>
class NegatedIntersectionFilter;


} // namespace internal


/**
 * \brief Interface for element filters to restrict the range of elements a functor is applied on by the Walker.
 *
 * \note  Filters are provided within the ApplyOn namespace.
 * \note  See also the Walker for typical use cases.
 *
 * \sa Walker
 */
template <class GL>
class ElementFilter
{
  static_assert(is_layer<GL>::value, "");

public:
  using GridViewType = GL;
  using ElementType = extract_entity_t<GridViewType>;

  virtual ~ElementFilter() = default;

  virtual ElementFilter<GridViewType>* copy() const = 0;

  virtual bool contains(const GridViewType& /*grid_layer*/, const ElementType& /*element*/) const = 0;

  std::unique_ptr<ElementFilter<GridViewType>> operator!() const
  {
    return std::make_unique<internal::NegatedElementFilter<GridViewType>>(*this);
  }

  std::unique_ptr<ElementFilter<GridViewType>> operator&&(const ElementFilter<GridViewType>& other) const
  {
    return std::make_unique<internal::CombinedElementFilter<GridViewType>>(
        *this, other, [](const auto& left, const auto& right) { return left && right; });
  }

  std::unique_ptr<ElementFilter<GridViewType>> operator&&(ElementFilter<GridViewType>*&& other) const
  {
    return std::make_unique<internal::CombinedElementFilter<GridViewType>>(
        *this, std::move(other), [](const auto& left, const auto& right) { return left && right; });
  }

  std::unique_ptr<ElementFilter<GridViewType>> operator&&(std::unique_ptr<ElementFilter<GridViewType>> other) const
  {
    return std::make_unique<internal::CombinedElementFilter<GridViewType>>(
        *this, std::move(other), [](const auto& left, const auto& right) { return left && right; });
  }

  std::unique_ptr<ElementFilter<GridViewType>> operator||(const ElementFilter<GridViewType>& other) const
  {
    return std::make_unique<internal::CombinedElementFilter<GridViewType>>(
        *this, other, [](const auto& left, const auto& right) { return left || right; });
  }

  std::unique_ptr<ElementFilter<GridViewType>> operator||(ElementFilter<GridViewType>*&& other) const
  {
    return std::make_unique<internal::CombinedElementFilter<GridViewType>>(
        *this, std::move(other), [](const auto& left, const auto& right) { return left || right; });
  }

  std::unique_ptr<ElementFilter<GridViewType>> operator||(std::unique_ptr<ElementFilter<GridViewType>> other) const
  {
    return std::make_unique<internal::CombinedElementFilter<GridViewType>>(
        *this, std::move(other), [](const auto& left, const auto& right) { return left || right; });
  }
}; // class ElementFilter


/**
 * \brief Interface for intersection filters to restrict the range of intersections a functor is applied on by the
 *        Walker.
 *
 * \note Filters are provided within the ApplyOn namespace.
 * \note See also the Walker for typical use cases.
 *
 * \sa Walker
 */
template <class GL>
class IntersectionFilter
{
  static_assert(is_layer<GL>::value, "");

public:
  using GridViewType = GL;
  using IntersectionType = extract_intersection_t<GridViewType>;

  virtual ~IntersectionFilter<GridViewType>() = default;

  virtual IntersectionFilter<GridViewType>* copy() const = 0;

  virtual bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& /*intersection*/) const = 0;

  std::unique_ptr<IntersectionFilter<GridViewType>> operator!() const
  {
    return std::make_unique<internal::NegatedIntersectionFilter<GridViewType>>(*this);
  }

  std::unique_ptr<IntersectionFilter<GridViewType>> operator&&(const IntersectionFilter<GridViewType>& other) const
  {
    return std::make_unique<internal::CombinedIntersectionFilter<GridViewType>>(
        *this, other, [](const auto& left, const auto& right) { return left && right; });
  }

  std::unique_ptr<IntersectionFilter<GridViewType>> operator&&(IntersectionFilter<GridViewType>*&& other) const
  {
    return std::make_unique<internal::CombinedIntersectionFilter<GridViewType>>(
        *this, std::move(other), [](const auto& left, const auto& right) { return left && right; });
  }

  std::unique_ptr<IntersectionFilter<GridViewType>>
  operator&&(std::unique_ptr<IntersectionFilter<GridViewType>> other) const
  {
    return std::make_unique<internal::CombinedIntersectionFilter<GridViewType>>(
        *this, std::move(other), [](const auto& left, const auto& right) { return left && right; });
  }

  std::unique_ptr<IntersectionFilter<GridViewType>> operator||(const IntersectionFilter<GridViewType>& other) const
  {
    return std::make_unique<internal::CombinedIntersectionFilter<GridViewType>>(
        *this, other, [](const auto& left, const auto& right) { return left || right; });
  }

  std::unique_ptr<IntersectionFilter<GridViewType>> operator||(IntersectionFilter<GridViewType>*&& other) const
  {
    return std::make_unique<internal::CombinedIntersectionFilter<GridViewType>>(
        *this, std::move(other), [](const auto& left, const auto& right) { return left || right; });
  }

  std::unique_ptr<IntersectionFilter<GridViewType>>
  operator||(std::unique_ptr<IntersectionFilter<GridViewType>> other) const
  {
    return std::make_unique<internal::CombinedIntersectionFilter<GridViewType>>(
        *this, std::move(other), [](const auto& left, const auto& right) { return left || right; });
  }
}; // class IntersectionFilter


namespace ApplyOn {


/**
 * \brief A filter which selects all elements.
 *
 * \sa Walker
 * \sa ElementFilter
 */
template <class GL>
class AllElements : public ElementFilter<GL>
{
  using BaseType = ElementFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;

  explicit AllElements() = default;

  ElementFilter<GridViewType>* copy() const override final
  {
    return new AllElements<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const ElementType& /*element*/) const override final
  {
    return true;
  }
}; // class AllElements


/**
 * \brief A filter which selects no elements.
 *
 * \sa Walker
 * \sa ElementFilter
 */
template <class GL>
class NoElements : public ElementFilter<GL>
{
  using BaseType = ElementFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;

  explicit NoElements() = default;

  ElementFilter<GridViewType>* copy() const override final
  {
    return new NoElements<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const ElementType& /*element*/) const override final
  {
    return false;
  }
}; // class NoElements


/**
 * \brief A filter which selects entities which have a boundary intersection.
 *
 * \sa Walker
 * \sa ElementFilter
 */
template <class GL>
class BoundaryEntities : public ElementFilter<GL>
{
  using BaseType = ElementFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;

  explicit BoundaryEntities() = default;

  ElementFilter<GridViewType>* copy() const override final
  {
    return new BoundaryEntities<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const ElementType& element) const override final
  {
    return element.hasBoundaryIntersections();
  }
}; // class BoundaryEntities


/**
 * \brief A filter which selects entities based on a lambda expression.
 *
 * \sa Walker
 * \sa ElementFilter
 */
template <class GL>
class LambdaFilteredElements : public ElementFilter<GL>
{
  using BaseType = ElementFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;
  using LambdaType = std::function<bool(const GridViewType&, const ElementType&)>;

  explicit LambdaFilteredElements(LambdaType lambda)
    : lambda_(lambda)
  {
  }

  ElementFilter<GridViewType>* copy() const override final
  {
    return new LambdaFilteredElements<GridViewType>(lambda_);
  }

  bool contains(const GridViewType& grid_layer, const ElementType& element) const override final
  {
    return lambda_(grid_layer, element);
  }

private:
  const LambdaType lambda_;
}; // class LambdaFilteredElements


/**
 * \brief A filter which selects entities in a compatible PartitionSet.
 *
 * \sa Walker
 * \sa ElementFilter
 */
template <class GL, class PartitionSetType>
class PartitionSetEntities : public ElementFilter<GL>
{
  using BaseType = ElementFilter<GL>;

public:
  using typename BaseType::GridViewType;
  typedef typename BaseType::ElementType ElementType;

  explicit PartitionSetEntities() = default;

  ElementFilter<GridViewType>* copy() const override final
  {
    return new PartitionSetEntities<GridViewType, PartitionSetType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const ElementType& element) const override final
  {
    return PartitionSetType::contains(element.partitionType());
  }
}; // class PartitionSetEntities


/**
 * \brief A filter which selects all intersections.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GL>
class AllIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit AllIntersections() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new AllIntersections<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& /*intersection*/) const override final
  {
    return true;
  }
}; // class AllIntersections


/**
 * \brief A filter which selects no intersections.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GL>
class NoIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit NoIntersections() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new NoIntersections<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& /*intersection*/) const override final
  {
    return false;
  }
}; // class NoIntersections


/**
 * \brief A filter which selects each inner intersection.
 *
 * \note To decide if this in an inner intersection,
\code
intersection.neighbor() && !intersection.boundary()
\endcode
 *       is used.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GL>
class InnerIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit InnerIntersections() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new InnerIntersections<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return intersection.neighbor() && !intersection.boundary();
  }
}; // class InnerIntersections


/**
 * \brief A filter which selects each inner intersection only once.
 *
 * \note To decide if this in an inner intersection,
\code
intersection.neighbor() && !intersection.boundary()
\endcode
 *       is used, and true is returned, if the index of the inside() element is smaller than the index of the outside()
 *       element.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GL>
class InnerIntersectionsOnce : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit InnerIntersectionsOnce() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new InnerIntersectionsOnce<GridViewType>();
  }

  bool contains(const GridViewType& grid_layer, const IntersectionType& intersection) const override final
  {
    if (intersection.neighbor() && !intersection.boundary()) {
      const auto inside_element = intersection.inside();
      const auto outside_element = intersection.outside();
      return grid_layer.indexSet().index(inside_element) < grid_layer.indexSet().index(outside_element);
    } else
      return false;
  }
}; // class InnerIntersectionsOnce


/**
 * \brief Selects each inner intersection in given partition only once.
 *
 * \sa InnerIntersectionsOnce
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GL, class PartitionSetType>
class PartitionSetInnerIntersectionsOnce : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit PartitionSetInnerIntersectionsOnce() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new PartitionSetInnerIntersectionsOnce<GridViewType, PartitionSetType>();
  }

  bool apply_on(const GridViewType& grid_layer, const IntersectionType& intersection) const override final
  {
    if (intersection.neighbor() && !intersection.boundary()
        && PartitionSetType::contains(intersection.inside().partitionType())) {
      const auto inside_element = intersection.inside();
      const auto outside_element = intersection.outside();
      if (!PartitionSetType::contains(intersection.outside().partitionType()))
        return true;
      return grid_layer.indexSet().index(inside_element) < grid_layer.indexSet().index(outside_element);
    } else
      return false;
  }
}; // class PartitionSetInnerIntersectionsOnce


/**
 * \brief A filter which selects boundary intersections.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GL>
class BoundaryIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit BoundaryIntersections() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new BoundaryIntersections<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return intersection.boundary();
  }
}; // class BoundaryIntersections


/**
 * \brief A filter which selects intersections on a non-periodic boundary.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GL>
class NonPeriodicBoundaryIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit NonPeriodicBoundaryIntersections() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new NonPeriodicBoundaryIntersections<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return intersection.boundary() && !intersection.neighbor();
  }
}; // class BoundaryIntersections


/**
 * \brief A filter which selects intersections on a periodic boundary only once.
 *
 *        To decide if this a periodic intersection,
\code
intersection.neighbor() && intersection.boundary()
\endcode
 *        is used.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GL>
class PeriodicBoundaryIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit PeriodicBoundaryIntersections() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new PeriodicBoundaryIntersections<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return intersection.neighbor() && intersection.boundary();
  }
}; // class PeriodicBoundaryIntersections


/**
 * \brief A filter which selects intersections on a periodic boundary only once.
 *
 *        To decide if this in an periodic intersection,
\code
intersection.neighbor() && intersection.boundary()
\endcode
 *        is used, and true is returned, if the index of the inside() element is smaller than the index of the outside()
 *        element.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GL>
class PeriodicBoundaryIntersectionsOnce : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit PeriodicBoundaryIntersectionsOnce() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new PeriodicBoundaryIntersectionsOnce<GridViewType>();
  }

  bool contains(const GridViewType& grid_layer, const IntersectionType& intersection) const override final
  {
    if (intersection.neighbor() && intersection.boundary()) {
      const auto inside_element = intersection.inside();
      const auto outside_element = intersection.outside();
      return grid_layer.indexSet().index(inside_element) < grid_layer.indexSet().index(outside_element);
    } else {
      return false;
    }
  }
}; // class PeriodicBoundaryIntersectionsOnce


/**
 * \brief A filter which selects intersections based on a lambda expression.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GL>
class LambdaFilteredIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;
  using LambdaType = std::function<bool(const GridViewType&, const IntersectionType&)>;

  explicit LambdaFilteredIntersections(LambdaType lambda)
    : lambda_(lambda)
  {
  }

  virtual IntersectionFilter<GridViewType>* copy() const override final
  {
    return new LambdaFilteredIntersections<GridViewType>(lambda_);
  }

  bool contains(const GridViewType& grid_layer, const IntersectionType& intersection) const override final
  {
    return lambda_(grid_layer, intersection);
  }

private:
  const LambdaType lambda_;
}; // class BoundaryIntersections


/**
 * \brief A filter which selects intersections on a given part of the boundary.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GL>
class CustomBoundaryIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  /**
   * \attention Takes ownership of boundary_type, do not delete manually!
   */
  explicit CustomBoundaryIntersections(const BoundaryInfo<IntersectionType>& boundary_info,
                                       BoundaryType*&& boundary_type)
    : boundary_info_(boundary_info)
    , boundary_type_(boundary_type)
  {
  }

  explicit CustomBoundaryIntersections(const BoundaryInfo<IntersectionType>& boundary_info,
                                       const std::shared_ptr<BoundaryType>& boundary_type)
    : boundary_info_(boundary_info)
    , boundary_type_(boundary_type)
  {
  }

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new CustomBoundaryIntersections<GridViewType>(boundary_info_, boundary_type_);
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return boundary_info_.type(intersection) == *boundary_type_;
  }

protected:
  const BoundaryInfo<IntersectionType>& boundary_info_;
  const std::shared_ptr<BoundaryType> boundary_type_;
}; // class CustomBoundaryIntersections


/**
 * \brief A filter which selects intersections on a given part of the physical boundary and the process boundary.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GL>
class CustomBoundaryAndProcessIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  /**
   * \attention Takes ownership of boundary_type, do not delete manually!
   */
  explicit CustomBoundaryAndProcessIntersections(const BoundaryInfo<IntersectionType>& boundary_info,
                                                 BoundaryType*&& boundary_type)
    : boundary_info_(boundary_info)
    , boundary_type_(std::move(boundary_type))
  {
  }

  explicit CustomBoundaryAndProcessIntersections(const BoundaryInfo<IntersectionType>& boundary_info,
                                                 const std::shared_ptr<BoundaryType>& boundary_type)
    : boundary_info_(boundary_info)
    , boundary_type_(boundary_type)
  {
  }

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new CustomBoundaryAndProcessIntersections<GridViewType>(boundary_info_, boundary_type_);
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    const bool process_boundary = !intersection.neighbor() && !intersection.boundary();
    const bool physical_boundary = boundary_info_.type(intersection) == *boundary_type_;
    return process_boundary || physical_boundary;
  }

protected:
  const BoundaryInfo<IntersectionType>& boundary_info_;
  const std::shared_ptr<BoundaryType> boundary_type_;
}; // class CustomBoundaryAndProcessIntersections


/**
 * \brief A filter which selects intersections on a the process boundary.
 *
 * \sa Walker
 * \sa IntersectionFilter
 */
template <class GL>
class ProcessIntersections : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit ProcessIntersections() = default;

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new ProcessIntersections<GridViewType>();
  }

  bool contains(const GridViewType& /*grid_layer*/, const IntersectionType& intersection) const override final
  {
    return (!intersection.neighbor() && !intersection.boundary());
  }
}; // class ProcessIntersections


} // namespace ApplyOn
namespace internal {


template <class GL>
class CombinedIntersectionFilter : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  CombinedIntersectionFilter(const BaseType& left,
                             const BaseType& right,
                             std::function<bool(const bool&, const bool&)> combine_lambda)
    : left_(left.copy())
    , right_(right.copy())
    , combine_lambda_(combine_lambda)
  {
  }

  CombinedIntersectionFilter(const BaseType& left,
                             BaseType*&& right,
                             std::function<bool(const bool&, const bool&)> combine_lambda)
    : left_(left.copy())
    , right_(std::move(right))
    , combine_lambda_(combine_lambda)
  {
  }

  CombinedIntersectionFilter(const BaseType& left,
                             std::unique_ptr<BaseType> right,
                             std::function<bool(const bool&, const bool&)> combine_lambda)
    : left_(left.copy())
    , right_(std::move(right))
    , combine_lambda_(combine_lambda)
  {
  }

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new CombinedIntersectionFilter<GridViewType>(*left_, *right_, combine_lambda_);
  }

  bool contains(const GridViewType& grid_layer, const IntersectionType& intersection) const override final
  {
    return combine_lambda_(left_->contains(grid_layer, intersection), right_->contains(grid_layer, intersection));
  }

private:
  const std::unique_ptr<BaseType> left_;
  const std::unique_ptr<BaseType> right_;
  const std::function<bool(const bool&, const bool&)> combine_lambda_;
}; // class CombinedIntersectionFilter


template <class GL>
class CombinedElementFilter : public ElementFilter<GL>
{
  using BaseType = ElementFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;

  CombinedElementFilter(const BaseType& left,
                        const BaseType& right,
                        std::function<bool(const bool&, const bool&)> combine_lambda)
    : left_(left.copy())
    , right_(right.copy())
    , combine_lambda_(combine_lambda)
  {
  }

  CombinedElementFilter(const BaseType& left,
                        BaseType*&& right,
                        std::function<bool(const bool&, const bool&)> combine_lambda)
    : left_(left.copy())
    , right_(std::move(right))
    , combine_lambda_(combine_lambda)
  {
  }

  CombinedElementFilter(const BaseType& left,
                        std::unique_ptr<BaseType> right,
                        std::function<bool(const bool&, const bool&)> combine_lambda)
    : left_(left.copy())
    , right_(std::move(right))
    , combine_lambda_(combine_lambda)
  {
  }

  ElementFilter<GridViewType>* copy() const override final
  {
    return new CombinedElementFilter<GridViewType>(*left_, *right_, combine_lambda_);
  }

  bool contains(const GridViewType& grid_layer, const ElementType& element) const override final
  {
    return combine_lambda_(left_->contains(grid_layer, element), right_->contains(grid_layer, element));
  }

private:
  const std::unique_ptr<BaseType> left_;
  const std::unique_ptr<BaseType> right_;
  const std::function<bool(const bool&, const bool&)> combine_lambda_;
}; // class CombinedElementFilter


template <class GL>
class NegatedIntersectionFilter : public IntersectionFilter<GL>
{
  using BaseType = IntersectionFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  NegatedIntersectionFilter(const BaseType& filter)
    : filter_(filter.copy())
  {
  }

  NegatedIntersectionFilter(BaseType*&& filter)
    : filter_(std::move(filter))
  {
  }

  IntersectionFilter<GridViewType>* copy() const override final
  {
    return new NegatedIntersectionFilter<GridViewType>(*filter_);
  }

  bool contains(const GridViewType& grid_layer, const IntersectionType& intersection) const override final
  {
    return !filter_->contains(grid_layer, intersection);
  }

private:
  const std::unique_ptr<BaseType> filter_;
}; // class NegatedIntersectionFilter


template <class GL>
class NegatedElementFilter : public ElementFilter<GL>
{
  using BaseType = ElementFilter<GL>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;

  NegatedElementFilter(const BaseType& filter)
    : filter_(filter.copy())
  {
  }

  NegatedElementFilter(BaseType*&& filter)
    : filter_(std::move(filter))
  {
  }

  ElementFilter<GridViewType>* copy() const override final
  {
    return new NegatedElementFilter<GridViewType>(*filter_);
  }

  bool contains(const GridViewType& grid_layer, const ElementType& element) const override final
  {
    return !filter_->contains(grid_layer, element);
  }

private:
  const std::unique_ptr<BaseType> filter_;
}; // class NegatedElementFilter


} // namespace internal
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_WALKER_FILTERS_HH
