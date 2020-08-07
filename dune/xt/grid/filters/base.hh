// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2019 - 2020)

#ifndef DUNE_XT_GRID_FILTERS_BASE_HH
#define DUNE_XT_GRID_FILTERS_BASE_HH

#include <functional>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune::XT::Grid {
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
  static_assert(is_layer<GL>::value);

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
  static_assert(is_layer<GL>::value);

public:
  using GridViewType = GL;
  using IntersectionType = extract_intersection_t<GridViewType>;

  mutable Common::DefaultLogger logger;

  IntersectionFilter(const std::string& logging_prefix = "xt.grid.intersectionfilter",
                     const bool logging_disabled = true)
    : logger(logging_prefix, logging_disabled)
  {}

  IntersectionFilter(const IntersectionFilter<GL>&) = default;

  IntersectionFilter(IntersectionFilter<GL>&&) = default;

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
  {}

  CombinedIntersectionFilter(const BaseType& left,
                             BaseType*&& right,
                             std::function<bool(const bool&, const bool&)> combine_lambda)
    : left_(left.copy())
    , right_(std::move(right))
    , combine_lambda_(combine_lambda)
  {}

  CombinedIntersectionFilter(const BaseType& left,
                             std::unique_ptr<BaseType> right,
                             std::function<bool(const bool&, const bool&)> combine_lambda)
    : left_(left.copy())
    , right_(std::move(right))
    , combine_lambda_(combine_lambda)
  {}

  IntersectionFilter<GridViewType>* copy() const final
  {
    return new CombinedIntersectionFilter<GridViewType>(*left_, *right_, combine_lambda_);
  }

  bool contains(const GridViewType& grid_layer, const IntersectionType& intersection) const final
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
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;

  CombinedElementFilter(const BaseType& left,
                        const BaseType& right,
                        std::function<bool(const bool&, const bool&)> combine_lambda)
    : left_(left.copy())
    , right_(right.copy())
    , combine_lambda_(combine_lambda)
  {}

  CombinedElementFilter(const BaseType& left,
                        BaseType*&& right,
                        std::function<bool(const bool&, const bool&)> combine_lambda)
    : left_(left.copy())
    , right_(std::move(right))
    , combine_lambda_(combine_lambda)
  {}

  CombinedElementFilter(const BaseType& left,
                        std::unique_ptr<BaseType> right,
                        std::function<bool(const bool&, const bool&)> combine_lambda)
    : left_(left.copy())
    , right_(std::move(right))
    , combine_lambda_(combine_lambda)
  {}

  ElementFilter<GridViewType>* copy() const final
  {
    return new CombinedElementFilter<GridViewType>(*left_, *right_, combine_lambda_);
  }

  bool contains(const GridViewType& grid_layer, const ElementType& element) const final
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
  {}

  NegatedIntersectionFilter(BaseType*&& filter)
    : filter_(std::move(filter))
  {}

  IntersectionFilter<GridViewType>* copy() const final
  {
    return new NegatedIntersectionFilter<GridViewType>(*filter_);
  }

  bool contains(const GridViewType& grid_layer, const IntersectionType& intersection) const final
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
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;

  NegatedElementFilter(const BaseType& filter)
    : filter_(filter.copy())
  {}

  NegatedElementFilter(BaseType*&& filter)
    : filter_(std::move(filter))
  {}

  ElementFilter<GridViewType>* copy() const final
  {
    return new NegatedElementFilter<GridViewType>(*filter_);
  }

  bool contains(const GridViewType& grid_layer, const ElementType& element) const final
  {
    return !filter_->contains(grid_layer, element);
  }

private:
  const std::unique_ptr<BaseType> filter_;
}; // class NegatedElementFilter


} // namespace internal
} // namespace Dune::XT::Grid

#endif // DUNE_XT_GRID_FILTERS_BASE_HH
