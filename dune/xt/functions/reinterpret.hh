// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2016, 2018)
//   Tobias Leibner  (2016 - 2018)

#ifndef DUNE_XT_FUNCTION_REINTERPRET_HH
#define DUNE_XT_FUNCTION_REINTERPRET_HH
#if 0
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#if DUNE_VERSION_NEWER(DUNE_COMMON, 3, 9) // EXADUNE
#include <dune/geometry/referenceelements.hh>
#else
#include <dune/geometry/referenceelements.hh>
#endif

#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/search.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/xt/functions/exceptions.hh>
#include <dune/xt/functions/interfaces.hh>

namespace Dune {
namespace XT {
namespace Functions {


/**
 * \brief Allows to reinterpret a given LocalizableFunctionInterface, associated with a given grid layer, on a different
 *        grid layer.
 *
 *        Therefore, we search for the correct element in the original grid layer and use the corresponding
 *        local_function to provide an evaluation for a point on the new grid layer. The physical domain covered by the
 *        new grid layer should thus be contained in the physical domain of the original grid layer. This is mainly used
 *        in the context of prolongations.
 *
 * \note  The current implementation is not thread safe (due to the element search).
 *
 * \note  There is no way to reliably obtain the local polynomial order of the source, and we thus use the order of the
 *        local_function corresponding to the first element.
 */
template <class SourceType, class GridLayerType>
class ReinterpretFunction : public LocalizableFunctionInterface<XT::Grid::extract_element_t<GridLayerType>,
                                                                SourceType::dimRange,
                                                                SourceType::dimRangeCols,
                                                                typename SourceType::RangeFieldType>
{
  static_assert(is_localizable_function<SourceType>::value, "");
  static_assert(Grid::is_layer<GridLayerType>::value, "");
  typedef LocalizableFunctionInterface<XT::Grid::extract_element_t<GridLayerType>,
                                       SourceType::dimRange,
                                       SourceType::dimRangeCols,
                                       typename SourceType::RangeFieldType>
      BaseType;
  typedef ReinterpretFunction<SourceType, GridLayerType> ThisType;

public:
  using typename BaseType::ElementType;
  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;
  using typename BaseType::RangeFieldType;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using typename BaseType::LocalFunctionType;

private:
  class ReinterpretLocalfunction
      : public LocalFunctionInterface<ElementType, dimRange, dimRangeCols, RangeFieldType>
  {
    typedef LocalFunctionInterface<ElementType, dimRange, dimRangeCols, RangeFieldType>
        BaseType;

  public:
    using typename BaseType::DomainType;
    using typename BaseType::RangeType;
    using typename BaseType::DerivativeRangeType;

    ReinterpretLocalfunction(const ElementType& element, const size_t order, const ThisType& func)
      : BaseType(element)
      , order_(order)
      , func_(func)
      , points_(1)
    {
    }

    size_t order(const Common::Parameter& /*mu*/ = {}) const override final
    {
      return order_;
    }

    RangeType evaluate(const DomainType& xx, const Common::Parameter& mu = {}) const override final
    {
      RangeType ret;
      points_[0] = this->element().geometry().global(xx);
      const auto source_element_ptr_unique_ptrs = func_.element_search_(points_);
      if (source_element_ptr_unique_ptrs.size() != 1)
        DUNE_THROW(Exceptions::reinterpretation_error,
                   "It was not possible to find a source element for this point:\n\n"
                       << points_[0]);
      const auto& source_element = *source_element_ptr_unique_ptrs[0];
      const auto source_local_function = func_.source_.local_function(source_element);
      source_local_function->evaluate(source_element.geometry().local(points_[0]), ret, mu);
    } // ... evaluate(...)

    DerivativeRangeType jacobian(const DomainType& xx, const Common::Parameter& mu = {}) const override final
    {
      points_[0] = this->element().geometry().global(xx);
      const auto source_element_ptr_unique_ptrs = func_.element_search_(points_);
      if (source_element_ptr_unique_ptrs.size() != 1)
        DUNE_THROW(Exceptions::reinterpretation_error,
                   "It was not possible to find a source element for this point:\n\n"
                       << points_[0]);
      const auto& source_element = *source_element_ptr_unique_ptrs[0];
      const auto source_local_function = func_.source_.local_function(source_element);
      source_local_function->jacobian(source_element.geometry().local(points_[0]), ret, mu);
    } // ... jacobian(...)

  private:
    const size_t order_;
    const ThisType& func_;
    mutable std::vector<DomainType> points_;
  }; // class ReinterpretLocalfunction

public:
  static std::string static_id()
  {
    return BaseType::static_id() + ".reinterpret";
  }

  ReinterpretFunction(const SourceType& source, const GridLayerType& source_grid_layer)
    : source_(source)
    , source_grid_layer_(source_grid_layer)
    , element_search_(source_grid_layer_)
    , guessed_source_order_(source_.local_function(*source_grid_layer_.template begin<0>())->order())
  {
  }

  virtual ~ReinterpretFunction() = default;

  virtual std::unique_ptr<LocalFunctionType> local_function(const ElementType& element) const
  {
    return Common::make_unique<ReinterpretLocalfunction>(element, guessed_source_order_, *this);
  }

  virtual std::string type() const
  {
    return static_id() + "_of_" + source_.type();
  }

  virtual std::string name() const
  {
    return "reinterpretation of " + source_.name();
  }

private:
  friend class ReinterpretLocalfunction;

  const SourceType& source_;
  const GridLayerType& source_grid_layer_;
  mutable XT::Grid::EntityInlevelSearch<GridLayerType> element_search_;
  const size_t guessed_source_order_;
}; // class ReinterpretFunction


} // namespace Functions
} // namespace XT
} // namespace Dune
#endif

#endif // DUNE_XT_FUNCTION_REINTERPRET_HH
