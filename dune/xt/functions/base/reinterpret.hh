// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2016, 2018)
//   TiKeil          (2018)
//   Tobias Leibner  (2016 - 2018)

#ifndef DUNE_XT_FUNCTION_BASE_REINTERPRET_HH
#define DUNE_XT_FUNCTION_BASE_REINTERPRET_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/xt/grid/search.hh>

#include <dune/xt/functions/exceptions.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

namespace Dune {
namespace XT {
namespace Functions {


/**
 * \brief Reinterprets a given localizable function (associated with a given grid layer), on a different grid layer.
 *
 *        Therefore, we search for the correct element in the original grid layer and use the corresponding
 *        local_function to provide an evaluation for a point on the new grid layer. The physical domain covered by the
 *        new grid layer should thus be contained in the physical domain of the original grid layer. This is mainly used
 *        in the context of prolongations.
 */
template <class SourceGridView,
          class TargetElement = XT::Grid::extract_entity_t<SourceGridView>,
          size_t range_dim = 1,
          size_t range_dim_cols = 1,
          class RangeField = double>
class ReinterpretLocalizableFunction
    : public GridFunctionInterface<TargetElement, range_dim, range_dim_cols, RangeField>
{
  static_assert(XT::Grid::is_layer<SourceGridView>::value, "");
  using ThisType = ReinterpretLocalizableFunction<SourceGridView, TargetElement, range_dim, range_dim_cols, RangeField>;
  using BaseType = GridFunctionInterface<TargetElement, range_dim, range_dim_cols, RangeField>;

public:
  using BaseType::r;
  using BaseType::rC;
  using typename BaseType::R;
  using typename BaseType::LocalFunctionType;

  using SourceType = GridFunctionInterface<XT::Grid::extract_entity_t<SourceGridView>, r, rC, R>;

  static std::string static_id()
  {
    return BaseType::static_id() + ".reinterpret";
  }

  ReinterpretLocalizableFunction(const SourceType& source, const SourceGridView& source_grid_view)
    : BaseType(source.parameter_type())
    , source_(source)
    , source_grid_view_(source_grid_view)
  {
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<ReinterpretLocalfunction>(source_, source_grid_view_);
  }

  std::string name() const
  {
    return BaseType::name() + ".reinterpret";
  }

private:
  class ReinterpretLocalfunction : public ElementFunctionInterface<TargetElement, r, rC, R>
  {
    using BaseType = ElementFunctionInterface<TargetElement, r, rC, R>;

  public:
    using typename BaseType::D;
    using BaseType::d;
    using typename BaseType::DomainType;
    using typename BaseType::RangeReturnType;
    using typename BaseType::DerivativeRangeReturnType;

    ReinterpretLocalfunction(const SourceType& source, const SourceGridView& source_grid_view)
      : BaseType(source.parameter_type())
      , source_(source)
      , source_grid_view_(source_grid_view)
      , source_element_search_(source_grid_view_)
      , local_source_(source_.local_function())
      , source_element_which_contains_complete_target_element_(nullptr)
      , source_element_which_contains_some_point_of_target_element_(nullptr)
    {
    }

  protected:
    void post_bind(const TargetElement& target_element)
    {
      // See if we find a source element which contais target_element completely. Therefore
      // * collect all vertices
      const auto& reference_element = ReferenceElements<D, d>::general(target_element.geometry().type());
      const auto num_vertices = reference_element.size(1);
      if (vertices_.size() != num_vertices)
        vertices_.resize(num_vertices);
      for (int ii = 0; ii < num_vertices; ++ii)
        vertices_[ii] = target_element.geometry().global(reference_element.position(ii, 1));
      // * search for a source element for each vertex
      auto source_element_ptrs = source_element_search_(vertices_);
      DUNE_THROW_IF(
          source_element_ptrs.size() != num_vertices,
          Exceptions::reinterpretation_error,
          "source_element_ptrs.size() = " << source_element_ptrs.size() << "\n   num_vertices = " << num_vertices);
      // * and check if these are all the same
      if (num_vertices == 1) {
        source_element_which_contains_complete_target_element_ = std::move(source_element_ptrs[0]);
        source_element_which_contains_some_point_of_target_element_ = nullptr;
        local_source_->bind(*source_element_which_contains_complete_target_element_);
      } else {
        int is_same_as_first_one = 0;
        for (int ii = 1; ii < num_vertices; ++ii)
          if (source_element_ptrs[ii] == source_element_ptrs[0])
            ++is_same_as_first_one;
        if (is_same_as_first_one == num_vertices - 1) {
          source_element_which_contains_complete_target_element_ = std::move(source_element_ptrs[0]);
          source_element_which_contains_some_point_of_target_element_ = nullptr;
          local_source_->bind(*source_element_which_contains_complete_target_element_);
        } else {
          // We could not find a single source element which contains target_element completely.
          source_element_which_contains_complete_target_element_ = nullptr;
          source_element_which_contains_some_point_of_target_element_ = std::move(source_element_ptrs[0]);
          local_source_->bind(*source_element_which_contains_some_point_of_target_element_);
        }
      }
    } // ... post_bind(...)

  public:
    /**
      * \note In some special situations (e.g., if the target element is not completely contained in one source
      *       element), this may give inaccurate results.
      **/
    int order(const Common::Parameter& param = {}) const override final
    {
      DUNE_THROW_IF(!source_element_which_contains_complete_target_element_
                        && !source_element_which_contains_some_point_of_target_element_,
                    Exceptions::not_bound_to_an_element_yet,
                    "");
      return local_source_->order(param);
    }

    RangeReturnType evaluate(const DomainType& point_in_target_reference_element,
                             const Common::Parameter& param = {}) const
    {
      ensure_local_source_is_bound_for_this_point(point_in_target_reference_element);
      const auto point_in_global_coordinates = this->element().geometry().global(point_in_target_reference_element);
      const auto point_in_source_reference_element =
          local_source_->element().geometry().local(point_in_global_coordinates);
      return local_source_->evaluate(point_in_source_reference_element, param);
    }

    DerivativeRangeReturnType jacobian(const DomainType& point_in_target_reference_element,
                                       const Common::Parameter& param = {}) const
    {
      ensure_local_source_is_bound_for_this_point(point_in_target_reference_element);
      const auto point_in_global_coordinates = this->element().geometry().global(point_in_target_reference_element);
      const auto point_in_source_reference_element =
          local_source_->element().geometry().local(point_in_global_coordinates);
      return local_source_->jacobian(point_in_source_reference_element, param);
    }

    DerivativeRangeReturnType derivative(const std::array<size_t, d>& alpha,
                                         const DomainType& point_in_target_reference_element,
                                         const Common::Parameter& param = {}) const
    {
      ensure_local_source_is_bound_for_this_point(point_in_target_reference_element);
      const auto point_in_global_coordinates = this->element().geometry().global(point_in_target_reference_element);
      const auto point_in_source_reference_element =
          local_source_->element().geometry().local(point_in_global_coordinates);
      return local_source_->derivative(alpha, point_in_source_reference_element, param);
    }

  private:
    void ensure_local_source_is_bound_for_this_point(const DomainType& point_in_target_reference_element) const
    {
      if (source_element_which_contains_complete_target_element_)
        return;
      if (single_point_.size() != 1)
        single_point_.resize(1);
      single_point_[0] = this->element().geometry().global(point_in_target_reference_element);
      auto source_element_ptrs = source_element_search_(single_point_);
      DUNE_THROW_IF(source_element_ptrs.size() != 1,
                    Exceptions::reinterpretation_error,
                    "source_element_ptrs.size() = " << source_element_ptrs.size());
      source_element_which_contains_some_point_of_target_element_ = std::move(source_element_ptrs[0]);
      local_source_->bind(*source_element_which_contains_some_point_of_target_element_);
    } // ... ensure_local_source_is_bound_for_this_point(...)

    const SourceType& source_;
    const SourceGridView& source_grid_view_;
    mutable XT::Grid::EntityInlevelSearch<SourceGridView> source_element_search_;
    mutable std::unique_ptr<typename SourceType::LocalFunctionType> local_source_;
    mutable std::unique_ptr<XT::Grid::extract_entity_t<SourceGridView>>
        source_element_which_contains_complete_target_element_;
    mutable std::unique_ptr<XT::Grid::extract_entity_t<SourceGridView>>
        source_element_which_contains_some_point_of_target_element_;
    mutable std::vector<DomainType> vertices_;
    mutable std::vector<DomainType> single_point_;
  }; // class ReinterpretLocalfunction

  const SourceType& source_;
  const SourceGridView& source_grid_view_;
}; // class ReinterpretLocalizableFunction


/**
 * \brief Reinterprets a given localizable function (associated with a given grid layer), on a different grid layer
 *        [TargetElement has to be provided].
 *
 * Use as in
\code
auto reinterpreted_source = reinterpret<TargetElement>(source, source_grid_view);
\endcode
 *
 * \sa ReinterpretLocalizableFunction
 */
template <class TargetElement, class SourceGridView, size_t r, size_t rC, class R>
std::enable_if_t<XT::Grid::is_layer<SourceGridView>::value,
                 ReinterpretLocalizableFunction<SourceGridView, TargetElement, r, rC, R>>
reinterpret(const GridFunctionInterface<XT::Grid::extract_entity_t<SourceGridView>, r, rC, R>& source,
            const SourceGridView& source_grid_view)
{
  return ReinterpretLocalizableFunction<SourceGridView, TargetElement, r, rC, R>(source, source_grid_view);
}


/**
 * \brief Reinterprets a given localizable function (associated with a given grid layer), on a different grid layer
 *        [TargetElement is automatically deduced].
 *
 * \sa ReinterpretLocalizableFunction
 */
template <class SourceGridView, size_t r, size_t rC, class R, class TargetGridView>
std::enable_if_t<XT::Grid::is_layer<SourceGridView>::value && XT::Grid::is_layer<TargetGridView>::value,
                 ReinterpretLocalizableFunction<SourceGridView, XT::Grid::extract_entity_t<TargetGridView>, r, rC, R>>
reinterpret(const GridFunctionInterface<XT::Grid::extract_entity_t<SourceGridView>, r, rC, R>& source,
            const SourceGridView& source_grid_view,
            const TargetGridView& /*target_grid_view*/)
{
  return reinterpret<XT::Grid::extract_entity_t<TargetGridView>>(source, source_grid_view);
}


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTION_REINTERPRET_HH
