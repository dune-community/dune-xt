// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_XT_FUNCTIONS_COMPOSITION_HH
#define DUNE_XT_FUNCTIONS_COMPOSITION_HH

#include <dune/grid/yaspgrid.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/parallel/threadstorage.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/grid/search.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace internal {


template <class InnerType, class OuterType, class OuterGridViewType>
struct GeneralElementFunctionChooser
{
  using GridFunctionInterfaceType = GridFunctionInterface<typename InnerType::ElementType,
                                                          OuterType::range_dim,
                                                          OuterType::range_dim_cols,
                                                          typename OuterType::RangeFieldType>;
  using ElementType = typename GridFunctionInterfaceType::ElementType;
  using DomainFieldType = typename GridFunctionInterfaceType::DomainFieldType;
  using RangeFieldType = typename GridFunctionInterfaceType::RangeFieldType;
  using DomainType = typename GridFunctionInterfaceType::DomainType;
  using RangeType = typename GridFunctionInterfaceType::RangeType;
  using DerivativeRangeType = typename GridFunctionInterfaceType::DerivativeRangeType;
  static const size_t domain_dim = GridFunctionInterfaceType::domain_dim;
  static const size_t range_dim = GridFunctionInterfaceType::range_dim;
  static const size_t range_dim_cols = GridFunctionInterfaceType::range_dim_cols;

  class ElementFunction : public ElementFunctionInterface<ElementType, range_dim, range_dim_cols, RangeFieldType>
  {
    using BaseType = ElementFunctionInterface<ElementType, range_dim, range_dim_cols, RangeFieldType>;

  public:
    ElementFunction(const InnerType& inner_function,
                    const OuterType& outer_function,
                    std::shared_ptr<typename Grid::EntityInlevelSearch<OuterGridViewType>>& element_search)
      : BaseType()
      , inner_function_(inner_function)
      , outer_function_(outer_function)
      , element_search_(element_search)
      , local_inner_function_(inner_function_.local_function())
      , local_outer_function_(outer_function_.local_function())
    {
    }

    ElementFunction(const ElementFunction& /*other*/) = delete;

    ElementFunction& operator=(const ElementFunction& /*other*/) = delete;

  protected:
    void post_bind(const ElementType& element)
    {
      local_inner_function_->bind(element);
    }

  public:
    int order(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return 2;
    }

    RangeType evaluate(const DomainType& xx) const override final
    {
      // evaluate inner function
      const auto inner_value = local_inner_function_->evaluate(xx);
      // find element on outer grid the value of inner function belongs to
      const auto element_ptrs = (*element_search_)(std::vector<typename OuterType::DomainType>(1, inner_value));
      const auto& element_ptr = element_ptrs[0];
      if (element_ptr == nullptr)
        DUNE_THROW(Dune::InvalidStateException,
                   "Could not find element, maybe inner function does not map to the domain of outer function");
      const auto& outer_element = *element_ptr;
      // evaluate outer function
      local_outer_function_->bind(outer_element);
      return local_outer_function_->evaluate(outer_element.geometry().local(inner_value));
    }

    DerivativeRangeType jacobian(const DomainType& /*xx*/) const override final
    {
      DUNE_THROW(Dune::NotImplemented, "");
    }

  private:
    const InnerType& inner_function_;
    const OuterType& outer_function_;
    std::shared_ptr<typename Grid::EntityInlevelSearch<OuterGridViewType>>& element_search_;
    std::unique_ptr<typename InnerType::LocalFunctionType>& local_inner_function_;
    std::unique_ptr<typename OuterType::LocalFunctionType>& local_outer_function_;
    const ElementType& element_;
  }; // class ElementFunction
}; // GeneralElementFunctionChooser


template <class InnerType, class OuterType, class OuterGridViewType>
struct ElementFunctionForGlobalChooser
{
  using GridFunctionInterfaceType = GridFunctionInterface<typename InnerType::ElementType,
                                                          OuterType::range_dim,
                                                          OuterType::range_dim_cols,
                                                          typename OuterType::RangeFieldType>;
  using ElementType = typename GridFunctionInterfaceType::ElementType;
  using DomainFieldType = typename GridFunctionInterfaceType::DomainFieldType;
  using RangeFieldType = typename GridFunctionInterfaceType::RangeFieldType;
  using DomainType = typename InnerType::DomainType;
  using RangeType = typename OuterType::RangeType;
  using DerivativeRangeType = typename OuterType::DerivativeRangeType;
  static const size_t domain_dim = GridFunctionInterfaceType::domain_dim;
  static const size_t range_dim = GridFunctionInterfaceType::range_dim;
  static const size_t range_dim_cols = GridFunctionInterfaceType::range_dim_cols;

  class ElementFunction : public ElementFunctionInterface<ElementType, range_dim, range_dim_cols, RangeFieldType>
  {
    using BaseType = ElementFunctionInterface<ElementType, range_dim, range_dim_cols, RangeFieldType>;

  public:
    ElementFunction(const InnerType& localizable_function,
                    const OuterType& global_function,
                    std::shared_ptr<typename Grid::EntityInlevelSearch<OuterGridViewType>>& /*element_search*/)
      : BaseType()
      , localizable_function_(localizable_function)
      , global_function_(global_function)
    {
    }

    ElementFunction(const ElementFunction& /*other*/) = delete;

    ElementFunction& operator=(const ElementFunction& /*other*/) = delete;

    int order(const XT::Common::Parameter& param = {}) const override final
    {
      return global_function_.order(param) * localizable_function_.local_function(element_)->order(param);
    }

    RangeType evaluate(const DomainType& xx, const XT::Common::Parameter& param = {}) const override final
    {
      return global_function_.evaluate(localizable_function_.local_function(element_)->evaluate(xx, param), param);
    }

    DerivativeRangeType jacobian(const DomainType& /*xx*/,
                                 const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      DUNE_THROW(Dune::NotImplemented, "");
    }

  private:
    const InnerType& localizable_function_;
    const OuterType& global_function_;
    const ElementType& element_;
  }; // class ElementFunction
}; // ElementFunctionForGlobalChooser

template <class InnerType, class OuterType, class OuterGridViewType>
struct ElementFunctionChooser
{
  using ElementFunctionType = typename std::
      conditional<std::is_base_of<Functions::FunctionInterface<OuterType::domain_dim,
                                                               OuterType::range_dim,
                                                               OuterType::range_dim_cols,
                                                               typename OuterType::RangeFieldType>,
                                  OuterType>::value,
                  typename ElementFunctionForGlobalChooser<InnerType, OuterType, OuterGridViewType>::ElementFunction,
                  typename GeneralElementFunctionChooser<InnerType, OuterType, OuterGridViewType>::ElementFunction>::
          type;
};


} // namespace internal


template <class InnerType, class OuterType, class OuterGridViewType = typename Dune::YaspGrid<1>::LeafGridView>
class CompositionFunction : public GridFunctionInterface<typename InnerType::ElementType,
                                                         OuterType::range_dim,
                                                         OuterType::range_dim_cols,
                                                         typename OuterType::RangeFieldType>
{
  using BaseType = GridFunctionInterface<typename InnerType::ElementType,
                                         OuterType::range_dim,
                                         OuterType::range_dim_cols,
                                         typename OuterType::RangeFieldType>;
  using ThisType = CompositionFunction<InnerType, OuterType, OuterGridViewType>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using BaseType::domain_dim;
  using BaseType::range_dim;
  using BaseType::range_dim_cols;

private:
  using ElementFunction =
      typename internal::ElementFunctionChooser<InnerType, OuterType, OuterGridViewType>::ElementFunctionType;

public:
  using typename BaseType::LocalFunctionType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".composition";
  }

  CompositionFunction(const InnerType inner_function,
                      const OuterType outer_function,
                      const OuterGridViewType outer_grid_view,
                      const std::string nm = static_id())
    : inner_function_(inner_function)
    , outer_function_(outer_function)
    , element_search_(std::make_shared<typename Grid::EntityInlevelSearch<OuterGridViewType>>(outer_grid_view))
    , name_(nm)
  {
  }

  // constructor without grid view, only makes sense if OuterType is derived from FunctionInterface
  CompositionFunction(const InnerType local_func, const OuterType global_func, const std::string nm = static_id())
    : inner_function_(local_func)
    , outer_function_(global_func)
    , name_(nm)
  {
    static_assert(std::is_base_of<XT::Functions::FunctionInterface<OuterType::domain_dim,
                                                                   OuterType::range_dim,
                                                                   OuterType::range_dim_cols,
                                                                   typename OuterType::RangeFieldType>,
                                  OuterType>::value,
                  "OuterType has to be derived from FunctionInterface if no GridView is provided");
  }

  CompositionFunction(const ThisType& other) = default;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  virtual std::string type() const override
  {
    return BaseType::static_id() + ".composition";
  }

  virtual std::string name() const override
  {
    return name_;
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return Common::make_unique<ElementFunction>(inner_function_, outer_function_, element_search_);
  } // ... local_function(...)


private:
  const InnerType inner_function_;
  const OuterType outer_function_;
  mutable std::shared_ptr<typename Grid::EntityInlevelSearch<OuterGridViewType>> element_search_;
  std::string name_;
}; // class Composition


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_COMPOSITION_HH
