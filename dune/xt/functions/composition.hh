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
#if 0
#include <dune/grid/yaspgrid.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/parallel/threadstorage.hh>
#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/grid/search.hh>


namespace Dune {
namespace XT {
namespace Functions {
namespace internal {

template <class InnerType, class OuterType, class OuterGridViewType>
struct GeneralLocalfunctionChooser
{
  typedef LocalizableFunctionInterface<typename InnerType::EntityType,
                                       typename InnerType::DomainFieldType,
                                       InnerType::dimDomain,
                                       typename OuterType::RangeFieldType,
                                       OuterType::dimRange,
                                       OuterType::dimRangeCols>
      LocalizableFunctionInterfaceType;
  typedef typename LocalizableFunctionInterfaceType::EntityType EntityType;
  typedef typename LocalizableFunctionInterfaceType::DomainFieldType DomainFieldType;
  typedef typename LocalizableFunctionInterfaceType::RangeFieldType RangeFieldType;
  typedef typename LocalizableFunctionInterfaceType::DomainType DomainType;
  typedef typename LocalizableFunctionInterfaceType::RangeType RangeType;
  typedef typename LocalizableFunctionInterfaceType::JacobianRangeType JacobianRangeType;
  static const size_t dimDomain = LocalizableFunctionInterfaceType::dimDomain;
  static const size_t dimRange = LocalizableFunctionInterfaceType::dimRange;
  static const size_t dimRangeCols = LocalizableFunctionInterfaceType::dimRangeCols;

  class Localfunction
      : public LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols>
  {
    typedef LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols>
        BaseType;

  public:
    Localfunction(const EntityType& entity,
                  const InnerType& inner_function,
                  const OuterType& outer_function,
                  std::shared_ptr<typename Grid::EntityInlevelSearch<OuterGridViewType>>& entity_search)
      : BaseType(entity)
      , inner_function_(inner_function)
      , outer_function_(outer_function)
      , entity_search_(entity_search)
      , entity_(entity)
    {
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override
    {
      return 2;
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const override
    {
      // evaluate inner function
      const auto inner_value = inner_function_.local_function(entity_)->evaluate(xx);
      // find entity on outer grid the value of inner function belongs to
      const auto entity_ptrs = (*entity_search_)(std::vector<typename OuterType::DomainType>(1, inner_value));
      const auto& entity_ptr = entity_ptrs[0];
      if (entity_ptr == nullptr)
        DUNE_THROW(Dune::InvalidStateException,
                   "Could not find entity, maybe inner function does not map to the domain of outer function");
      const auto& outer_entity = *entity_ptr;
      // evaluate outer function
      outer_function_.local_function(outer_entity)->evaluate(outer_entity.geometry().local(inner_value), ret);
    }

    virtual void jacobian(const DomainType& /*xx*/, JacobianRangeType& /*ret*/) const override
    {
      DUNE_THROW(Dune::NotImplemented, "");
    }

  private:
    const InnerType& inner_function_;
    const OuterType& outer_function_;
    std::shared_ptr<typename Grid::EntityInlevelSearch<OuterGridViewType>>& entity_search_;
    const EntityType& entity_;
  }; // class Localfunction
}; // GeneralLocalfunctionChooser

template <class InnerType, class OuterType, class OuterGridViewType>
struct LocalfunctionForGlobalChooser
{
  typedef LocalizableFunctionInterface<typename InnerType::EntityType,
                                       typename InnerType::DomainFieldType,
                                       InnerType::dimDomain,
                                       typename OuterType::RangeFieldType,
                                       OuterType::dimRange,
                                       OuterType::dimRangeCols>
      LocalizableFunctionInterfaceType;
  typedef typename LocalizableFunctionInterfaceType::EntityType EntityType;
  typedef typename LocalizableFunctionInterfaceType::DomainFieldType DomainFieldType;
  typedef typename LocalizableFunctionInterfaceType::RangeFieldType RangeFieldType;
  typedef typename LocalizableFunctionInterfaceType::DomainType DomainType;
  typedef typename LocalizableFunctionInterfaceType::RangeType RangeType;
  typedef typename LocalizableFunctionInterfaceType::JacobianRangeType JacobianRangeType;
  static const size_t dimDomain = LocalizableFunctionInterfaceType::dimDomain;
  static const size_t dimRange = LocalizableFunctionInterfaceType::dimRange;
  static const size_t dimRangeCols = LocalizableFunctionInterfaceType::dimRangeCols;

  class Localfunction
      : public LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols>
  {
    typedef LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols>
        BaseType;

  public:
    Localfunction(const EntityType& entity,
                  const InnerType& localizable_function,
                  const OuterType& global_function,
                  std::shared_ptr<typename Grid::EntityInlevelSearch<OuterGridViewType>>& /*entity_search*/)
      : BaseType(entity)
      , localizable_function_(localizable_function)
      , global_function_(global_function)
      , entity_(entity)
    {
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order(const XT::Common::Parameter& mu = {}) const override
    {
      return global_function_.order(mu) * localizable_function_.local_function(entity_)->order(mu);
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret, const XT::Common::Parameter& mu = {}) const override
    {
      global_function_.evaluate(localizable_function_.local_function(entity_)->evaluate(xx, mu), ret, mu);
    }

    virtual void jacobian(const DomainType& /*xx*/,
                          JacobianRangeType& /*ret*/,
                          const XT::Common::Parameter& /*mu*/ = {}) const override
    {
      DUNE_THROW(Dune::NotImplemented, "");
    }

  private:
    const InnerType& localizable_function_;
    const OuterType& global_function_;
    const EntityType& entity_;
  }; // class Localfunction
}; // LocalfunctionForGlobalChooser

template <class InnerType, class OuterType, class OuterGridViewType>
struct LocalfunctionChooser
{
  typedef typename std::
      conditional<std::is_base_of<Functions::GlobalFunctionInterface<typename OuterType::EntityType,
                                                                     typename OuterType::DomainFieldType,
                                                                     OuterType::dimDomain,
                                                                     typename OuterType::RangeFieldType,
                                                                     OuterType::dimRange,
                                                                     OuterType::dimRangeCols>,
                                  OuterType>::value,
                  typename LocalfunctionForGlobalChooser<InnerType, OuterType, OuterGridViewType>::Localfunction,
                  typename GeneralLocalfunctionChooser<InnerType, OuterType, OuterGridViewType>::Localfunction>::type
          LocalfunctionType;
};


} // namespace internal


template <class InnerType, class OuterType, class OuterGridViewType = typename Dune::YaspGrid<1>::LeafGridView>
class CompositionFunction : public LocalizableFunctionInterface<typename InnerType::EntityType,
                                                                typename InnerType::DomainFieldType,
                                                                InnerType::dimDomain,
                                                                typename OuterType::RangeFieldType,
                                                                OuterType::dimRange,
                                                                OuterType::dimRangeCols>
{
  typedef LocalizableFunctionInterface<typename InnerType::EntityType,
                                       typename InnerType::DomainFieldType,
                                       InnerType::dimDomain,
                                       typename OuterType::RangeFieldType,
                                       OuterType::dimRange,
                                       OuterType::dimRangeCols>
      BaseType;
  typedef CompositionFunction<InnerType, OuterType, OuterGridViewType> ThisType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;

private:
  typedef
      typename internal::LocalfunctionChooser<InnerType, OuterType, OuterGridViewType>::LocalfunctionType Localfunction;

public:
  using typename BaseType::LocalfunctionType;

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
    , entity_search_(std::make_shared<typename Grid::EntityInlevelSearch<OuterGridViewType>>(outer_grid_view))
    , name_(nm)
  {
  }

  // constructor without grid view, only makes sense if OuterType is derived from GlobalFunctionInterface
  CompositionFunction(const InnerType local_func, const OuterType global_func, const std::string nm = static_id())
    : inner_function_(local_func)
    , outer_function_(global_func)
    , name_(nm)
  {
    static_assert(std::is_base_of<XT::Functions::GlobalFunctionInterface<typename OuterType::EntityType,
                                                                         typename OuterType::DomainFieldType,
                                                                         OuterType::dimDomain,
                                                                         typename OuterType::RangeFieldType,
                                                                         OuterType::dimRange,
                                                                         OuterType::dimRangeCols>,
                                  OuterType>::value,
                  "OuterType has to be derived from GlobalFunctionInterface if no GridView is provided");
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

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override
  {
    return Common::make_unique<Localfunction>(entity, inner_function_, outer_function_, entity_search_);
  } // ... local_function(...)

private:
  const InnerType inner_function_;
  const OuterType outer_function_;
  mutable std::shared_ptr<typename Grid::EntityInlevelSearch<OuterGridViewType>> entity_search_;
  std::string name_;
}; // class Composition


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif
#endif // DUNE_XT_FUNCTIONS_COMPOSITION_HH
