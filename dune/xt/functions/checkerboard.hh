// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2012 - 2017)
//   Rene Milk       (2013 - 2016)
//   Tobias Leibner  (2014 - 2016)

#ifndef DUNE_XT_FUNCTIONS_CHECKERBOARD_HH
#define DUNE_XT_FUNCTIONS_CHECKERBOARD_HH

#include <cmath>
#include <memory>
#include <vector>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/debug.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/fvector.hh>

#include <dune/xt/functions/affine.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/expression.hh>

#include "interfaces.hh"

namespace Dune {
namespace XT {
namespace Functions {


// forward
template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols,
          class LocalizableFunctionImp =
              ConstantFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>>
class CheckerboardFunction;


namespace internal {


template <class FunctionImp,
          bool is_localizable = std::is_base_of<LocalizableFunctionInterface<typename FunctionImp::EntityType,
                                                                             typename FunctionImp::DomainFieldType,
                                                                             FunctionImp::dimDomain,
                                                                             typename FunctionImp::RangeFieldType,
                                                                             FunctionImp::dimRange,
                                                                             FunctionImp::dimRangeCols>,
                                                FunctionImp>::value>
struct CheckerboardInterfaceChooser;

template <class FunctionImp>
struct CheckerboardInterfaceChooser<FunctionImp, true>
{
  typedef LocalizableFunctionInterface<typename FunctionImp::EntityType,
                                       typename FunctionImp::DomainFieldType,
                                       FunctionImp::dimDomain,
                                       typename FunctionImp::RangeFieldType,
                                       FunctionImp::dimRange,
                                       FunctionImp::dimRangeCols>
      Type;
};

template <class FunctionImp>
struct CheckerboardInterfaceChooser<FunctionImp, false>
{

  typedef LocalizableFluxFunctionInterface<typename FunctionImp::EntityType,
                                           typename FunctionImp::DomainFieldType,
                                           FunctionImp::dimDomain,
                                           typename FunctionImp::StateType,
                                           0,
                                           typename FunctionImp::RangeFieldType,
                                           FunctionImp::dimRange,
                                           FunctionImp::dimRangeCols>
      Type;
};

template <class CheckerboardFunctionType>
class CheckerboardFunctionFactory
{
public:
  static std::unique_ptr<CheckerboardFunctionType> create(const Common::Configuration /*config*/,
                                                          const std::string /*sub_name*/)
  {
    DUNE_THROW(Dune::NotImplemented,
               "CheckerboardFunctionFactory is not implemented for this LocalizableFunctionType!");
    return std::unique_ptr<CheckerboardFunctionType>();
  }
}; // class CheckerboardFunctionFactory< ... >

template <>
template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols>
class CheckerboardFunctionFactory<CheckerboardFunction<EntityImp,
                                                       DomainFieldImp,
                                                       domainDim,
                                                       RangeFieldImp,
                                                       rangeDim,
                                                       rangeDimCols,
                                                       ConstantFunction<EntityImp,
                                                                        DomainFieldImp,
                                                                        domainDim,
                                                                        RangeFieldImp,
                                                                        rangeDim,
                                                                        rangeDimCols>>>
{
  typedef ConstantFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      ConstantFunctionType;
  typedef CheckerboardFunction<EntityImp,
                               DomainFieldImp,
                               domainDim,
                               RangeFieldImp,
                               rangeDim,
                               rangeDimCols,
                               ConstantFunctionType>
      CheckerboardFunctionType;

public:
  static std::unique_ptr<CheckerboardFunctionType> create(const Common::Configuration config,
                                                          const std::string sub_name)
  {
    // get correct config
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = CheckerboardFunctionType::default_config();
    // calculate number of values and get values
    auto num_elements =
        cfg.get("num_elements", default_cfg.get<Common::FieldVector<size_t, domainDim>>("num_elements"), domainDim);
    size_t num_values = 1;
    for (size_t ii = 0; ii < num_elements.size(); ++ii)
      num_values *= num_elements[ii];
    std::vector<ConstantFunctionType> values;
    auto values_range = cfg.get(
        "values", default_cfg.get<std::vector<typename CheckerboardFunctionType::RangeType>>("values"), num_values);
    for (size_t ii = 0; ii < values_range.size(); ++ii)
      values.emplace_back(ConstantFunctionType(values_range[ii], "constant in tile " + Common::to_string(ii)));
    // create
    return Common::make_unique<CheckerboardFunctionType>(
        cfg.get("lower_left", default_cfg.get<FieldVector<DomainFieldImp, domainDim>>("lower_left"), domainDim),
        cfg.get("upper_right", default_cfg.get<FieldVector<DomainFieldImp, domainDim>>("upper_right"), domainDim),
        std::move(num_elements),
        std::move(values),
        cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)
}; // class CheckerboardFunctionFactory< ... >


} // namespace internal


template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols,
          class LocalizableFunctionImp>
class CheckerboardFunction : public internal::CheckerboardInterfaceChooser<LocalizableFunctionImp>::Type
{
  //  static_assert(is_localizable_function<LocalizableFunctionImp>::value
  //                    || is_localizable_flux_function <LocalizableFunctionImp::value,
  //                "LocalizableFunctionImp needs to be derived from XT::Localizable(Flux)FunctionInterface!");
  static_assert(domainDim <= 3, "Not implemented for dimDomain > 3 (see find_subdomain method)!");
  typedef typename internal::CheckerboardInterfaceChooser<LocalizableFunctionImp>::Type BaseType;
  typedef CheckerboardFunction<EntityImp,
                               DomainFieldImp,
                               domainDim,
                               RangeFieldImp,
                               rangeDim,
                               rangeDimCols,
                               LocalizableFunctionImp>
      ThisType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::LocalfunctionType;

  using BaseType::dimDomain;

  static const bool available = true;

  typedef LocalizableFunctionImp LocalizableFunctionType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".checkerboard";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["lower_left"] = "[0.0 0.0 0.0]";
    config["upper_right"] = "[1.0 1.0 1.0]";
    config["num_elements"] = "[2 2 2]";
    config["values"] = "[1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0]";
    config["name"] = static_id();
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  }

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    return internal::CheckerboardFunctionFactory<ThisType>::create(config, sub_name);
  } // ... create(...)

  CheckerboardFunction(const DomainType& lowerLeft,
                       const DomainType& upperRight,
                       const FieldVector<size_t, dimDomain>& numElements,
                       const std::vector<LocalizableFunctionType>& values,
                       const std::string nm = static_id())
    : lowerLeft_(new DomainType(lowerLeft))
    , upperRight_(new DomainType(upperRight))
    , numElements_(new FieldVector<size_t, dimDomain>(numElements))
    , name_(nm)
  {
    for (size_t ii = 0; ii < values.size(); ++ii)
      values_.emplace_back(new LocalizableFunctionType(values[ii]));
#ifndef NDEBUG
    // checks
    size_t totalSubdomains = 1;
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      const auto& ll = (*lowerLeft_)[dd];
      const auto& ur = (*upperRight_)[dd];
      const auto& ne = (*numElements_)[dd];
      if (!(ll < ur))
        DUNE_THROW(Dune::RangeError, "lowerLeft has to be elementwise smaller than upperRight!");
      totalSubdomains *= ne;
    }
    if (values_.size() < totalSubdomains)
      DUNE_THROW(Dune::RangeError,
                 "values too small (is " << values_.size() << ", should be " << totalSubdomains << ")");
#endif
  } // Checkerboard(...)

  CheckerboardFunction(const ThisType& other) = default;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  virtual std::string type() const override
  {
    return BaseType::static_id() + ".checkerboard";
  }

  virtual std::string name() const override
  {
    return name_;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override
  {
    const size_t subdomain = find_subdomain(entity);
    return values_[subdomain]->local_function(entity);
  }

private:
  size_t find_subdomain(const EntityType& entity) const
  {
    // decide on the subdomain the center of the entity belongs to
    const auto center = entity.geometry().center();
    std::vector<size_t> whichPartition(dimDomain, 0);
    const auto& ll = *lowerLeft_;
    const auto& ur = *upperRight_;
    const auto& ne = *numElements_;
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      // for points that are on upperRight_[d], this selects one partition too much
      // so we need to cap this
      whichPartition[dd] =
          std::min(size_t(std::floor(ne[dd] * ((center[dd] - ll[dd]) / (ur[dd] - ll[dd])))), ne[dd] - 1);
    }
    size_t subdomain = 0;
    if (dimDomain == 1)
      subdomain = whichPartition[0];
    else if (dimDomain == 2)
      subdomain = whichPartition[0] + whichPartition[1] * ne[0];
    else
      subdomain = whichPartition[0] + whichPartition[1] * ne[0] + whichPartition[2] * ne[1] * ne[0];
    return subdomain;
  } // ... find_subdomain(...)

  std::shared_ptr<const DomainType> lowerLeft_;
  std::shared_ptr<const DomainType> upperRight_;
  std::shared_ptr<const FieldVector<size_t, dimDomain>> numElements_;
  std::vector<std::shared_ptr<const LocalizableFunctionType>> values_;
  std::string name_;
}; // class CheckerboardFunction


} // namespace Functions
} // namespace XT
} // namespace Dune


#include "checkerboard.lib.hh"


#endif // DUNE_XT_FUNCTIONS_CHECKERBOARD_HH
