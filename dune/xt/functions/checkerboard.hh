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
#include <dune/xt/functions/expression.hh>

#include "interfaces.hh"

namespace Dune {
namespace XT {
namespace Functions {

template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class CheckerboardFunction
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      BaseType;
  typedef CheckerboardFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;

  class Localfunction
      : public LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
  {
    typedef LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
        BaseType;

  public:
    typedef typename BaseType::EntityType EntityType;

    typedef typename BaseType::DomainType DomainType;
    typedef typename BaseType::RangeFieldType RangeFieldType;
    typedef typename BaseType::RangeType RangeType;
    typedef typename BaseType::JacobianRangeType JacobianRangeType;

    Localfunction(const EntityType& ent, const RangeType value)
      : BaseType(ent)
      , value_(value)
    {
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order() const override
    {
      return 0;
    }

    virtual void evaluate(const DomainType& DXTC_DEBUG_ONLY(xx),
                          RangeType& ret,
                          const Common::Parameter& /*mu*/ = Common::Parameter()) const override
    {
      assert(this->is_a_valid_point(xx));
      ret = value_;
    }

    virtual void jacobian(const DomainType& DXTC_DEBUG_ONLY(xx),
                          JacobianRangeType& ret,
                          const Common::Parameter& /*mu*/ = Common::Parameter()) const override
    {
      assert(this->is_a_valid_point(xx));
      clear_jacobian<rangeDim, rangeDimCols>()(ret);
    }

  private:
    template <size_t _r, size_t _rC, class Anything = void>
    struct clear_jacobian
    {
      void operator()(JacobianRangeType& ret)
      {
        for (auto& col_jacobian : ret)
          col_jacobian *= 0.;
      }
    };

    template <class Anything>
    struct clear_jacobian<1, 1, Anything>
    {
      void operator()(JacobianRangeType& ret)
      {
        ret *= 0.;
      }
    };

    const RangeType value_;
  }; // class Localfunction

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = BaseType::dimDomain;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::RangeType RangeType;

  static const bool available = true;

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
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    // get correct config
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    // calculate number of values and get values
    auto num_elements =
        cfg.get("num_elements", default_cfg.get<Common::FieldVector<size_t, dimDomain>>("num_elements"), dimDomain);
    size_t num_values = 1;
    for (size_t ii = 0; ii < num_elements.size(); ++ii)
      num_values *= num_elements[ii];
    std::vector<RangeType> values(num_values);
    auto values_rf = cfg.get("values", default_cfg.get<std::vector<RangeFieldType>>("values"), num_values);
    for (size_t ii = 0; ii < values_rf.size(); ++ii)
      values[ii] = RangeType(values_rf[ii]);
    // create
    return Common::make_unique<ThisType>(
        cfg.get(
            "lower_left", default_cfg.get<Common::FieldVector<DomainFieldType, dimDomain>>("lower_left"), dimDomain),
        cfg.get(
            "upper_right", default_cfg.get<Common::FieldVector<DomainFieldType, dimDomain>>("upper_right"), dimDomain),
        std::move(num_elements),
        std::move(values),
        cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)

  CheckerboardFunction(const Common::FieldVector<DomainFieldType, dimDomain>& lowerLeft,
                       const Common::FieldVector<DomainFieldType, dimDomain>& upperRight,
                       const Common::FieldVector<size_t, dimDomain>& numElements,
                       const std::vector<RangeType>& values,
                       const std::string nm = static_id())
    : lowerLeft_(new Common::FieldVector<DomainFieldType, dimDomain>(lowerLeft))
    , upperRight_(new Common::FieldVector<DomainFieldType, dimDomain>(upperRight))
    , numElements_(new Common::FieldVector<size_t, dimDomain>(numElements))
    , values_(new std::vector<RangeType>(values))
    , name_(nm)
  {
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
    if (values_->size() < totalSubdomains)
      DUNE_THROW(Dune::RangeError,
                 "values too small (is " << values_->size() << ", should be " << totalSubdomains << ")");
  } // Checkerboard(...)

#if !DUNE_XT_WITH_PYTHON_BINDINGS
  CheckerboardFunction(const ThisType& other) = default;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;
#endif // DUNE_XT_WITH_PYTHON_BINDINGS

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
    // return the component that belongs to the subdomain
    return std::unique_ptr<Localfunction>(new Localfunction(entity, (*values_)[subdomain]));
  } // ... local_function(...)

private:
  std::shared_ptr<const Common::FieldVector<DomainFieldType, dimDomain>> lowerLeft_;
  std::shared_ptr<const Common::FieldVector<DomainFieldType, dimDomain>> upperRight_;
  std::shared_ptr<const Common::FieldVector<size_t, dimDomain>> numElements_;
  std::shared_ptr<const std::vector<RangeType>> values_;
  std::string name_;
}; // class CheckerboardFunction

#if HAVE_DUNE_XT_LA

// forward ////////////////////////
template <class LocalizableFunctionImp,
          class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class FunctionCheckerboardFunction;

// configs /////////////////////////
template <class LocalizableFunctionImp>
struct CheckerboardConfigProvider
{
  static inline Common::Configuration default_config(const std::string /*sub_name*/ = "")
  {
    DUNE_THROW(Dune::NotImplemented, "No default config available for this LocalizableFunctionImp");
    return Common::Configuration();
  }
};

template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols>
struct CheckerboardConfigProvider<ExpressionFunction<EntityImp,
                                                     DomainFieldImp,
                                                     domainDim,
                                                     RangeFieldImp,
                                                     rangeDim,
                                                     rangeDimCols>>
{
  static inline Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["lower_left"] = "[0.0 0.0 0.0]";
    config["upper_right"] = "[1.0 1.0 1.0]";
    config["num_elements"] = "[2 2 2]";
    config["variable"] = "u";
    config["values"] = "[1.0*u[0] 2.0*u[0] 3.0*sin(u[0]) 4.0 5.0 6.0*cos(u[0]) 7.0 8.0]";
    config["name"] = "expressioncheckerboard";
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... checkerboard_default_config< ExpressionFunction< ... > >(...)
};

template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols>
struct CheckerboardConfigProvider<AffineFunction<EntityImp,
                                                 DomainFieldImp,
                                                 domainDim,
                                                 RangeFieldImp,
                                                 rangeDim,
                                                 rangeDimCols>>
{
  static inline Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["lower_left"] = "[0.0 0.0 0.0]";
    config["upper_right"] = "[1.0 1.0 1.0]";
    config["num_elements"] = "[2 1 1]";
    config["A.0"] = "[1 0 0; 0 1 0; 0 0 1]";
    config["A.1"] = "[1 2 3; 4 5 6; 7 8 9]";
    config["b.0"] = "[1 2 3]";
    config["b.1"] = "[0 0 0]";
    config["name"] = "affinecheckerboard";
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... checkerboard_default_config< Affine< ... > >(...)
};

// factory methods ////////////////////
template <class LocalizableFunctionImp,
          class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class FunctionCheckerboardFunctionFactory
{
  typedef FunctionCheckerboardFunction<LocalizableFunctionImp,
                                       EntityImp,
                                       DomainFieldImp,
                                       domainDim,
                                       RangeFieldImp,
                                       rangeDim,
                                       rangeDimCols>
      FunctionCheckerboardType;

public:
  static std::unique_ptr<FunctionCheckerboardType> create(const Common::Configuration /*config*/,
                                                          const std::string /*sub_name*/)
  {
    DUNE_THROW(Dune::NotImplemented,
               "FunctionCheckerboardFactory is not implemented for this LocalizableFunctionType!");
    return std::unique_ptr<FunctionCheckerboardType>();
  }
}; // class FunctionCheckerboardFactory< ... >


template <class ExpressionDomainFieldImp,
          size_t expressionDomainDim,
          class ExpressionEntityImp,
          class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols>
class FunctionCheckerboardFunctionFactory<ExpressionFunction<ExpressionEntityImp,
                                                             ExpressionDomainFieldImp,
                                                             expressionDomainDim,
                                                             RangeFieldImp,
                                                             rangeDim,
                                                             rangeDimCols>,
                                          EntityImp,
                                          DomainFieldImp,
                                          domainDim,
                                          RangeFieldImp,
                                          rangeDim,
                                          rangeDimCols>
{
  typedef ExpressionFunction<ExpressionEntityImp,
                             ExpressionDomainFieldImp,
                             expressionDomainDim,
                             RangeFieldImp,
                             rangeDim,
                             rangeDimCols>
      ExpressionFunctionType;
  typedef FunctionCheckerboardFunction<ExpressionFunctionType,
                                       EntityImp,
                                       DomainFieldImp,
                                       domainDim,
                                       RangeFieldImp,
                                       rangeDim,
                                       rangeDimCols>
      FunctionCheckerboardType;

public:
  static std::unique_ptr<FunctionCheckerboardType>
  create(const Common::Configuration config = CheckerboardConfigProvider<ExpressionFunctionType>::default_config(),
         const std::string sub_name = FunctionCheckerboardType::static_id())
  {
    // get correct config
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = CheckerboardConfigProvider<ExpressionFunctionType>::default_config();
    // calculate number of values and get values
    auto num_elements =
        cfg.get("num_elements", default_cfg.get<Common::FieldVector<size_t, domainDim>>("num_elements"), domainDim);
    size_t num_values = 1;
    for (size_t ii = 0; ii < num_elements.size(); ++ii)
      num_values *= num_elements[ii];
    const std::string variable = cfg.get("variable", default_cfg.get<std::string>("variable"));
    std::vector<ExpressionFunctionType> values;
    if (config.has_key("values.0")) { // get every value from its own config entry
      try { // get value as matrix
        std::vector<std::string> row_as_std_vector(rangeDimCols);
        for (size_t ii = 0; ii < num_values; ++ii) {
          const auto values_matrix =
              cfg.get<Dune::DynamicMatrix<std::string>>("values." + Common::to_string(ii), rangeDim, rangeDimCols);
          typename ExpressionFunctionType::ExpressionStringVectorType expression_vector;
          for (size_t rr = 0; rr < rangeDim; ++rr) {
            const auto row = values_matrix[rr];
            for (size_t cc = 0; cc < rangeDimCols; ++cc)
              row_as_std_vector[cc] = row[cc];
            expression_vector.emplace_back(row_as_std_vector);
          }
          values.emplace_back(ExpressionFunctionType(variable, expression_vector));
        }
      } catch (const std::exception& e) {
        if (rangeDimCols == 1) { // get value as vector
          for (size_t ii = 0; ii < num_values; ++ii) {
            const auto values_vector = cfg.get<std::vector<std::string>>("values." + Common::to_string(ii), rangeDim);
            values.emplace_back(ExpressionFunctionType(variable, values_vector));
          }
        } else {
          DUNE_THROW(Common::Exceptions::conversion_error, "");
          std::cout << e.what() << std::endl;
        }
      }
    } else {
      // get values as a vector of scalars
      auto values_rf = cfg.get("values", default_cfg.get<std::vector<std::string>>("values"), num_values);
      for (size_t ii = 0; ii < values_rf.size(); ++ii)
        values.emplace_back(ExpressionFunctionType(variable, values_rf[ii]));
    }
    // create
    return Common::make_unique<FunctionCheckerboardType>(
        cfg.get("lower_left", default_cfg.get<Common::FieldVector<DomainFieldImp, domainDim>>("lower_left"), domainDim),
        cfg.get(
            "upper_right", default_cfg.get<Common::FieldVector<DomainFieldImp, domainDim>>("upper_right"), domainDim),
        num_elements,
        values,
        cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)
}; // class FunctionCheckerboardFactory< Expression, ... >

template <class AffineDomainFieldImp,
          size_t affineDomainDim,
          class AffineEntityImp,
          class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols>
class FunctionCheckerboardFunctionFactory<AffineFunction<AffineEntityImp,
                                                         AffineDomainFieldImp,
                                                         affineDomainDim,
                                                         RangeFieldImp,
                                                         rangeDim,
                                                         rangeDimCols>,
                                          EntityImp,
                                          DomainFieldImp,
                                          domainDim,
                                          RangeFieldImp,
                                          rangeDim,
                                          rangeDimCols>
{
  typedef AffineFunction<AffineEntityImp, AffineDomainFieldImp, affineDomainDim, RangeFieldImp, rangeDim, rangeDimCols>
      AffineFunctionType;
  typedef typename AffineFunctionType::MatrixType MatrixType;
  typedef typename AffineFunctionType::RangeType RangeType;
  typedef FunctionCheckerboardFunction<AffineFunctionType,
                                       EntityImp,
                                       DomainFieldImp,
                                       domainDim,
                                       RangeFieldImp,
                                       rangeDim,
                                       rangeDimCols>
      FunctionCheckerboardType;

public:
  static std::unique_ptr<FunctionCheckerboardType>
  create(const Common::Configuration config = CheckerboardConfigProvider<AffineFunctionType>::default_config(),
         const std::string sub_name = FunctionCheckerboardType::static_id())
  {
    // get correct config
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = CheckerboardConfigProvider<AffineFunctionType>::default_config();
    // calculate number of values and get values
    auto num_elements =
        cfg.get("num_elements", default_cfg.get<Common::FieldVector<size_t, domainDim>>("num_elements"), domainDim);
    size_t num_values = 1;
    for (size_t ii = 0; ii < num_elements.size(); ++ii)
      num_values *= num_elements[ii];
    std::vector<AffineFunctionType> values;
    for (size_t ii = 0; ii < num_values; ++ii) {
      const MatrixType A_ii = cfg.get<MatrixType>("A." + Common::to_string(ii), rangeDim, affineDomainDim);
      const RangeType b_ii = cfg.get<RangeType>("b." + Common::to_string(ii), rangeDim);
      values.emplace_back(AffineFunctionType(A_ii, b_ii));
    }
    // create
    return Common::make_unique<FunctionCheckerboardType>(
        cfg.get("lower_left", default_cfg.get<Common::FieldVector<DomainFieldImp, domainDim>>("lower_left"), domainDim),
        cfg.get(
            "upper_right", default_cfg.get<Common::FieldVector<DomainFieldImp, domainDim>>("upper_right"), domainDim),
        num_elements,
        values,
        cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)
}; // class FunctionCheckerboardFactory< Affine, ... >


/** TODO: replace Checkerboard completely (and rename FunctionCheckerboard to Checkerboard)?
 * FunctionCheckerboard< Constant, ... > should give the same results as Checkerboard. What about performance? */
template <class LocalizableFunctionImp,
          class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols>
class FunctionCheckerboardFunction
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  static_assert(is_localizable_function<LocalizableFunctionImp>::value,
                "LocalizableFunctionImp needs to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(domainDim <= 3, "Not implemented for dimDomain > 3 (see find_subdomain method)!");
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      BaseType;
  typedef FunctionCheckerboardFunction<LocalizableFunctionImp,
                                       EntityImp,
                                       DomainFieldImp,
                                       domainDim,
                                       RangeFieldImp,
                                       rangeDim,
                                       rangeDimCols>
      ThisType;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;
  using typename BaseType::EntityType;
  using typename BaseType::LocalfunctionType;

  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;
  using typename BaseType::RangeFieldType;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;

  static const bool available = true;

  typedef LocalizableFunctionImp LocalizableFunctionType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".functioncheckerboard";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    return CheckerboardConfigProvider<LocalizableFunctionType>::default_config(sub_name);
  }

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    return FunctionCheckerboardFunctionFactory<LocalizableFunctionImp,
                                               EntityImp,
                                               DomainFieldImp,
                                               domainDim,
                                               RangeFieldImp,
                                               rangeDim,
                                               rangeDimCols>::create(config, sub_name);
  } // ... create(...)

  FunctionCheckerboardFunction(const Common::FieldVector<DomainFieldType, dimDomain>& lowerLeft,
                               const Common::FieldVector<DomainFieldType, dimDomain>& upperRight,
                               const Common::FieldVector<size_t, dimDomain>& numElements,
                               const std::vector<LocalizableFunctionType>& values,
                               const std::string nm = static_id())
    : lowerLeft_(new Common::FieldVector<DomainFieldType, dimDomain>(lowerLeft))
    , upperRight_(new Common::FieldVector<DomainFieldType, dimDomain>(upperRight))
    , numElements_(new Common::FieldVector<size_t, dimDomain>(numElements))
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

  FunctionCheckerboardFunction(const ThisType& other) = default;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  virtual std::string type() const override
  {
    return BaseType::static_id() + ".functioncheckerboard";
  }

  virtual std::string name() const override
  {
    return name_;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override
  {
    return LocalFunctionChooser<LocalfunctionType,
                                LocalizableFunctionType,
                                std::is_same<EntityType, typename LocalizableFunctionType::EntityType>::value>::
        local_function(entity, localizable_function(entity));
  } // ... local_function(...)

  const LocalizableFunctionType& localizable_function(const EntityType& entity) const
  {
    const size_t subdomain = find_subdomain(entity);
    // return the component that belongs to the subdomain
    return *(values_[subdomain]);
  } // ... localizable_function(...)

private:
  template <class LocalfunctionType, class LocalizableFunctionType, bool entity_types_match = true>
  struct LocalFunctionChooser
  {
    static std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity,
                                                             const LocalizableFunctionType& localizable_function_entity)
    {
      return localizable_function_entity.local_function(entity);
    }
  };

  template <class LocalfunctionType, class LocalizableFunctionType>
  struct LocalFunctionChooser<LocalfunctionType, LocalizableFunctionType, false>
  {
    static std::unique_ptr<LocalfunctionType>
    local_function(const EntityType& /*entity*/, const LocalizableFunctionType& /*localizable_function_entity*/)
    {
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong,
                 "Checkerboard cannot provide a local function because LocalizableFunctionImp and Checkerboard are "
                 "defined on different grids. Use localizable_function(entity) instead!");
      return std::unique_ptr<LocalfunctionType>();
    }
  };

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

  std::shared_ptr<const Common::FieldVector<DomainFieldType, dimDomain>> lowerLeft_;
  std::shared_ptr<const Common::FieldVector<DomainFieldType, dimDomain>> upperRight_;
  std::shared_ptr<const Common::FieldVector<size_t, dimDomain>> numElements_;
  std::vector<std::shared_ptr<const LocalizableFunctionType>> values_;
  std::string name_;
}; // class FunctionCheckerboardFunction

#endif // HAVE_DUNE_XT_LA

} // namespace Functions
} // namespace XT
} // namespace Dune


#include "checkerboard.lib.hh"


#endif // DUNE_XT_FUNCTIONS_CHECKERBOARD_HH
