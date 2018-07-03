// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2012 - 2017)
//   Rene Milk       (2013 - 2018)
//   Tobias Leibner  (2014 - 2017)

#ifndef DUNE_XT_FUNCTIONS_CHECKERBOARD_HH
#define DUNE_XT_FUNCTIONS_CHECKERBOARD_HH

#include <dune/xt/common/configuration.hh>

#include <dune/xt/functions/interfaces.hh>

namespace Dune {
namespace XT {
namespace Functions {

/**
 * Note: This function does not allow for functions on the subdomains anymore. Only constant values are possible.
 */
template <class E, size_t r = 1, size_t rC = 1, class R = double>
class CheckerboardFunction : public GridFunctionInterface<E, r, rC, R>
{
  using BaseType = GridFunctionInterface<E, r, rC, R>;
  using ThisType = CheckerboardFunction<E, r, rC, R>;
  using BaseType::dimDomain;
  static_assert(dimDomain <= 3, "Not implemented for dimDomain > 3 (see find_subdomain method)!");

  class LocalCheckerboardFunction : public ElementFunctionInterface<E, r, rC, R>
  {
    typedef ElementFunctionInterface<E, r, rC, R> InterfaceType;

  public:
    using typename InterfaceType::ElementType;
    using typename InterfaceType::DomainType;
    using typename InterfaceType::RangeType;
    using typename InterfaceType::DerivativeRangeType;
    using GeometryType = typename ElementType::Geometry;

    LocalCheckerboardFunction(const ElementType& element,
                              const DomainType& lower_left,
                              const DomainType& upper_right,
                              const FieldVector<size_t, dimDomain>& num_elements,
                              const std::vector<RangeType>& values)
      : InterfaceType(element)
      , lower_left_(lower_left)
      , upper_right_(upper_right)
      , num_elements_(num_elements)
      , values_(values)
    {
      post_bind(element);
    }

    LocalCheckerboardFunction(const DomainType& lower_left,
                              const DomainType& upper_right,
                              const FieldVector<size_t, dimDomain>& num_elements,
                              std::vector<RangeType>& values)
      : InterfaceType()
      , lower_left_(lower_left)
      , upper_right_(upper_right)
      , num_elements_(num_elements)
      , values_(values)
    {
    }

    void post_bind(const ElementType& element) override final
    {
      current_value_ = 0;
      if (is_in_checkerboard(element)) {
        const size_t subdomain = find_subdomain(element);
        current_value_ = values_[subdomain];
      }
    } // ... post_bind(...)

    int order(const Common::Parameter& /*param*/ = {}) const override final
    {
      return 0;
    }

    RangeType evaluate(const DomainType& point_in_reference_element,
                       const Common::Parameter& /*param*/ = {}) const override final
    {
      this->assert_inside_reference_element(point_in_reference_element);
      return current_value_;
    }

    DerivativeRangeType jacobian(const DomainType& point_in_reference_element,
                                 const Common::Parameter& /*param*/ = {}) const override final
    {
      this->assert_inside_reference_element(point_in_reference_element);
      return DerivativeRangeType();
    }

  private:
    bool is_in_checkerboard(const ElementType& element) const
    {
      const auto center = element.geometry().center();
      if (Common::FloatCmp::le(lower_left_, center) && Common::FloatCmp::lt(center, upper_right_))
        return 1;
      else
        return 0;
    }

    size_t find_subdomain(const ElementType& element) const
    {
      // decide to which subdomain the center of the element belongs to
      const auto center = element.geometry().center();
      std::vector<size_t> which_partition(dimDomain, 0);
      const auto& ll = lower_left_;
      const auto& ur = upper_right_;
      const auto& ne = num_elements_;
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        // for points that are on upper_right_[d], this selects one partition too much
        // so we need to cap this
        which_partition[dd] =
            std::min(size_t(std::floor(ne[dd] * ((center[dd] - ll[dd]) / (ur[dd] - ll[dd])))), ne[dd] - 1);
      }
      size_t subdomain = 0;
      if (dimDomain == 1)
        subdomain = which_partition[0];
      else if (dimDomain == 2)
        subdomain = which_partition[0] + which_partition[1] * ne[0];
      else
        subdomain = which_partition[0] + which_partition[1] * ne[0] + which_partition[2] * ne[1] * ne[0];
      return subdomain;
    } // ... find_subdomain(...)


    const DomainType lower_left_;
    const DomainType upper_right_;
    const FieldVector<size_t, dimDomain> num_elements_;
    const std::vector<RangeType>& values_;
    RangeType current_value_;
  }; // class LocalCheckerboardFunction

public:
  using typename BaseType::ElementType;
  using typename BaseType::LocalFunctionType;

  using RangeType = typename LocalFunctionType::RangeType;
  using DomainType = typename LocalFunctionType::DomainType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".checkerboard";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["type"] = static_id();
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
                                          const std::string sub_name = "")
  {
    // get correct config
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration def_cfg = default_config();
    // calculate number of values and get values
    auto num_elements =
        cfg.get("num_elements", def_cfg.get<Common::FieldVector<size_t, dimDomain>>("num_elements"), dimDomain);
    size_t num_values = 1;
    for (size_t ii = 0; ii < num_elements.size(); ++ii)
      num_values *= num_elements[ii];
    std::vector<RangeType> values;
    auto values_range = cfg.get("values", def_cfg.get<std::vector<R>>("values"), num_values);
    for (size_t ii = 0; ii < values_range.size(); ++ii)
      values.emplace_back(values_range[ii]);
    // create
    return Common::make_unique<ThisType>(cfg.get("lower_left", def_cfg.get<DomainType>("lower_left"), dimDomain),
                                         cfg.get("upper_right", def_cfg.get<DomainType>("upper_right"), dimDomain),
                                         std::move(num_elements),
                                         std::move(values),
                                         cfg.get("name", def_cfg.get<std::string>("name")));
  } // ... create(...)

  CheckerboardFunction(const DomainType& lower_left,
                       const DomainType& upper_right,
                       const FieldVector<size_t, dimDomain>& num_elements,
                       const std::vector<RangeType>& values,
                       const std::string nm = "checkerboard")
    : lower_left_(lower_left)
    , upper_right_(upper_right)
    , num_elements_(num_elements)
    , values_(new std::vector<RangeType>(values))
    , name_(nm)
  {
#ifndef NDEBUG
    // checks
    size_t total_subdomains = 1;
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      const auto& ll = (lower_left_)[dd];
      const auto& ur = (upper_right_)[dd];
      const auto& ne = (num_elements_)[dd];
      if (!(ll < ur))
        DUNE_THROW(Dune::RangeError, "lower_left has to be elementwise smaller than upper_right!");
      total_subdomains *= ne;
    }
    if (values_->size() < total_subdomains)
      DUNE_THROW(Dune::RangeError,
                 "values too small (is " << values_->size() << ", should be " << total_subdomains << ")");
#endif
  } // CheckerboardFunction(...)

  CheckerboardFunction(const ThisType& other) = default;
  CheckerboardFunction(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  std::string type() const override
  {
    return BaseType::static_id() + ".checkerboard";
  }

  std::string name() const override
  {
    return name_;
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<LocalCheckerboardFunction>(lower_left_, upper_right_, num_elements_, *values_);
  }

  std::unique_ptr<LocalFunctionType> local_function(const ElementType& element) const override final
  {
    return std::make_unique<LocalCheckerboardFunction>(element, lower_left_, upper_right_, num_elements_, *values_);
  }

  size_t subdomain(const ElementType& element) const
  {
    return find_subdomain(element);
  }

  size_t subdomains() const
  {
    return values_->size();
  }

  const std::vector<std::shared_ptr<const RangeType>>& values() const
  {
    return values_;
  }

private:
  size_t find_subdomain(const ElementType& element) const
  {
    // decide on the subdomain the center of the element belongs to
    const auto center = element.geometry().center();
    std::vector<size_t> which_partition(dimDomain, 0);
    const auto& ll = lower_left_;
    const auto& ur = upper_right_;
    const auto& ne = num_elements_;
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      // for points that are on upper_right_[d], this selects one partition too much
      // so we need to cap this
      which_partition[dd] =
          std::min(size_t(std::floor(ne[dd] * ((center[dd] - ll[dd]) / (ur[dd] - ll[dd])))), ne[dd] - 1);
    }
    size_t subdomain = 0;
    if (dimDomain == 1)
      subdomain = which_partition[0];
    else if (dimDomain == 2)
      subdomain = which_partition[0] + which_partition[1] * ne[0];
    else
      subdomain = which_partition[0] + which_partition[1] * ne[0] + which_partition[2] * ne[1] * ne[0];
    return subdomain;
  } // ... find_subdomain(...)

  const DomainType lower_left_;
  const DomainType upper_right_;
  const FieldVector<size_t, dimDomain> num_elements_;
  std::shared_ptr<std::vector<RangeType>> values_;
  std::string name_;
}; // class CheckerboardFunction

} // namespace Functions
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_FUNCTIONS_CHECKERBOARD_HH
