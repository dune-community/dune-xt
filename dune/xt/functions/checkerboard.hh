// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2012 - 2018)
//   René Fritze     (2013 - 2018)
//   TiKeil          (2018)
//   Tobias Leibner  (2014 - 2017)

#ifndef DUNE_XT_FUNCTIONS_CHECKERBOARD_HH
#define DUNE_XT_FUNCTIONS_CHECKERBOARD_HH

#include <dune/xt/common/configuration.hh>

#include <dune/xt/functions/interfaces/grid-function.hh>

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
  using BaseType::domain_dim;
  static_assert(domain_dim <= 3, "Not implemented for domain_dim > 3 (see find_subdomain method)!");

  class LocalCheckerboardFunction : public ElementFunctionInterface<E, r, rC, R>
  {
    using InterfaceType = ElementFunctionInterface<E, r, rC, R>;

  public:
    using typename InterfaceType::ElementType;
    using typename InterfaceType::DomainType;
    using typename InterfaceType::RangeType;
    using typename InterfaceType::RangeReturnType;
    using typename InterfaceType::DerivativeRangeType;
    using typename InterfaceType::DerivativeRangeReturnType;
    using GeometryType = typename ElementType::Geometry;

    LocalCheckerboardFunction(const DomainType& lower_left,
                              const DomainType& upper_right,
                              const FieldVector<size_t, domain_dim>& num_elements,
                              std::vector<RangeType>& values)
      : InterfaceType()
      , lower_left_(lower_left)
      , upper_right_(upper_right)
      , num_elements_(num_elements)
      , values_(values)
    {
    }

  protected:
    void post_bind(const ElementType& element) override final
    {
      current_value_ = 0;
      if (is_in_checkerboard(element)) {
        const size_t subdomain = find_subdomain(element);
        current_value_ = values_[subdomain];
      }
    }

  public:
    int order(const Common::Parameter& /*param*/ = {}) const override final
    {
      return 0;
    }

    RangeReturnType evaluate(const DomainType& point_in_reference_element,
                             const Common::Parameter& /*param*/ = {}) const override final
    {
      this->assert_inside_reference_element(point_in_reference_element);
      return current_value_;
    }

    DerivativeRangeReturnType jacobian(const DomainType& point_in_reference_element,
                                       const Common::Parameter& /*param*/ = {}) const override final
    {
      this->assert_inside_reference_element(point_in_reference_element);
      return DerivativeRangeReturnType();
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
      std::vector<size_t> which_partition(domain_dim, 0);
      const auto& ll = lower_left_;
      const auto& ur = upper_right_;
      const auto& ne = num_elements_;
      for (size_t dd = 0; dd < domain_dim; ++dd) {
        // for points that are on upper_right_[d], this selects one partition too much
        // so we need to cap this
        which_partition[dd] =
            std::min(size_t(std::floor(ne[dd] * ((center[dd] - ll[dd]) / (ur[dd] - ll[dd])))), ne[dd] - 1);
      }
      size_t subdomain = 0;
      if (domain_dim == 1)
        subdomain = which_partition[0];
      else if (domain_dim == 2)
        subdomain = which_partition[0] + which_partition[1] * ne[0];
      else
        subdomain = which_partition[0] + which_partition[1] * ne[0] + which_partition[2] * ne[1] * ne[0];
      return subdomain;
    } // ... find_subdomain(...)


    const DomainType lower_left_;
    const DomainType upper_right_;
    const FieldVector<size_t, domain_dim> num_elements_;
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

  static Common::Configuration defaults()
  {
    Common::Configuration config;
    config["lower_left"] = "[0.0 0.0 0.0]";
    config["upper_right"] = "[1.0 1.0 1.0]";
    config["num_elements"] = "[2 2 2]";
    config["values"] = "[1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0]";
    config["name"] = static_id();
    return config;
  } // ... defaults(...)

  CheckerboardFunction(const DomainType& lower_left,
                       const DomainType& upper_right,
                       const FieldVector<size_t, domain_dim>& num_elements,
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
    for (size_t dd = 0; dd < domain_dim; ++dd) {
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

  std::string name() const override
  {
    return name_;
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<LocalCheckerboardFunction>(lower_left_, upper_right_, num_elements_, *values_);
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
    std::vector<size_t> which_partition(domain_dim, 0);
    const auto& ll = lower_left_;
    const auto& ur = upper_right_;
    const auto& ne = num_elements_;
    for (size_t dd = 0; dd < domain_dim; ++dd) {
      // for points that are on upper_right_[d], this selects one partition too much
      // so we need to cap this
      which_partition[dd] =
          std::min(size_t(std::floor(ne[dd] * ((center[dd] - ll[dd]) / (ur[dd] - ll[dd])))), ne[dd] - 1);
    }
    size_t subdomain = 0;
    if (domain_dim == 1)
      subdomain = which_partition[0];
    else if (domain_dim == 2)
      subdomain = which_partition[0] + which_partition[1] * ne[0];
    else
      subdomain = which_partition[0] + which_partition[1] * ne[0] + which_partition[2] * ne[1] * ne[0];
    return subdomain;
  } // ... find_subdomain(...)

  const DomainType lower_left_;
  const DomainType upper_right_;
  const FieldVector<size_t, domain_dim> num_elements_;
  std::shared_ptr<std::vector<RangeType>> values_;
  std::string name_;
}; // class CheckerboardFunction


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_CHECKERBOARD_HH
