// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2018)
//   René Fritze     (2014 - 2016, 2018)
//   Tim Keil        (2018)
//   Tobias Leibner  (2015, 2017)

#ifndef DUNE_XT_FUNCTIONS_FLATTOP_HH
#define DUNE_XT_FUNCTIONS_FLATTOP_HH

#include <dune/xt/common/configuration.hh>

#include <dune/xt/functions/interfaces/function.hh>

namespace Dune {
namespace XT {
namespace Functions {


/**
 * Based on: Brenner, S. C. and Davis, C. B. and Sung, L.-y.
 *           A partition of unity method for the displacement obstacle problem of clamped Kirchhoff plates
 *           http://dx.doi.org/10.1016/j.cam.2013.09.033
 *           Subsection 2.1.1
 */
template <size_t d, size_t r, size_t rC = 1, class R = double>
class FlatTopFunction : public FunctionInterface<d, r, rC, R>
{
  FlatTopFunction()
  {
    static_assert(AlwaysFalse<R>::value, "Not available for these dimensions!");
  }
};

template <size_t d, class R>
class FlatTopFunction<d, 1, 1, R> : public FunctionInterface<d, 1, 1, R>
{
  using BaseType = FunctionInterface<d, 1, 1, R>;
  using ThisType = FlatTopFunction<d, 1, 1, R>;

public:
  using DomainFieldType = typename BaseType::DomainFieldType;
  using RangeFieldType = typename BaseType::RangeFieldType;
  using DomainType = typename BaseType::DomainType;
  using RangeReturnType = typename BaseType::RangeReturnType;

  static const size_t domain_dim = BaseType::domain_dim;
  static const size_t range_dim = BaseType::range_dim;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".flattop";
  }

  static Common::Configuration defaults()
  {
    Common::Configuration config;
    config["lower_left"] = "[0.0 0.0 0.0]";
    config["upper_right"] = "[1.0 1.0 1.0]";
    config["boundary_layer"] = "[1e-1 1e-1 1e-1]";
    config["value"] = "1.0";
    config["name"] = static_id();
    return config;
  } // ... defaults(...)

  FlatTopFunction(const DomainType& lower_left,
                  const DomainType& upper_right,
                  const DomainType& boundary_layer,
                  const RangeReturnType& value = defaults().template get<RangeReturnType>("value"),
                  const std::string name_in = defaults().template get<std::string>("name"))
    : lower_left_(lower_left)
    , upper_right_(upper_right)
    , boundary_layer_(boundary_layer)
    , value_(value)
    , name_(name_in)
  {
    check_input();
  }

  FlatTopFunction(const ThisType& other) = default;

  ThisType& operator=(const ThisType& other) = delete;

  virtual ~FlatTopFunction() {}

  std::string name() const override final
  {
    return name_;
  }

  int order(const XT::Common::Parameter& /*param*/ = {}) const override
  {
    return 3 * domain_dim;
  }

  RangeReturnType evaluate(const DomainType& point_in_reference_element,
                           const Common::Parameter& /*param*/ = {}) const override final
  {
    RangeReturnType ret = value_;
    for (size_t dd = 0; dd < domain_dim; ++dd) {
      const auto& left = lower_left_[dd];
      const auto& right = upper_right_[dd];
      const auto& point = point_in_reference_element[dd];
      const auto& delta = boundary_layer_[dd];
      if (point < left - delta) {
        // outside
        ret[0] = 0.0;
        break;
      } else if (point < left + delta) {
        // left boundary layer
        ret[0] *= phi_left((point - (left + delta)) / (2.0 * delta));
      } else if (point < right - delta) {
        // inside
        // do nothing (keep value)
      } else if (point < right + delta) {
        // right boundary layer
        ret[0] *= phi_right((point - (right - delta)) / (2.0 * delta));
      } else {
        // outside
        ret[0] = 0.0;
        break;
      }
    }
    return ret;
  } // ... evaluate(...)

private:
  void check_input() const
  {
    if (!(Common::FloatCmp::gt(upper_right_, lower_left_)))
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "upper_right has to be greater than lower_left!\n"
                     << "lower_left = [" << lower_left_ << "]\n"
                     << "upper_right = [" << upper_right_ << "]");
    if (!(Common::FloatCmp::gt(boundary_layer_, DomainType(0))))
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "boundary_layer has to be strictly positive!\n"
                     << "boundary_layer = [" << boundary_layer_ << "]");
    if (Common::FloatCmp::gt(boundary_layer_ * 2.0, upper_right_ - lower_left_))
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "boundary_layer has to be thin enough!\n"
                 "2*boundary_layer = ["
                     << boundary_layer_ * 2.0 << "]\n"
                     << "upper_right - lower_left = [" << upper_right_ - lower_left_ << "]");
  } // .. check_input(...)

  RangeFieldType phi_left(const RangeFieldType& point) const
  {
    if (point < -1.0)
      return 0.0;
    else if (point > 0.0)
      return 1.0;
    else
      return std::pow(1.0 + point, 2) * (1.0 - 2.0 * point);
  } // ... phi_left(...)

  RangeFieldType phi_right(const RangeFieldType& point) const
  {
    if (point < 0.0)
      return 1.0;
    else if (point > 1.0)
      return 0.0;
    else
      return std::pow(1.0 - point, 2) * (1.0 + 2.0 * point);
  } // ... phi_right(...)

  const DomainType lower_left_;
  const DomainType upper_right_;
  const DomainType boundary_layer_;
  const RangeReturnType value_;
  const std::string name_;
}; // class FlatTopFunction< ..., 1, 1 >


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_FLATTOP_HH
