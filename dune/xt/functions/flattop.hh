// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014 - 2016)
//   Tobias Leibner  (2015)

#ifndef DUNE_XT_FUNCTIONS_FLATTOP_HH
#define DUNE_XT_FUNCTIONS_FLATTOP_HH

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/type_traits.hh>

#include <dune/xt/functions/interfaces.hh>

namespace Dune {
namespace XT {
namespace Functions {


/**
 * Based on: Brenner, S. C. and Davis, C. B. and Sung, L.-y.
 *           A partition of unity method for the displacement obstacle problem of clamped Kirchhoff plates
 *           http://dx.doi.org/10.1016/j.cam.2013.09.033
 *           Subsection 2.1.1
 */
template <class E, class D, size_t d, class R, size_t r, size_t rC = 1>
class FlatTopFunction : public LocalizableFunctionInterface<E, D, d, R, r, rC>
{
  FlatTopFunction()
  {
    static_assert(AlwaysFalse<E>::value, "Not available for these dimensions!");
  }
};

template <class E, class D, size_t d, class R>
class FlatTopFunction<E, D, d, R, 1, 1> : public GlobalFunctionInterface<E, D, d, R, 1, 1>
{
  typedef GlobalFunctionInterface<E, D, d, R, 1, 1> BaseType;
  typedef FlatTopFunction<E, D, d, R, 1, 1> ThisType;

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const size_t dimRange = BaseType::dimRange;
  typedef typename BaseType::RangeType RangeType;

  typedef Common::FieldVector<DomainFieldType, dimDomain> StuffDomainType;
  typedef Common::FieldVector<RangeFieldType, dimRange> StuffRangeType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".flattop";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["lower_left"] = "[0.0 0.0 0.0]";
    config["upper_right"] = "[1.0 1.0 1.0]";
    config["boundary_layer"] = "[1e-1 1e-1 1e-1]";
    config["value"] = "1.0";
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
    return Common::make_unique<ThisType>(cfg.get("lower_left", default_cfg.get<DomainType>("lower_left")),
                                         cfg.get("upper_right", default_cfg.get<DomainType>("upper_right")),
                                         cfg.get("boundary_layer", default_cfg.get<DomainType>("boundary_layer")),
                                         cfg.get("value", default_cfg.get<RangeType>("value")),
                                         cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)

  FlatTopFunction(const StuffDomainType& lower_left,
                  const StuffDomainType& upper_right,
                  const StuffDomainType& boundary_layer,
                  const StuffRangeType& value = default_config().template get<StuffRangeType>("value"),
                  const std::string name_in = default_config().template get<std::string>("name"))
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

  virtual ~FlatTopFunction()
  {
  }

  virtual std::string type() const override final
  {
    return BaseType::static_id() + ".flattop";
  }

  virtual std::string name() const override final
  {
    return name_;
  }

  virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override
  {
    return 3 * dimDomain;
  }

  virtual void evaluate(const DomainType& xx, RangeType& ret, const Common::Parameter& /*mu*/ = {}) const override
  {
    ret = value_;
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      const auto& left = lower_left_[dd];
      const auto& right = upper_right_[dd];
      const auto& point = xx[dd];
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
  } // ... evaluate(...)

private:
  void check_input() const
  {
    if (!(Common::FloatCmp::gt(upper_right_, lower_left_)))
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "upper_right has to be greater than lower_left!\n"
                     << "lower_left = ["
                     << lower_left_
                     << "]\n"
                     << "upper_right = ["
                     << upper_right_
                     << "]");
    if (!(Common::FloatCmp::gt(boundary_layer_, StuffDomainType(0))))
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "boundary_layer has to be strictly positive!\n"
                     << "boundary_layer = ["
                     << boundary_layer_
                     << "]");
    if (Common::FloatCmp::gt(boundary_layer_ * 2.0, upper_right_ - lower_left_))
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "boundary_layer has to be thin enough!\n"
                 "2*boundary_layer = ["
                     << boundary_layer_ * 2.0
                     << "]\n"
                     << "upper_right - lower_left = ["
                     << upper_right_ - lower_left_
                     << "]");
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

  const StuffDomainType lower_left_;
  const StuffDomainType upper_right_;
  const StuffDomainType boundary_layer_;
  const StuffRangeType value_;
  const std::string name_;
}; // class FlatTopFunction< ..., 1, 1 >


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_FLATTOP_HH
