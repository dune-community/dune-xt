// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTIONS_FLATTOP_HH
#define DUNE_STUFF_FUNCTIONS_FLATTOP_HH

#include <dune/common/static_assert.hh>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/functions/interfaces.hh>

namespace Dune {
namespace Stuff {
namespace Functions {


/**
 * Based on: Brenner, S. C. and Davis, C. B. and Sung, L.-y.
 *           A partition of unity method for the displacement obstacle problem of clamped Kirchhoff plates
 *           http://dx.doi.org/10.1016/j.cam.2013.09.033
 *           Subsection 2.1.1
 */
template <class E, class D, int d, class R, int r, int rC = 1>
class FlatTop
{
  static_assert(AlwaysFalse<E>::value, "Not available for these dimensions!");
};


template <class E, class D, int d, class R>
class FlatTop<E, D, d, R, 1, 1> : public GlobalFunctionInterface<E, D, d, R, 1, 1>
{
  typedef GlobalFunctionInterface<E, D, d, R, 1, 1> BaseType;
  typedef FlatTop<E, D, d, R, 1, 1> ThisType;

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimDomain = BaseType::dimDomain;

  static std::string static_id()
  {
    return BaseType::static_id() + ".flattop";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["lower_left"]     = "[0.0 0.0 0.0]";
    config["upper_right"]    = "[1.0 1.0 1.0]";
    config["boundary_layer"] = "[1e-1 1e-1 1e-1]";
    config["value"]          = "1.0";
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
    const Common::Configuration cfg         = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    return Common::make_unique<ThisType>(cfg.get("lower_left", default_cfg.get<DomainType>("lower_left")),
                                         cfg.get("upper_right", default_cfg.get<DomainType>("upper_right")),
                                         cfg.get("boundary_layer", default_cfg.get<DomainType>("boundary_layer")),
                                         cfg.get("value", default_cfg.get<RangeType>("value")),
                                         cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)

  FlatTop(const DomainType& lower_left = default_config().get<DomainType>("lower_left"),
          const DomainType& upper_right    = default_config().get<DomainType>("upper_right"),
          const DomainType& boundary_layer = default_config().get<DomainType>("boundary_layer"),
          const RangeType& value           = default_config().get<RangeType>("value"),
          const std::string name = default_config().get<std::string>("name"))
    : lower_left_(lower_left)
    , upper_right_(upper_right)
    , boundary_layer_(boundary_layer)
    , value_(value)
    , name_(name)
  {
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      if (!(upper_right_[dd] > lower_left_[dd]))
        DUNE_THROW(
            Exceptions::wrong_input_given,
            "lower_left_[" << dd << "] = " << lower_left_[dd] << ", upper_right_[" << dd << "] = " << upper_right_);
      if (!(boundary_layer_[dd] > 0))
        DUNE_THROW(Exceptions::wrong_input_given, "boundary_layer_[" << dd << "] = " << boundary_layer_[dd]);
    }
  } // ... FlatTop(...)

  FlatTop(const ThisType& other)
    : lower_left_(other.lower_left_)
    , upper_right_(other.upper_right_)
    , boundary_layer_(other.boundary_layer_)
    , value_(other.value_)
    , name_(other.name_)
  {
  }

  ThisType& operator=(const ThisType& other) = delete;

  virtual ~FlatTop()
  {
  }

  virtual ThisType* copy() const /*DS_OVERRIDE DS_FINAL*/
  {
    return new ThisType(*this);
  }

  virtual std::string type() const /*DS_OVERRIDE DS_FINAL*/
  {
    return BaseType::static_id() + ".flattop";
  }

  virtual std::string name() const /*DS_OVERRIDE DS_FINAL*/
  {
    return name_;
  }

  virtual size_t order() const
  {
    // could be 3 but I am to lazy
    return 4;
  } // ... order(...)

  virtual void evaluate(const DomainType& xx, RangeType& ret) const
  {
    ret = value_;
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      const auto& left  = lower_left_[dd];
      const auto& right = upper_right_[dd];
      const auto& point = xx[dd];
      const auto& delta = boundary_layer_[dd];
      if (point < left - delta) {
        // outside
        ret[0] = 0.0;
        break;
      } else if (point < left) {
        // left boundary layer
        ret[0] *= phi_left((point - left) / delta);
      } else if (point < right) {
        // inside
        // do nothing (keep value)
      } else if (point < right + delta) {
        // right boundary layer
        ret[0] *= phi_right((point - right) / delta);
      } else {
        // outside
        ret[0] = 0.0;
        break;
      }
    }
  } // ... evaluate(...)

private:
  RangeFieldType phi_left(const RangeFieldType& point) const
  {
    assert(!(point < -1.0));
    assert(!(point > 0.0));
    return std::pow(1.0 + point, 2) * (1.0 - 2.0 * point);
  } // ... phi_left(...)

  RangeFieldType phi_right(const RangeFieldType& point) const
  {
    assert(!(point < 0.0));
    assert(!(point > 1.0));
    return std::pow(1.0 - point, 2) * (1.0 + 2.0 * point);
  } // ... phi_right(...)

  const DomainType lower_left_;
  const DomainType upper_right_;
  const DomainType boundary_layer_;
  const RangeType value_;
  const std::string name_;
}; // class FlatTop< ..., 1, 1 >


} // namespace Functions
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTIONS_FLATTOP_HH
