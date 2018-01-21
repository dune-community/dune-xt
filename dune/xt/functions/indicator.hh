// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2015 - 2016, 2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_FUNCTIONS_INDICATOR_HH
#define DUNE_XT_FUNCTIONS_INDICATOR_HH

#include <memory>
#include <utility>
#include <vector>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace XT {
namespace Functions {

template <class E, class D, size_t d, class R, size_t r, size_t rC = 1>
class IndicatorFunction : public LocalizableFunctionInterface<E, D, d, R, r, rC>
{
  IndicatorFunction()
  {
    static_assert(AlwaysFalse<E>::value, "Not available for these dimensions!");
  }
};

template <class E, class D, size_t d, class R>
class IndicatorFunction<E, D, d, R, 1> : public LocalizableFunctionInterface<E, D, d, R, 1>
{
  typedef LocalizableFunctionInterface<E, D, d, R, 1> BaseType;
  typedef IndicatorFunction<E, D, d, R, 1> ThisType;

  class Localfunction : public LocalfunctionInterface<E, D, d, R, 1>
  {
    typedef LocalfunctionInterface<E, D, d, R, 1> InterfaceType;

  public:
    using typename InterfaceType::EntityType;
    using typename InterfaceType::DomainType;
    using typename InterfaceType::RangeType;
    using typename InterfaceType::JacobianRangeType;

    Localfunction(const EntityType& entity, const RangeType& value)
      : InterfaceType(entity)
      , value_(value)
    {
    }

    virtual size_t order(const Common::Parameter& /*mu*/ = {}) const override final
    {
      return 0;
    }

    virtual void
    evaluate(const DomainType& xx, RangeType& ret, const Common::Parameter& /*mu*/ = {}) const override final
    {
      assert(this->is_a_valid_point(xx));
      ret = value_;
    }

    virtual void
    jacobian(const DomainType& xx, JacobianRangeType& ret, const Common::Parameter& /*mu*/ = {}) const override final
    {
      assert(this->is_a_valid_point(xx));
      ret *= 0.0;
    }

  private:
    const RangeType value_;
  }; // class Localfunction

public:
  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::LocalfunctionType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".indicator";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration cfg;
    cfg["name"] = static_id();
    if (d == 1)
      cfg["0.domain"] = "[0.25 0.75]";
    else if (d == 2)
      cfg["0.domain"] = "[0.25 0.75; 0.25 0.75]";
    else if (d == 3)
      cfg["0.domain"] = "[0.25 0.75; 0.25 0.75; 0.25 0.75]";
    else
      DUNE_THROW(NotImplemented, "Indeed!");
    cfg["0.value"] = "1";
    if (sub_name.empty())
      return cfg;
    else {
      Common::Configuration tmp;
      tmp.add(cfg, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration def_cfg = default_config();
    std::vector<std::tuple<DomainType, DomainType, RangeFieldType>> values;
    DomainType tmp_lower;
    DomainType tmp_upper;
    size_t cc = 0;
    while (cfg.has_sub(Common::to_string(cc))) {
      const Common::Configuration local_cfg = cfg.sub(Common::to_string(cc));
      if (local_cfg.has_key("domain") && local_cfg.has_key("value")) {
        auto domains = local_cfg.get<FieldMatrix<DomainFieldType, d, 2>>("domain");
        for (size_t dd = 0; dd < d; ++dd) {
          tmp_lower[dd] = domains[dd][0];
          tmp_upper[dd] = domains[dd][1];
        }
        auto val = local_cfg.get<RangeFieldType>("value");
        values.emplace_back(tmp_lower, tmp_upper, val);
      } else
        break;
      ++cc;
    }
    return Common::make_unique<ThisType>(values, cfg.get("name", def_cfg.get<std::string>("name")));
  } // ... create(...)

  IndicatorFunction(const std::vector<std::tuple<DomainType, DomainType, R>>& values,
                    const std::string name_in = "indicator")
    : values_(values)
    , name_(name_in)
  {
  }

  IndicatorFunction(
      const std::vector<std::pair<std::pair<Common::FieldVector<D, d>, Common::FieldVector<D, d>>, R>>& values,
      const std::string name_in = "indicator")
    : values_(convert(values))
    , name_(name_in)
  {
  }

  virtual ~IndicatorFunction()
  {
  }

  virtual std::string name() const override final
  {
    return name_;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override final
  {
    const auto center = entity.geometry().center();
    for (const auto& element : values_)
      if (Common::FloatCmp::le(std::get<0>(element), center) && Common::FloatCmp::lt(center, std::get<1>(element)))
        return Common::make_unique<Localfunction>(entity, std::get<2>(element));
    return Common::make_unique<Localfunction>(entity, 0.0);
  } // ... local_function(...)

private:
  static std::vector<std::tuple<DomainType, DomainType, R>>
  convert(const std::vector<std::pair<std::pair<Common::FieldVector<D, d>, Common::FieldVector<D, d>>, R>>& values)
  {
    std::vector<std::tuple<DomainType, DomainType, R>> ret;
    for (const auto& element : values)
      ret.emplace_back(element.first.first, element.first.second, element.second);
    return ret;
  } // convert(...)

  const std::vector<std::tuple<DomainType, DomainType, R>> values_;
  const std::string name_;
}; // class IndicatorFunction

} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_INDICATOR_HH
