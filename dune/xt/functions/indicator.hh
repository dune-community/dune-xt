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

#include <utility>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/numeric_cast.hh>

#include <dune/xt/functions/interfaces/grid-function.hh>

namespace Dune {
namespace XT {
namespace Functions {


template <class E, size_t r, size_t rC = 1, class R = double>
class IndicatorFunction : public GridFunctionInterface<E, r, rC, R>
{
  using BaseType = GridFunctionInterface<E, r, rC, R>;
  using ThisType = IndicatorFunction<E, r, rC, R>;

  class LocalIndicatorFunction : public ElementFunctionInterface<E, r, rC, R>
  {
    using InterfaceType = ElementFunctionInterface<E, r, rC, R>;

  public:
    using typename InterfaceType::ElementType;
    using typename InterfaceType::DomainType;
    using typename InterfaceType::RangeReturnType;
    using typename InterfaceType::RangeType;
    using typename InterfaceType::DerivativeRangeReturnType;
    using GeometryType = typename ElementType::Geometry;

    LocalIndicatorFunction(const std::vector<std::tuple<DomainType, DomainType, RangeType>>& subdomain_and_value_tuples)
      : InterfaceType()
      , subdomain_and_value_tuples_(subdomain_and_value_tuples)
    {
    }

    void post_bind(const ElementType& element) override final
    {
      current_value_ = 0.;
      const auto center = element.geometry().center();
      for (const auto& subdomain_and_value_tuple : subdomain_and_value_tuples_) {
        const auto& subdomain_ll = std::get<0>(subdomain_and_value_tuple);
        const auto& subdomain_ur = std::get<1>(subdomain_and_value_tuple);
        if (Common::FloatCmp::le(subdomain_ll, center) && Common::FloatCmp::lt(center, subdomain_ur))
          current_value_ += std::get<2>(subdomain_and_value_tuple);
      }
    } // ... post_bind(...)

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
    const std::vector<std::tuple<DomainType, DomainType, RangeType>>& subdomain_and_value_tuples_;
    RangeType current_value_;
  }; // class LocalIndicatorFunction


public:
  using DomainType = typename LocalIndicatorFunction::DomainType;
  using RangeType = typename LocalIndicatorFunction::RangeType;

  using typename BaseType::ElementType;
  using typename BaseType::LocalFunctionType;
  using BaseType::d;
  using typename BaseType::D;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".indicator";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["type"] = static_id();
    config["0.domain"] = "[0 1; 0 1; 0 1]";
    config["0.value"] = "[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]";
    config["1.domain"] = "[-1 0.5; -1 0.5; -1 0.5]";
    config["1.value"] = "[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]";
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
                                          const std::string sub_name = "")
  {
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration def_cfg = default_config();
    std::vector<std::tuple<DomainType, DomainType, RangeType>> values;
    DomainType tmp_lower;
    DomainType tmp_upper;
    size_t cc = 0;
    while (cfg.has_sub(Common::to_string(cc))) {
      const Common::Configuration local_cfg = cfg.sub(Common::to_string(cc));
      if (local_cfg.has_key("domain") && local_cfg.has_key("value")) {
        auto domains = local_cfg.get<FieldMatrix<D, d, 2>>("domain");
        for (size_t dd = 0; dd < d; ++dd) {
          tmp_lower[dd] = domains[dd][0];
          tmp_upper[dd] = domains[dd][1];
        }
        auto val = local_cfg.get<RangeType>("value");
        values.emplace_back(tmp_lower, tmp_upper, val);
      } else
        break;
      ++cc;
    }
    return Common::make_unique<ThisType>(values, cfg.get("name", def_cfg.get<std::string>("name")));
  } // ... create(...)

  /**
   *        Can be used for declaration of lower left corner and upper right corner of the desired domains.
\code
FunctionType function({{lowerleft_1, upperright_1, value_1}, {lowerleft_2, upperright_2, value_2}});
\endcode
   */

  IndicatorFunction(const std::vector<std::tuple<DomainType, DomainType, RangeType>>& values,
                    const std::string name_in = "indicator")
    : subdomain_and_value_tuples_(values)
    , name_(name_in)
  {
  }

  /**
   * \brief Convenience ctor.
   *
   *        Can be used as in
\code
FunctionType function({{{{0., 1.}, {0., 1.}}, 0.7}, {{{6., 10.}, {8., 10.}}, 0.9}});
\endcode
   * if you want to set indicator intervals [0,1] x [0,1] with value 0.7 and [6,10] x [8,10] with value 0.9.
   */
  IndicatorFunction(const std::vector<std::pair<Common::FieldMatrix<D, d, 2>, RangeType>>& values,
                    const std::string name_in = "indicator")
    : subdomain_and_value_tuples_(convert_from_domains(values))
    , name_(name_in)
  {
  }

  std::string name() const override final
  {
    return name_;
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<LocalIndicatorFunction>(subdomain_and_value_tuples_);
  }

private:
  static std::vector<std::tuple<DomainType, DomainType, RangeType>>
  convert_from_domains(const std::vector<std::pair<Common::FieldMatrix<D, d, 2>, RangeType>>& values)
  {
    std::vector<std::tuple<DomainType, DomainType, RangeType>> ret;
    for (const auto& element : values) {
      DomainType tmp_coordinates_ll(0.);
      DomainType tmp_coordinates_ur(0.);
      for (size_t dd = 0; dd < d; ++dd) {
        tmp_coordinates_ll[dd] = element.first[dd][0];
        tmp_coordinates_ur[dd] = element.first[dd][1];
      }
      ret.emplace_back(tmp_coordinates_ll, tmp_coordinates_ur, element.second);
    }
    return ret;
  } // convert_from_tuples(...)

  const std::vector<std::tuple<DomainType, DomainType, RangeType>> subdomain_and_value_tuples_;
  const std::string name_;
}; // class IndicatorFunction


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_INDICATOR_HH
