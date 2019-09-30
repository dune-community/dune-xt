// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2018)
//   René Fritze     (2015 - 2016, 2018)
//   Tim Keil        (2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_FUNCTIONS_INDICATOR_HH
#define DUNE_XT_FUNCTIONS_INDICATOR_HH

#include <utility>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/numeric_cast.hh>

#include <dune/xt/functions/interfaces/function.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

namespace Dune {
namespace XT {
namespace Functions {


template <class E, size_t r, size_t rC = 1, class R = double>
class IndicatorGridFunction : public GridFunctionInterface<E, r, rC, R>
{
  using BaseType = GridFunctionInterface<E, r, rC, R>;
  using ThisType = IndicatorGridFunction;

  class LocalIndicatorGridFunction : public ElementFunctionInterface<E, r, rC, R>
  {
    using InterfaceType = ElementFunctionInterface<E, r, rC, R>;

  public:
    using typename InterfaceType::DerivativeRangeReturnType;
    using typename InterfaceType::DomainType;
    using typename InterfaceType::ElementType;
    using typename InterfaceType::RangeReturnType;
    using typename InterfaceType::RangeType;
    using GeometryType = typename ElementType::Geometry;

    LocalIndicatorGridFunction(
        const std::vector<std::tuple<DomainType, DomainType, RangeType>>& subdomain_and_value_tuples)
      : InterfaceType()
      , subdomain_and_value_tuples_(subdomain_and_value_tuples)
    {}

  protected:
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
    const std::vector<std::tuple<DomainType, DomainType, RangeType>>& subdomain_and_value_tuples_;
    RangeType current_value_;
  }; // class LocalIndicatorGridFunction


public:
  using DomainType = typename LocalIndicatorGridFunction::DomainType;
  using RangeType = typename LocalIndicatorGridFunction::RangeType;

  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::ElementType;
  using typename BaseType::LocalFunctionType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".indicator";
  }

  static Common::Configuration defaults()
  {
    Common::Configuration config;
    config["type"] = static_id();
    config["0.domain"] = "[0 1; 0 1; 0 1]";
    config["0.value"] = "[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]";
    config["1.domain"] = "[-1 0.5; -1 0.5; -1 0.5]";
    config["1.value"] = "[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]";
    config["name"] = static_id();
    return config;
  } // ... defaults(...)

  /**
   *        Can be used for declaration of lower left corner and upper right corner of the desired domains.
\code
FunctionType function({{lowerleft_1, upperright_1, value_1}, {lowerleft_2, upperright_2, value_2}});
\endcode
   */
  IndicatorGridFunction(const std::vector<std::tuple<DomainType, DomainType, RangeType>>& values,
                        const std::string name_in = "indicator")
    : subdomain_and_value_tuples_(values)
    , name_(name_in)
  {}

  /**
   * \brief Convenience ctor.
   *
   *        Can be used as in
\code
FunctionType function({{{{0., 1.}, {0., 1.}}, 0.7}, {{{6., 10.}, {8., 10.}}, 0.9}});
\endcode
   * if you want to set indicator intervals [0,1] x [0,1] with value 0.7 and [6,10] x [8,10] with value 0.9.
   */
  IndicatorGridFunction(const std::vector<std::pair<Common::FieldMatrix<D, d, 2>, RangeType>>& values,
                        const std::string name_in = "indicator")
    : subdomain_and_value_tuples_(convert_from_domains(values))
    , name_(name_in)
  {}

  std::string name() const override final
  {
    return name_;
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<LocalIndicatorGridFunction>(subdomain_and_value_tuples_);
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
}; // class IndicatorGridFunction


template <size_t d, size_t r = 1, size_t rC = 1, class R = double>
class IndicatorFunction : public FunctionInterface<d, r, rC, R>
{
  using BaseType = FunctionInterface<d, r, rC, R>;

public:
  using typename BaseType::D;
  using typename BaseType::DerivativeRangeReturnType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeReturnType;

  IndicatorFunction(const std::vector<std::tuple<DomainType, DomainType, RangeReturnType>>& values,
                    const std::string nm = "indicator")
    : subdomain_and_value_tuples_(values)
    , name_(nm)
  {}

  IndicatorFunction(const std::vector<std::pair<Common::FieldMatrix<D, d, 2>, RangeReturnType>>& values,
                    const std::string nm = "indicator")
    : subdomain_and_value_tuples_(convert_from_domains(values))
    , name_(nm)
  {}

  int order(const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    return 0;
  }

  static std::string static_id()
  {
    return "dune.xt.functions.indicatorfunction";
  }

  static Common::Configuration defaults()
  {
    Common::Configuration config;
    config["type"] = static_id();
    config["0.domain"] = "[0 1; 0 1; 0 1]";
    config["0.value"] = "[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]";
    config["1.domain"] = "[-1 0.5; -1 0.5; -1 0.5]";
    config["1.value"] = "[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]";
    config["name"] = static_id();
    return config;
  } // ... defaults(...)

  std::string name() const override final
  {
    return "dune.xt.functions.indicatorfunction";
  }

  using BaseType::evaluate;

  RangeReturnType evaluate(const DomainType& point_in_global_coordinates,
                           const Common::Parameter& /*param*/ = {}) const override final
  {
    RangeReturnType value(0.);
    for (const auto& subdomain_and_value_tuple : subdomain_and_value_tuples_) {
      const auto& subdomain_ll = std::get<0>(subdomain_and_value_tuple);
      const auto& subdomain_ur = std::get<1>(subdomain_and_value_tuple);
      if (Common::FloatCmp::le(subdomain_ll, point_in_global_coordinates)
          && Common::FloatCmp::le(point_in_global_coordinates, subdomain_ur))
        value += std::get<2>(subdomain_and_value_tuple);
    }
    return value;
  } // ... evaluate(...)

  using BaseType::jacobian;

  DerivativeRangeReturnType jacobian(const DomainType& /*point_in_global_coordinates*/,
                                     const Common::Parameter& /*param*/ = {}) const override final
  {
    return DerivativeRangeReturnType(); // <- defaults to 0
  }

private:
  static std::vector<std::tuple<DomainType, DomainType, RangeReturnType>>
  convert_from_domains(const std::vector<std::pair<Common::FieldMatrix<D, d, 2>, RangeReturnType>>& values)
  {
    std::vector<std::tuple<DomainType, DomainType, RangeReturnType>> ret;
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

  const std::vector<std::tuple<DomainType, DomainType, RangeReturnType>> subdomain_and_value_tuples_;
  const std::string name_;
}; // class IndicatorFunction


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_INDICATOR_HH
