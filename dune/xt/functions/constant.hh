// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Kirsten Weber   (2013)
//   Rene Milk       (2013 - 2018)
//   Tobias Leibner  (2014 - 2015, 2017)

#ifndef DUNE_XT_FUNCTIONS_CONSTANT_HH
#define DUNE_XT_FUNCTIONS_CONSTANT_HH

#include <memory>

#include <dune/xt/common/configuration.hh>

#include "interfaces.hh"

namespace Dune {
namespace XT {
namespace Functions {
namespace internal {


//! TODO move
template <class K, int rows, int cols>
struct EyeMatrix
{
  typedef FieldMatrix<K, rows, cols> type;

  static type value()
  {
    type ret(0);
    for (size_t dd = 0; dd < std::min(rows, cols); ++dd)
      ret[dd][dd] = 1;
    return ret;
  }

  static std::string value_str()
  {
    std::string str = "[";
    for (size_t rr = 0; rr < rows; ++rr) {
      if (rr > 0)
        str += "; ";
      for (size_t cc = 0; cc < cols; ++cc) {
        if (cc > 0)
          str += " ";
        if (cc == rr)
          str += "1";
        else
          str += "0";
      }
    }
    str += "]";
    return str;
  }
}; // struct EyeMatrix

template <class K>
struct EyeMatrix<K, 1, 1>
{
  typedef FieldVector<K, 1> type;

  static type value()
  {
    return type(1);
  }

  static std::string value_str()
  {
    return "1";
  }
}; // struct EyeMatrix< K, 1, 1 >

template <class K, int dim>
struct UnitMatrix : EyeMatrix<K, dim, dim>
{
};

template <class K, int dim>
typename UnitMatrix<K, dim>::type unit_matrix()
{
  return UnitMatrix<K, dim>::value();
}


} // namespace internal


template <size_t d, size_t r = 1, size_t rC = 1, class RangeFieldImp = double>
class ConstantFunction : public SmoothFunctionInterface<d, r, rC, RangeFieldImp>
{
  using BaseType = SmoothFunctionInterface<d, r, rC, RangeFieldImp>;
  using ThisType = ConstantFunction<d, r, rC, RangeFieldImp>;

public:
  using DerivativeRangeType = typename BaseType::DerivativeRangeType;
  using DomainType = typename BaseType::DomainType;
  using RangeType = typename BaseType::RangeType;
  using SingleDerivativeRangeType = typename BaseType::SingleDerivativeRangeType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".constant";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["type"] = static_id();
    config["value"] = "[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]";
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
    // get correct config
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    return Common::make_unique<ThisType>(cfg.get<RangeType>("value"),
                                         cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)

  explicit ConstantFunction(const RangeType& constant, const std::string name_in = static_id())
    : constant_(constant)
    , name_(name_in)
  {
  }

  explicit ConstantFunction(const RangeFieldImp& constant, const std::string name_in = static_id())
    : constant_(constant)
    , name_(name_in)
  {
  }

#if !DUNE_XT_WITH_PYTHON_BINDINGS
  ConstantFunction(const ThisType& other) = default;
#endif

  int order(const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    return 0;
  }

  RangeType evaluate(const DomainType& /*x*/, const Common::Parameter& /*mu*/ = {}) const override final
  {
    return constant_;
  }

  DerivativeRangeType jacobian(const DomainType& /*x*/, const Common::Parameter& /*mu*/ = {}) const override final
  {
    return DerivativeRangeType(0.);
  }

  std::string name() const override final
  {
    return name_;
  }

  std::string type() const
  {
    return BaseType::static_id() + ".constant";
  }

  const RangeType constant_;
  const std::string name_;
};


} // namespace Functions
} // namespace XT
} // namespace Dune


#include "constant.lib.hh"


#endif // DUNE_XT_FUNCTIONS_CONSTANT_HH
