// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Kirsten Weber   (2013)
//   Rene Milk       (2013 - 2016)
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


template <class R, size_t r, size_t rC>
struct Get
{
  static std::string value_str()
  {
    std::string str = "[";
    for (size_t rr = 0; rr < r; ++rr) {
      if (rr > 0)
        str += "; ";
      for (size_t cc = 0; cc < rC; ++cc) {
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
};

template <class R, size_t rC>
struct Get<R, 1, rC>
{
  static std::string value_str()
  {
    std::string str = "[";
    for (size_t cc = 0; cc < rC; ++cc) {
      if (cc > 0)
        str += " ";
      str += "1";
    }
    str += "]";
    return str;
  }
};

template <class R, size_t r>
struct Get<R, r, 1>
{
  static std::string value_str()
  {
    return Get<R, 1, r>::value_str();
  }
};

template <class R>
struct Get<R, 1, 1>
{
  static std::string value_str()
  {
    return "1";
  }
};

} // namespace internal

template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class ConstantFunction
    : public GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> BaseType;
  typedef ConstantFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;

public:
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  using typename BaseType::LocalfunctionType;

  static_assert(std::is_same<typename LocalfunctionType::RangeType, RangeType>::value, "RangeType mismatch");
  static_assert(std::is_same<typename LocalfunctionType::DomainType, DomainType>::value, "DomainType mismatch");
  static_assert(std::is_same<Dune::FieldVector<DomainFieldImp, domainDim>, DomainType>::value, "RangeType mismatch");

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".constant";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["value"] = internal::Get<RangeFieldImp, rangeDim, rangeDimCols>::value_str();
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
    return Common::make_unique<ThisType>(cfg.get("value", default_cfg.get<RangeType>("value")),
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

  virtual std::string type() const override final
  {
    return BaseType::static_id() + ".constant";
  }

  virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    return 0;
  }

  virtual void
  evaluate(const DomainType& /*x*/, RangeType& ret, const Common::Parameter& /*mu*/ = {}) const override final
  {
    ret = constant_;
  }

  virtual void
  jacobian(const DomainType& /*x*/, JacobianRangeType& ret, const Common::Parameter& /*mu*/ = {}) const override final
  {
    //! TODO: why not ```ret *=0```?
    clear_jacobian<rangeDim, rangeDimCols>()(ret);
  }

  virtual std::string name() const override final
  {
    return name_;
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

  const RangeType constant_;
  const std::string name_;
};


} // namespace Functions
} // namespace XT
} // namespace Dune


#include "constant.lib.hh"


#endif // DUNE_XT_FUNCTIONS_CONSTANT_HH
