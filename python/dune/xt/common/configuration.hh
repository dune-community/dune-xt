// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_COMMON_CONFIGURATION_PBH
#define DUNE_XT_COMMON_CONFIGURATION_PBH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/cast.h>

#include <dune/xt/common/configuration.hh>

NAMESPACE_BEGIN(pybind11)
NAMESPACE_BEGIN(detail)


template <>
struct type_caster<Dune::XT::Common::Configuration>
{
  typedef Dune::XT::Common::Configuration type;
  typedef std::string S;
  using key_conv = make_caster<S>;
  using value_conv = make_caster<S>;

  bool load(handle src, bool convert)
  {
    if (!isinstance<dict>(src))
      return false;
    auto d = reinterpret_borrow<dict>(src);
    key_conv kconv;
    value_conv vconv;
    value = type();
    for (auto it : d) {
      if (!kconv.load(it.first.ptr(), convert) || !vconv.load(it.second.ptr(), convert))
        return false;
      value[cast_op<S>(kconv)] = cast_op<S>(vconv);
    }
    return true;
  } // ... load(...)

  static handle cast(const type& src, return_value_policy policy, handle parent)
  {
    dict d;
    for (auto const& kv : src.flatten()) {
      auto key = reinterpret_steal<object>(key_conv::cast(kv.first, policy, parent));
      auto value = reinterpret_steal<object>(value_conv::cast(kv.second, policy, parent));
      if (!key || !value)
        return handle();
      d[key] = value;
    }
    return d.release();
  } // ...handle(...)

  PYBIND11_TYPE_CASTER(type, _("Dict[") + key_conv::name + _(", ") + value_conv::name + _("]"));
}; // struct type_caster<Dune::XT::Common::Configuration>


NAMESPACE_END(detail)
NAMESPACE_END(pybind11)

#endif // DUNE_XT_COMMON_CONFIGURATION_PBH
