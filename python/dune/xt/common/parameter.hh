// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_XT_COMMON_PARAMETER_HH
#define PYTHON_DUNE_XT_COMMON_PARAMETER_HH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/cast.h>

#include <dune/xt/common/parameter.hh>


PYBIND11_NAMESPACE_BEGIN(pybind11)
PYBIND11_NAMESPACE_BEGIN(detail)


template <typename Type, typename Value>
struct simple_dict_caster
{
  using Key = std::string;
  using key_conv = make_caster<Key>;
  using value_conv = make_caster<Value>;

  bool load(handle src, bool convert)
  {
    if (!isinstance<dict>(src))
      return false;
    auto d = reinterpret_borrow<dict>(src);
    value.clear();
    for (auto it : d) {
      key_conv kconv;
      value_conv vconv;
      if (!kconv.load(it.first.ptr(), convert) || !vconv.load(it.second.ptr(), convert))
        return false;
      value.set(cast_op<Key&&>(std::move(kconv)), cast_op<Value&&>(std::move(vconv)));
    }
    return true;
  } // ... load(...)

  static handle cast(const Type& src, return_value_policy policy, handle parent)
  {
    dict d;
    return_value_policy policy_key = policy;
    return_value_policy policy_value = policy;
    if (!std::is_lvalue_reference<Type>::value) {
      policy_key = return_value_policy_override<Key>::policy(policy_key);
      policy_value = return_value_policy_override<Value>::policy(policy_value);
    }
    for (auto&& key : src.keys()) {
      auto key_obj = reinterpret_steal<object>(key_conv::cast(key, policy_key, parent));
      auto value_obj = reinterpret_steal<object>(value_conv::cast(src.get(key), policy_value, parent));
      if (!key_obj || !value_obj)
        return handle();
      d[key_obj] = value_obj;
    }
    return d.release();
  }

  PYBIND11_TYPE_CASTER(Type, _("Dict[") + key_conv::name + _(", ") + value_conv::name + _("]"));
}; // struct simple_dict_caster


template <>
struct type_caster<Dune::XT::Common::ParameterType> : public simple_dict_caster<Dune::XT::Common::ParameterType, size_t>
{};

template <>
struct type_caster<Dune::XT::Common::Parameter>
  : public simple_dict_caster<Dune::XT::Common::Parameter, std::vector<double>>
{};


PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(pybind11)

#endif // PYTHON_DUNE_XT_COMMON_PARAMETER_HH
