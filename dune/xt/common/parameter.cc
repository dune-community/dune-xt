// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2020)
//   René Fritze     (2018 - 2020)
//   Tobias Leibner  (2017 - 2018, 2020)

#include "config.h"

#include <algorithm>
#include <sstream>
#include <type_traits>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/numeric_cast.hh>

#include "parameter.hh"

namespace Dune::XT::Common {
namespace internal {


// ========================
// ===== SimpleDict =======
// ========================


template class SimpleDict<size_t>;
template class SimpleDict<std::vector<double>>;


} // namespace internal


// =========================
// ===== ParameterType =====
// =========================
ParameterType::ParameterType()
  : BaseType()
{}

ParameterType::ParameterType(const std::string& key)
  : BaseType(key, 1)
{}

ParameterType::ParameterType(const std::string& key, const size_t& sz)
  : BaseType(key, sz)
{}

ParameterType::ParameterType(const std::pair<std::string, size_t>& key_size_pair)
  : BaseType(key_size_pair.first, key_size_pair.second)
{}

ParameterType::ParameterType(const std::pair<const char*, int>& key_size_pair)
  : BaseType(std::string(key_size_pair.first), numeric_cast<size_t>(key_size_pair.second))
{}

ParameterType::ParameterType(const std::vector<std::pair<std::string, size_t>>& key_size_pairs)
  : BaseType(key_size_pairs)
{}

ParameterType::ParameterType(BaseType&& source)
  : BaseType(std::move(source))
{}

ParameterType ParameterType::operator+(const ParameterType& other) const
{
  return this->merge(
      other,
      [](const auto& left, const auto& right) { return left == right; },
      [](const auto& /*left*/, const auto& /*right*/) {
        return "cannot add parameter types which contain the same key with different sizes:";
      });
} // ... operator+(...)

bool ParameterType::operator==(const ParameterType& other) const
{
  if (this->size() == 1 && other.size() == 1) {
    const auto& this_single_element = *this->dict_.begin();
    const auto& other_single_element = *other.dict_.begin();
    if (this_single_element.first == other_single_element.first || this_single_element.first == "__unspecified__"
        || other_single_element.first == "__unspecified__") {
      return this_single_element.second == other_single_element.second;
    }
    return false;
  }
  return this->dict_ == other.dict_;
} // ... operator==(...)

bool ParameterType::operator!=(const ParameterType& other) const
{
  return !(*this == other);
}

bool ParameterType::operator<(const ParameterType& other) const
{
  if (this->dict_.size() >= other.dict_.size())
    return false;
  for (const auto& this_element : this->dict_) {
    const auto& this_key = this_element.first;
    const auto& this_keys_lenght = this_element.second;
    const auto other_key_search_result = other.dict_.find(this_key);
    if (other_key_search_result == other.dict_.end())
      return false;
    const auto& other_key_lenght = other_key_search_result->second;
    if (other_key_lenght != this_keys_lenght)
      return false;
  }
  // now we know that
  // - each of our keys is contained in other
  // - our length for each key equals the length of others keys
  // - other contains more keys
  return true;
} // ... operator<(...)

bool ParameterType::operator<=(const ParameterType& other) const
{
  return *this < other || *this == other;
}

std::string ParameterType::report() const
{
  return "ParameterType(" + BaseType::report("ParameterType(") + ")";
}


std::ostream& operator<<(std::ostream& out, const ParameterType& param_type)
{
  out << param_type.report();
  return out;
}


// =====================
// ===== Parameter =====
// =====================
Parameter::Parameter(const double& value)
  : BaseType("__unspecified__", {value})
{}

Parameter::Parameter(const std::vector<double>& value)
  : BaseType("__unspecified__", value)
{}

Parameter::Parameter(const std::string& key, const double& value)
  : BaseType(key, {value})
{}

Parameter::Parameter(const std::string& key, const ValueType& value)
  : BaseType(key, value)
{}

Parameter::Parameter(const std::vector<std::pair<std::string, ValueType>>& key_value_pairs)
  : BaseType(key_value_pairs)
{}

Parameter::Parameter(const std::initializer_list<std::pair<std::string, ValueType>>& key_value_pairs)
  : BaseType(key_value_pairs)
{}

Parameter::Parameter(BaseType&& source)
  : BaseType(std::move(source))
{}

Parameter Parameter::operator+(const Parameter& other) const
{
  return this->merge(
      other,
      [](const auto& left, const auto& right) {
        if (left.size() != right.size())
          return false;
        return FloatCmp::eq(left, right);
      },
      [](const auto& left, const auto& right) {
        if (left.size() != right.size())
          return "cannot add parameters which contain the same key with different sizes:";
        return "cannot add parameters which contain the same key with different values:";
      });
} // ... operator+(...)

bool Parameter::operator<(const Parameter& other) const
{
  return this->dict_ < other.dict_;
}

ParameterType Parameter::type() const
{
  ParameterType ret;
  for (const auto& key_value_pairs : dict_)
    ret.set(key_value_pairs.first, key_value_pairs.second.size());
  return ret;
}

std::string Parameter::report() const
{
  return "Parameter(" + BaseType::report("Parameter(") + ")";
}


std::ostream& operator<<(std::ostream& out, const Parameter& mu)
{
  out << mu.report();
  return out;
}


// ===============================
// ===== ParametricInterface =====
// ===============================
ParametricInterface::ParametricInterface(const ParameterType& param_type)
  : parameter_type_(param_type)
{}

bool ParametricInterface::is_parametric() const
{
  return !parameter_type().empty();
}

const ParameterType& ParametricInterface::parameter_type() const
{
  return parameter_type_;
}

Parameter ParametricInterface::parse_parameter(const Parameter& mu) const
{
  const auto& this_type = this->parameter_type();
  const auto& mus_type = mu.type();
  if (this_type.size() == 1 && mus_type.size() == 1) {
    // both have only one key, so either key might be '__unspecified__'
    const auto this_single_key = this_type.keys().at(0);
    const auto mus_single_key = mus_type.keys().at(0);
    if (this_single_key == mus_single_key)
      return mu;
    if (this_single_key == "__unspecified__") {
      return mu;
    }
    if (mus_single_key == "__unspecified__") {
      return Parameter(this_single_key, mu.get("__unspecified__"));
    }
    // both have only one key, but the keys don't match and neither is '__unspecified__'
    DUNE_THROW(Exceptions::parameter_error,
               "this->parameter_type() = " << this_type << "\n   "
                                           << "mu.type() = " << mus_type);
  }
  // one of them has more than one key, so '__unspecified__' does not play a role here
  if (this_type <= mus_type)
    return mu;
  DUNE_THROW(Exceptions::parameter_error,
             "this->parameter_type() = " << this_type << "\n   "
                                         << "mu.type() = " << mus_type);
  return Parameter();
} // ... parse_parameter(...)


void ParametricInterface::extend_parameter_type(const ParameterType& additional_parameter_type)
{
  try {
    this->parameter_type_ = this->parameter_type() + additional_parameter_type;
  } catch (Exceptions::parameter_error& ee) {
    DUNE_THROW(Exceptions::parameter_error, "Error extending parameter type:\n   " << ee.what());
  }
} // ... extend_parameter_type(...)


} // namespace Dune::XT::Common
