// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2019)
//   René Fritze     (2017 - 2019)
//   Tobias Leibner  (2017 - 2020)

#ifndef DUNE_XT_COMMON_PARAMETER_HH
#define DUNE_XT_COMMON_PARAMETER_HH

#include <algorithm>
#include <string>
#include <map>
#include <vector>
#include <iosfwd>

#include "float_cmp.hh"
#include "string.hh"

namespace Dune::XT::Common {
namespace internal {


template <class ValueType>
class SimpleDict
{
public:
  using ThisType = SimpleDict;

  SimpleDict() = default;

  SimpleDict(const ThisType& other) = default;

  SimpleDict(ThisType&& source) = default;

  SimpleDict(const std::string& key, const ValueType& value)
    : dict_({std::make_pair(key, value)})
  {
    update_keys();
  }

  SimpleDict(const std::vector<std::pair<std::string, ValueType>>& key_value_pairs)
  {
    for (const auto& key_value_pair : key_value_pairs)
      dict_.emplace(key_value_pair);
    update_keys();
  }

  ThisType& operator=(const ThisType& other) = default;

  ThisType& operator=(ThisType&& source) = default;

  const std::vector<std::string>& keys() const
  {
    return keys_;
  }

  void clear()
  {
    keys_.clear();
    dict_.clear();
  }

  bool empty() const
  {
    return dict_.size() == 0;
  }

  bool has_key(const std::string& key) const
  {
    return dict_.find(key) != dict_.end();
  }

  void set(const std::string& key, const ValueType& value, const bool overwrite = false)
  {
    if (key.empty())
      DUNE_THROW(Exceptions::parameter_error, "Given key must not be empty!");
    const bool key_was_present = has_key(key);
    if (!overwrite && key_was_present)
      DUNE_THROW(Exceptions::parameter_error,
                 "You are trying to overwrite the key '"
                     << key << "' (although a value is already set), and overwrite is false!");
    dict_[key] = value;
    if (!key_was_present)
      update_keys();
  } // ... set(...)

  const ValueType& get(const std::string& key) const
  {
    const auto result = dict_.find(key);
    if (result == dict_.end())
      DUNE_THROW(Exceptions::parameter_error, "Key '" << key << "' does not exist!");
    return result->second;
  }

  size_t size() const
  {
    return dict_.size();
  }

protected:
  std::string report(const std::string& prefix) const
  {
    if (dict_.empty())
      return "{}";
    DXT_ASSERT(keys_.size() > 0);
    const auto whitespaced_prefix = whitespaceify(prefix);
    std::stringstream ss;
    ss << "{" << keys_[0] << ": " << dict_.at(keys_[0]);
    for (size_t ii = 1; ii < keys_.size(); ++ii) {
      ss << ",\n" << whitespaced_prefix << " " << keys_[ii] << ": " << dict_.at(keys_[ii]);
    }
    ss << "}";
    return ss.str();
  } // ... report(...)

  SimpleDict merge(const SimpleDict& other,
                   std::function<bool(ValueType, ValueType)> value_comparator,
                   std::function<std::string(ValueType, ValueType)> error_msg_prefix) const
  {
    if (this->empty())
      return other;
    if (other.empty())
      return *this;
    SimpleDict<ValueType> ret = *this;
    for (const auto& other_element : other.dict_) {
      const auto& other_key = other_element.first;
      const auto& other_value = other_element.second;
      const auto this_key_search_result = dict_.find(other_key);
      if (this_key_search_result == dict_.end()) {
        // key of others is not contained in this, just add to ret
        ret.set(other_key, other_value);
      } else {
        // key of other is also present in this
        const auto& this_value = this_key_search_result->second;
        DUNE_THROW_IF(!value_comparator(this_value, other_value),
                      Exceptions::parameter_error,
                      error_msg_prefix(this_value, other_value)
                          << "\n   this->get(\"" << other_key << "\") = " << this_value << "\n   other.get(\""
                          << other_value << "\")");
        // and the respective values agree, so no need to do something
      }
    }
    return ret;
  } // ... merge(...)

  void update_keys()
  {
    keys_.resize(dict_.size());
    size_t ii = 0;
    for (const auto& key_value_pair : dict_) {
      keys_[ii] = key_value_pair.first;
      ++ii;
    }
    std::sort(keys_.begin(), keys_.end());
  } // ... update_keys(...)

  std::map<std::string, ValueType> dict_;
  std::vector<std::string> keys_;
}; // class SimpleDict


extern template class SimpleDict<size_t>;
extern template class SimpleDict<std::vector<double>>;


} // namespace internal


class ParameterType : public internal::SimpleDict<size_t>
{
  using BaseType = internal::SimpleDict<size_t>;

public:
  ParameterType();

  ParameterType(const ParameterType& other) = default;

  ParameterType(ParameterType&& source) = default;

private:
  ParameterType(BaseType&& source);

public:
  ParameterType(const std::string& key);

  ParameterType(const std::string& key, const size_t& sz);

  ParameterType(const std::pair<std::string, size_t>& key_size_pair);

  ParameterType(const std::pair<const char*, int>& key_size_pair);

  ParameterType(const std::vector<std::pair<std::string, size_t>>& key_size_pairs);

  ParameterType& operator=(const ParameterType& other) = default;

  ParameterType& operator=(ParameterType&& source) = default;

  ParameterType operator+(const ParameterType& other) const;

  /**
   * \note In the special case that this and other both have only a single key, and either of the keys is
   *       '__unspecified__', then they compare equal if the sizes corresponding to these keys compare equal.
   * \sa   Take a look at test/parameter.cc for examples.
   */
  bool operator==(const ParameterType& other) const;

  bool operator!=(const ParameterType& other) const;

  /**
   * Returns true, if each of our keys is contained in other (with equal respective lengths), and other contains
   * additional keys.
   */
  bool operator<(const ParameterType& other) const;

  bool operator<=(const ParameterType& other) const;

  std::string report() const;
}; // class ParameterType


std::ostream& operator<<(std::ostream& out, const ParameterType& param_type);


class Parameter : public internal::SimpleDict<std::vector<double>>
{
  using BaseType = internal::SimpleDict<std::vector<double>>;
  using ValueType = std::vector<double>;

public:
  Parameter(const Parameter& other) = default;

  Parameter(Parameter&& source) = default;

private:
  Parameter(BaseType&& source);

public:
  /// \note this is somehow necessary to make clang 3.8 happy (and cannot be defaulted)
  ~Parameter() = default;

  Parameter(const std::vector<std::pair<std::string, ValueType>>& key_value_pairs = {});

  Parameter(const std::initializer_list<std::pair<std::string, ValueType>>& key_value_pairs);

  /**
   * \brief Same as Parameter({"__unspecified__", value});
   */
  Parameter(const double& value);

  Parameter(const std::vector<double>& value);

  Parameter(const std::string& key, const double& value);

  Parameter(const std::string& key, const ValueType& value);

  Parameter& operator=(const Parameter& other) = default;

  Parameter& operator=(Parameter&& source) = default;

  Parameter operator+(const Parameter& other) const;

  bool operator<(const Parameter& other) const;

  ParameterType type() const;

  std::string report() const;
}; // class Parameter


std::ostream& operator<<(std::ostream& out, const Parameter& mu);


class ParametricInterface
{
public:
  ParametricInterface(const ParameterType& param_type = {});

  ParametricInterface(const ParametricInterface& other) = default;

  ParametricInterface(ParametricInterface&& source) = default;

  virtual ~ParametricInterface() = default;

  ParametricInterface& operator=(const ParametricInterface& other) = default;

  ParametricInterface& operator=(ParametricInterface&& source) = default;

  virtual bool is_parametric() const;

  virtual const ParameterType& parameter_type() const;

  Parameter parse_parameter(const Parameter& mu) const;

protected:
  void extend_parameter_type(const ParameterType& additional_parameter_type);

private:
  ParameterType parameter_type_;
}; // class ParametricInterface


} // namespace Dune::XT::Common

#endif // DUNE_XT_COMMON_PARAMETER_HH
