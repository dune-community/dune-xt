// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014, 2016 - 2017)
//   René Fritze     (2012 - 2013, 2015 - 2016, 2018 - 2019)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_COMMON_VALIDATION_HH
#define DUNE_XT_COMMON_VALIDATION_HH

#include <ostream>
#include <string>
#include <boost/format.hpp>

#include <dune/xt/common/misc.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/float_cmp.hh>

namespace Dune {
namespace XT {
namespace Common {

/** \brief Base class for all Validators
 *  the idea is from dune-fem, only our class is an actual interface
 **/
template <class T, class Derived>
class ValidatorInterface
{
public:
  bool operator()(const T& rhs) const
  {
    return as_imp()(rhs);
  }

  std::string msg(const T& rhs) const
  {
    return as_imp().msg(rhs);
  }

protected:
  inline const Derived& as_imp() const
  {
    return static_cast<const Derived&>(*this);
  }

  inline Derived& as_imp()
  {
    return static_cast<Derived&>(*this);
  }
};

//! a class usable as a default validator
template <class T>
class ValidateAny : public ValidatorInterface<T, ValidateAny<T>>
{
  typedef ValidateAny<T> ThisType;
  typedef ValidatorInterface<T, ThisType> BaseType;

public:
  inline ValidateAny() {}
  inline ValidateAny(const ThisType&) {}

  inline bool operator()(const T&) const
  {
    return true;
  }

  std::string msg(const T&) const
  {
    return "ValidateAny: all values should be valid... \n\n";
  }
};

//! validates arg iff in given list
template <class T, class ListImp = std::vector<T>>
class ValidateInList : public ValidatorInterface<T, ValidateInList<T, ListImp>>
{
  typedef ValidateInList<T, ListImp> ThisType;
  typedef ValidatorInterface<T, ThisType> BaseType;
  typedef ListImp ListType;
  ListType valid_list_;

public:
  explicit ValidateInList(const ListType& valid_list)
    : valid_list_(valid_list)
  {}

  inline bool operator()(const T& rhs) const
  {
    return std::find(valid_list_.begin(), valid_list_.end(), rhs) != valid_list_.end();
  }

  std::string msg(const T& value) const
  {
    return (boost::format("ValidateInList: checked Parameter %s was not in valid list %s\n") % to_string(value)
            % to_string(valid_list_))
        .str();
  }
};

//! validate arg iff lhs < rhs
template <class T>
class ValidateLess : public ValidatorInterface<T, ValidateLess<T>>
{
public:
  ValidateLess(const T& lhs)
    : lhs_(lhs)
  {}
  inline bool operator()(const T& rhs) const
  {
    return FloatCmp::lt(lhs_, rhs, 0., 0.);
  }

  std::string msg(const T& value) const
  {
    return (boost::format("given value %s was invalid: not less than %s") % to_string(value) % to_string(lhs_)).str();
  }

private:
  const T lhs_;
};

//! validate arg iff lhs > rhs, obviously
template <class T>
class ValidateGreater : public ValidatorInterface<T, ValidateGreater<T>>
{
public:
  ValidateGreater(const T& lhs)
    : lhs_(lhs)
  {}
  inline bool operator()(const T& rhs) const
  {
    return FloatCmp::gt(lhs_, rhs, 0., 0.);
  }

  std::string msg(const T& value) const
  {
    return (boost::format("given value %s was invalid: not greater than %s") % to_string(value) % to_string(lhs_))
        .str();
  }

private:
  const T lhs_;
};

//! validate iff not Validator(arg)
template <class T, class Validator>
class ValidateInverse : public ValidatorInterface<T, ValidateInverse<T, Validator>>
{
public:
  ValidateInverse(const Validator validator = Validator())
    : validator_(validator)
  {}
  ValidateInverse(const T arg)
    : validator_(Validator(arg))
  {}
  inline bool operator()(const T& val) const
  {
    return !validator_(val);
  }

  std::string msg(const T& value) const
  {
    return "Inverse condition failed: " + validator_.msg(value);
  }

private:
  const Validator validator_;
};

#define INVERSE_VALIDATE(V_NEW_NAME, V_BASE_NAME)                                                                      \
  template <class T>                                                                                                   \
  struct V_NEW_NAME : public ValidateInverse<T, V_BASE_NAME<T>>                                                        \
  {                                                                                                                    \
    template <typename... Types>                                                                                       \
    V_NEW_NAME(Types... args)                                                                                          \
      : ValidateInverse<T, V_BASE_NAME<T>>(args...)                                                                    \
    {}                                                                                                                 \
  }

INVERSE_VALIDATE(ValidateNone, ValidateAny);
INVERSE_VALIDATE(ValidateGreaterOrEqual, ValidateLess);
INVERSE_VALIDATE(ValidateLessOrEqual, ValidateGreater);

//! validate arg iff rhs \in [min,max]
template <class T>
class ValidateInterval : public ValidatorInterface<T, ValidateInterval<T>>
{
public:
  ValidateInterval(const T& min, const T& max)
    : min_(min)
    , max_(max)
  {}

  inline bool operator()(const T& rhs) const
  {
    const bool lowerOk = ValidateGreaterOrEqual<T>(rhs)(min_);
    const bool upperOk = ValidateGreaterOrEqual<T>(max_)(rhs);
    return lowerOk && upperOk;
  }

  std::string msg(const T& value) const
  {
    return (boost::format("given value %s was invalid: value not in set interval [%s %s]") % to_string(value) % min_
            % max_)
        .str();
  }

private:
  const T min_;
  const T max_;
};

//! example partial specialisation
template <typename T>
struct Typename<ValidateAny<T>>
{
  static std::string value()
  {
    return "Parameter::ValidateAny<T>";
  }
};

} // end namespace Common
} // end namespace XT
} // end namespace Dune

#endif // DUNE_XT_COMMON_VALIDATION_HH
