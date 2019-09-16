// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018 - 2019)
//   René Fritze     (2019)
//   Tobias Leibner  (2019)

/// \file Should be included before any other DUNE header.
///
/// Using e.g. the DynamicVector with a field type of long unsigned int failes in dune/common/densevector.hh in
/// one_norm() due to ambiguity of std::abs. We thus need the dirty hacks below and we need them to be included before
/// anyone else includes any DUNE header which pulls dune/common/densevector.hh.

#ifndef DUNE_XT_COMMON_FIX_AMBIGUOUS_STD_MATH_OVERLOADS_HH
#define DUNE_XT_COMMON_FIX_AMBIGUOUS_STD_MATH_OVERLOADS_HH

#include <stdexcept>
#include <cstdint>

#include <dune/common/bigunsignedint.hh>

#include <dune/xt/common/exceptions.hh>


namespace std {


long unsigned int abs(const long unsigned int& value);


unsigned char abs(unsigned char value);


template <int k>
Dune::bigunsignedint<k> abs(const Dune::bigunsignedint<k>& value)
{
  return value;
}


template <int k>
inline Dune::bigunsignedint<k> pow(Dune::bigunsignedint<k> /*value*/, std::uintmax_t /*n*/)
{
  DUNE_THROW(Dune::NotImplemented, "pow is not available for Dune::bigunisgnedint!");
  return Dune::bigunsignedint<k>();
}


template <int k>
inline Dune::bigunsignedint<k> sqrt(Dune::bigunsignedint<k> value)
{
  DUNE_THROW(Dune::NotImplemented, "sqrt is not available for Dune::bigunisgnedint!");
  return Dune::bigunsignedint<k>(std::sqrt(value.todouble()));
}


template <int k>
inline Dune::bigunsignedint<k> conj(Dune::bigunsignedint<k> value)
{
  return value;
}


template <int k>
inline bool isnan(Dune::bigunsignedint<k> /*value*/)
{
  return false;
}


template <int k>
inline bool isinf(Dune::bigunsignedint<k> /*value*/)
{
  return false;
}


} // namespace std


#endif // DUNE_XT_COMMON_FIX_AMBIGUOUS_STD_MATH_OVERLOADS_HH
