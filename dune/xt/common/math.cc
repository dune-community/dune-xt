// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2014, 2016 - 2018)
//   René Fritze     (2013, 2015 - 2016, 2018 - 2019)
//   Tobias Leibner  (2018 - 2020)

#include "config.h"

#include "math.hh"

namespace Dune {
namespace XT {
namespace Common {
namespace internal {


char abs(const char& val)
{
  return val < 0 ? static_cast<char>(-val) : val;
}


} // namespace internal


//! calculates binomial coefficient for arbitrary n
double binomial_coefficient(const double n, const size_t k)
{
  double ret(1);
  for (size_t ii = 1; ii <= k; ++ii)
    ret *= (n + 1 - ii) / ii;
  return ret;
}

// this is just for the test to compile
const std::string Epsilon<std::string, false>::value = "_";


} // namespace Common
} // namespace XT
} // namespace Dune
