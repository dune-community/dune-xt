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

#include <dune/xt/common/disable_warnings.hh>
#include <boost/numeric/conversion/cast.hpp>
#include <dune/xt/common/reenable_warnings.hh>

#include "math.hh"

namespace Dune {
namespace XT {
namespace Common {
namespace internal {


char abs(const char& val)
{
  return val < 0 ? static_cast<char>(-val) : val;
}

char absolute_difference(char a, char b)
{
  // calculating with chars returns an int, so we have to cast back to char
  return boost::numeric_cast<char>((a > b) ? a - b : b - a);
}


} // namespace internal


// this is just for the test to compile
const std::string Epsilon<std::string, false>::value = "_";


} // namespace Common
} // namespace XT
} // namespace Dune
