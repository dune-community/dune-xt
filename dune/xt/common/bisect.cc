// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2019)
//   Tobias Leibner  (2019 - 2020)

#include <functional>

#include "config.h"

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/float_cmp.hh>

namespace Dune {
namespace XT {
namespace Common {


double find_largest_by_bisection(const double& left,
                                 const double& right,
                                 std::function<bool(const double&)> condition,
                                 const double rel_error,
                                 const double abs_error)
{
  double ll = (left < right) ? left : right;
  double rr = (left < right) ? right : left;
  if (condition(rr))
    return rr;
  DUNE_THROW_IF(!condition(ll), Exceptions::bisection_error, "");
  // no we know that ll is good, rr is bad
  while (FloatCmp::gt(rr, ll, rel_error, abs_error)) {
    const double middle = 0.5 * (ll + rr);
    if (condition(middle))
      ll = middle;
    else
      rr = middle;
  }
  return ll;
} // ... find_largest_by_bisection(...)


} // namespace Common
} // namespace XT
} // namespace Dune
