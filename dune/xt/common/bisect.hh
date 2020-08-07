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

#ifndef DUNE_XT_COMMON_BISECT_HH
#define DUNE_XT_COMMON_BISECT_HH

#include <functional>

#include <dune/xt/common/float_cmp.hh>

namespace Dune::XT::Common {


/**
 * \brief Finds the largest number x between left and right, where condition(x) is true, but condition(y) is false for
 *        any y > x.
 *
 * \note Presumes that: if condition(x_1) == condition(x_2) == value,
 *                      then there does not exist a x_1 < y < x_2, s.t. condition(y) == !value;
 */
double find_largest_by_bisection(const double& left,
                                 const double& right,
                                 std::function<bool(const double&)> condition,
                                 const double rel_error = FloatCmp::DefaultEpsilon<double>::value(),
                                 const double abs_error = FloatCmp::DefaultEpsilon<double>::value());


} // namespace Dune::XT::Common

#endif // DUNE_XT_COMMON_BISECT_HH
