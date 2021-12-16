// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017, 2020)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2017, 2020)

#ifndef DUNE_XT_LA_PRINT_HH
#define DUNE_XT_LA_PRINT_HH

#include <dune/xt/common/print.hh>

namespace Dune::XT::LA {


// Without the NOLINTs, these declarations are removed by clang-tidy since they are unused used in some compilation
// units (in particular the headerchecks)
// NOLINTNEXTLINE(misc-unused-using-decls)
using XT::Common::print;
// NOLINTNEXTLINE(misc-unused-using-decls)
using XT::Common::repr;


} // namespace Dune::XT::LA

#endif // DUNE_XT_LA_PRINT_HH
