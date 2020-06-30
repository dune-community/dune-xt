// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2018)
//   René Fritze     (2015 - 2016, 2018 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2014, 2018 - 2020)

#ifndef DUNE_XT_COMMON_FMATRIX_HH
#define DUNE_XT_COMMON_FMATRIX_HH

#include <dune/common/version.hh>

#if DUNE_VERSION_GTE(DUNE_COMMON, 2, 7)
#  include "fmatrix-2.7.hh"
#else
#  include "fmatrix-2.6.hh"
#endif

#endif // DUNE_XT_COMMON_FMATRIX_HH
