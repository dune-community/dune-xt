// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2015 - 2016, 2018 - 2019)

#ifndef DUNE_XT_COMMON_COMPILER_HH
#define DUNE_XT_COMMON_COMPILER_HH

#if defined(__GNUC__)
#  if defined(__GNUC_PATCHLEVEL__)
#    define DXTC_GNUC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#  else
#    define DXTC_GNUC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100)
#  endif
#endif

#endif // DUNE_XT_COMMON_COMPILER_HH
