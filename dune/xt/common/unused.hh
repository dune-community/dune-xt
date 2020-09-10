// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2015 - 2016, 2018 - 2019)
//   Tobias Leibner  (2020)

#ifndef DXTC_UNUSED_HH
#define DXTC_UNUSED_HH

#ifdef HAS_WORKING_UNUSED_ATTRIBUTE
#  define DXTC_UNUSED(identifier) identifier __attribute__((unused))
#else
#  define DXTC_UNUSED(identifier) /*identifier*/
#endif

#ifndef NDEBUG
#  define DXTC_DEBUG_ONLY(param) param
#else
#  define DXTC_DEBUG_ONLY(param) DXTC_UNUSED(param)
#endif

#if defined(HAVE_MPI) && HAVE_MPI
#  define DXTC_MPI_ONLY(param) param
#else
#  define DXTC_MPI_ONLY(param) DXTC_UNUSED(param)
#endif

#endif // DXTC_UNUSED_HH
