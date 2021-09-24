// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2015 - 2016, 2018 - 2020)
//   Tobias Leibner  (2020 - 2021)

#ifndef DXTC_UNUSED_HH
#define DXTC_UNUSED_HH

#ifndef NDEBUG
#  define DXTC_DEBUG_ONLY
#else
#  if HAS_WORKING_UNUSED_ATTRIBUTE
#    define DXTC_DEBUG_ONLY [[maybe_unused]]
#  else
#    define DXTC_DEBUG_ONLY __attribute__((unused))
#  endif
#endif

#if defined(HAVE_MPI) && HAVE_MPI
#  define DXTC_MPI_ONLY
#else
#  if HAS_WORKING_UNUSED_ATTRIBUTE
#    define DXTC_MPI_ONLY [[maybe_unused]]
#  else
#    define DXTC_MPI_ONLY __attribute__((unused))
#  endif
#endif

#endif // DXTC_UNUSED_HH
