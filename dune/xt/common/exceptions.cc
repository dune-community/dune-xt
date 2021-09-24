// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017, 2020)
//   René Fritze     (2017 - 2020)
//   Tobias Leibner  (2017 - 2018, 2020)

#include <config.h>

#include <iostream>

#if HAVE_TBB && __has_include(<tbb/tbb_exception.h>)
#  include <tbb/tbb_exception.h>
#endif

#include <dune/xt/common/timings.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/parallel/helper.hh>

#include "exceptions.hh"

namespace Dune::XT::Common {


int handle_exception(const Dune::Exception& exp)
{
  std::cerr << "Failed with Dune::Exception: " << exp.what();
  DXTC_TIMINGS.output_per_rank("profiler");
  mem_usage();
  return abort_all_mpi_processes();
}


int handle_exception(const std::exception& exp)
{
  std::cerr << "Failed with std::exception: " << exp.what();
  DXTC_TIMINGS.output_per_rank("profiler");
  mem_usage();
  return abort_all_mpi_processes();
}


} // namespace Dune::XT::Common
