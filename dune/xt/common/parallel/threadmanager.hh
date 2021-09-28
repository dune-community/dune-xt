// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   René Fritze     (2013 - 2016, 2018 - 2020)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_COMMON_THREADMANAGER_HH
#define DUNE_XT_COMMON_THREADMANAGER_HH

#include <thread>

namespace Dune::XT::Common {

struct ThreadManager;
//! global singleton ThreadManager
ThreadManager& threadManager();

/** abstractions of threading functionality
 *  currently controls tbb, falls back to single-thread dummy imp
 **/
struct ThreadManager
{
  static size_t default_max_threads();

  //! return maximal number of threads possbile in the current run
  size_t max_threads();

  //! return number of current threads
  size_t current_threads();

  //! return thread number
  size_t thread();

  //! set maximal number of threads available during run
  void set_max_threads(const size_t count);

  ~ThreadManager() = default;

private:
  friend ThreadManager& threadManager();
  //! init tbb with given thread count, prepare Eigen for smp if possible
  ThreadManager();

  size_t max_threads_;
};

inline ThreadManager& threadManager()
{
  static ThreadManager tm;
  return tm;
}

} // namespace Dune::XT::Common

#endif // DUNE_XT_COMMON_THREADMANAGER_HH
