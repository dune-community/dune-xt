// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"
#include "threadmanager.hh"

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/common/configuration.hh>
#include <dune/common/exceptions.hh>

#include <dune/stuff/fem.hh>
#if HAVE_DUNE_FEM
#include <dune/fem/misc/threads/threadmanager.hh>
#endif

#if HAVE_EIGEN
#include <Eigen/Core>
#endif

// some assertions only make sense if dune-fem's threading manager is non-trivial
#if defined(USE_PTHREADS) || defined(_OPENMP)
#define WITH_DUNE_FEM_AND_THREADING(expr) WITH_DUNE_FEM(expr)
#else
#define WITH_DUNE_FEM_AND_THREADING(expr)
#endif

#if HAVE_TBB

#include <thread>
#include <tbb/task_scheduler_init.h>

size_t Dune::Stuff::ThreadManager::max_threads()
{
  const auto threads = DSC_CONFIG_GET("threading.max_count", 1);
  WITH_DUNE_FEM_AND_THREADING(assert(Dune::Fem::ThreadManager::maxThreads() == threads);)
  return threads;
}

size_t Dune::Stuff::ThreadManager::current_threads()
{
  const auto threads = max_threads();
  WITH_DUNE_FEM_AND_THREADING(assert(long(Dune::Fem::ThreadManager::currentThreads()) == long(threads));)
  return threads;
}

size_t Dune::Stuff::ThreadManager::thread()
{
  const auto tbb_id = std::this_thread::get_id();
  static std::map<decltype(tbb_id), size_t> thread_ids;
  const auto it = thread_ids.find(tbb_id);
  if (it == thread_ids.end())
    DSC::map_emplace(thread_ids, tbb_id, thread_ids.size());
  return thread_ids.at(tbb_id);
}

//! both std::hw_concur and intel's default_thread_count fail for mic
size_t Dune::Stuff::ThreadManager::default_max_threads()
{
#ifndef __MIC__
  return std::thread::hardware_concurrency();
#else
  return DS_MAX_MIC_THREADS;
#endif
}

void Dune::Stuff::ThreadManager::set_max_threads(const size_t count)
{
  DSC_CONFIG.set("threading.max_count", count, true);
  if (tbb_init_.is_active()) {
    DSC_LOG_DEBUG << (boost::format("Re-initializing TBB from %d to %d threads") % max_threads_ % count).str();
    tbb_init_.terminate();
  }
  max_threads_        = count;
  const int int_count = boost::numeric_cast<int>(count);
  WITH_DUNE_FEM(Dune::Fem::ThreadManager::setMaxNumberThreads(int_count);)
#if HAVE_EIGEN
  Eigen::setNbThreads(int_count);
#endif
  tbb_init_.initialize(int_count);
}

Dune::Stuff::ThreadManager::ThreadManager()
  : max_threads_(default_max_threads())
  , tbb_init_(tbb::task_scheduler_init::deferred)
{
#if HAVE_EIGEN
  // must be called before tbb threads are created via tbb::task_scheduler_init object ctor
  Eigen::initParallel();
  Eigen::setNbThreads(1);
#endif
  WITH_DUNE_FEM(Dune::Fem::ThreadManager::setMaxNumberThreads(1);)
}

#else // if HAVE_TBB

size_t Dune::Stuff::ThreadManager::max_threads()
{
  return 1;
}

size_t Dune::Stuff::ThreadManager::current_threads()
{
  return 1;
}

size_t Dune::Stuff::ThreadManager::thread()
{
  return 0;
}

void Dune::Stuff::ThreadManager::set_max_threads(const size_t count)
{
  if (count > 1)
    DUNE_THROW(InvalidStateException, "Trying to use more than one thread w/o TBB");
}

size_t Dune::Stuff::ThreadManager::default_max_threads()
{
  return 1;
}

Dune::Stuff::ThreadManager::ThreadManager()
  : max_threads_(1)
{
}

#endif // HAVE_DUNE_FEM
