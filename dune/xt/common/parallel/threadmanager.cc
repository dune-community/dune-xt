// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   René Fritze     (2013 - 2016, 2018 - 2019)
//   Tobias Leibner  (2014 - 2016, 2018, 2020)

#include "config.h"

#include <atomic>
#include <thread>

#if HAVE_TBB
#  include <tbb/concurrent_unordered_map.h>
#endif

#include <boost/numeric/conversion/cast.hpp>
#include <boost/config.hpp>

#if HAVE_EIGEN
#  include <dune/xt/common/disable_warnings.hh>
#  include <Eigen/Core>
#  include <dune/xt/common/reenable_warnings.hh>
#endif

#include <dune/common/exceptions.hh>

#include <dune/xt/common/configuration.hh>

#include "threadmanager.hh"


#if HAVE_TBB

struct ThreadIdHashCompare
{
  size_t hash(std::thread::id id) const
  {
    static std::hash<std::thread::id> hasher;
    return hasher(id);
  }

  size_t operator()(std::thread::id id) const
  {
    return hash(id);
  }

  bool equal(std::thread::id a, std::thread::id b) const
  {
    return a == b;
  }
};

size_t Dune::XT::Common::ThreadManager::max_threads()
{
  const auto threads = DXTC_CONFIG_GET("threading.max_count", 1);
  return threads;
}

size_t Dune::XT::Common::ThreadManager::current_threads()
{
  const auto threads = max_threads();
  return threads;
}

template <typename Key, typename T, typename MapType>
std::pair<typename MapType::iterator, bool> tbb_map_emplace(MapType& map_in, Key key, T value)
{
#  if defined(BOOST_CLANG) && BOOST_CLANG
  return map_in.insert(typename MapType::value_type(key, value));
#  else
  return map_in.emplace(key, value);
#  endif
}

size_t Dune::XT::Common::ThreadManager::thread()
{
  const std::thread::id tbb_id = std::this_thread::get_id();
  static tbb::concurrent_unordered_map<std::thread::id, size_t, ThreadIdHashCompare> thread_ids;
  static std::atomic<size_t> map_size(0);
  const auto it = thread_ids.find(tbb_id);
  if (it == thread_ids.end())
    return tbb_map_emplace(thread_ids, tbb_id, map_size++).first->second;
  return it->second;
}

//! both std::hw_concur and intel's default_thread_count fail for mic
size_t Dune::XT::Common::ThreadManager::default_max_threads()
{
#  ifndef __MIC__
  return std::thread::hardware_concurrency();
#  else
  return DS_MAX_MIC_THREADS;
#  endif
}

void Dune::XT::Common::ThreadManager::set_max_threads(const size_t count)
{
  DXTC_CONFIG.set("threading.max_count", count, true);
  max_threads_ = count;
#  if HAVE_EIGEN
  Eigen::setNbThreads(boost::numeric_cast<int>(count));
#  endif
}

Dune::XT::Common::ThreadManager::ThreadManager()
  : max_threads_(default_max_threads())
{
#  if HAVE_EIGEN
  Eigen::initParallel();
  Eigen::setNbThreads(1);
#  endif
}

#else // if HAVE_TBB

size_t Dune::XT::Common::ThreadManager::max_threads()
{
  return max_threads_;
}

size_t Dune::XT::Common::ThreadManager::current_threads()
{
  return 1;
}

size_t Dune::XT::Common::ThreadManager::thread()
{
  return 0;
}

void Dune::XT::Common::ThreadManager::set_max_threads(const size_t count)
{
  if (count > 1)
    DUNE_THROW(InvalidStateException, "Trying to use more than one thread w/o TBB");
}

size_t Dune::XT::Common::ThreadManager::default_max_threads()
{
  return 1;
}

Dune::XT::Common::ThreadManager::ThreadManager()
  : max_threads_(1)
{}

#endif // HAVE_TBB
