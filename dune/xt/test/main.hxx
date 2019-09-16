// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2018)
//   René Fritze     (2014 - 2016, 2018 - 2019)
//   Tobias Leibner  (2015 - 2018)

#include "config.h"

#ifndef DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS
#  define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 0
#endif
#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_DEBUG_LOGGING
#  define DUNE_XT_COMMON_TEST_MAIN_ENABLE_DEBUG_LOGGING 0
#endif
#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING
#  define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 0
#endif
#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING
#  define DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING 0
#endif

#define GTEST_DONT_DEFINE_TEST 1

#include <string>
#include <vector>
#include <map>
#include <random>
#include <fstream>
#include <limits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/xt/common/math.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/vector.hh>
#include <dune/xt/test/gtest/gtest.h>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/common/convergence-study.hh>
#include <dune/xt/common/parallel/threadmanager.hh>


#if HAVE_TBB
#  include <tbb/task_scheduler_init.h>
#endif

#include <dune/xt/test/common.hh>

int main(int argc, char** argv)
{
  using namespace Dune;
  using namespace Dune::XT::Common;
#if DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS
  try {
#endif

    testing::InitGoogleTest(&argc, argv);
    if (argc > 1)
      DXTC_CONFIG.read_command_line(argc, argv);

    MPIHelper::instance(argc, argv);

    Logger().create(
#if DUNE_XT_COMMON_TEST_MAIN_ENABLE_DEBUG_LOGGING
        LOG_CONSOLE | LOG_INFO | LOG_DEBUG | LOG_ERROR
#elif DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING
      LOG_CONSOLE | LOG_INFO | LOG_ERROR
#else
      LOG_CONSOLE | LOG_ERROR
#endif
        ,
        "",
        "",
        "");

    const ssize_t max_info_level = DXTC_CONFIG.get("timedlogging.max_info",
#if DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING && DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING
                                                   std::numeric_limits<ssize_t>::max()
#else
                                                 -1
#endif
    );
    const ssize_t max_debug_level = DXTC_CONFIG.get("timedlogging.max_debug",
#if DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING && DUNE_XT_COMMON_TEST_MAIN_ENABLE_DEBUG_LOGGING
                                                    std::numeric_limits<ssize_t>::max()
#else
                                                  -1
#endif
    );
    TimedLogger().create(max_info_level /*info*/, max_debug_level /*debug*/);

    const size_t threads = DXTC_CONFIG.has_key("threading.max_count") // <- doing this so complicated to
                               ? DXTC_CONFIG.get<size_t>("threading.max_count") //    silence the WARNING: ...
                               : 1;
#if HAVE_TBB
    tbb::task_scheduler_init tbb_init(boost::numeric_cast<int>(threads));
#endif
    threadManager().set_max_threads(threads);

    return RUN_ALL_TESTS();

#if DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS
  } catch (Dune::Exception& e) {
    std::cerr << "\nDune reported error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (std::exception& e) {
    std::cerr << "\n" << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return EXIT_FAILURE;
  } // try
#endif // DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS
} // ... main(...)
