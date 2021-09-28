// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   René Fritze     (2014 - 2016, 2018 - 2019)
//   Tobias Leibner  (2016, 2020)

#ifndef DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS
#  define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 0
#endif

#include "config.h"

#include <dune/xt/test/gtest/gtest.h>
#include <dune/xt/common/timedlogging.hh>

#include <dune/xt/test/common.hh>

using namespace Dune;
using namespace Dune::XT;
using namespace Dune::XT::Common;

void before_create()
{
  TimedLogger().get("before_create").info() << "this info should not be visible" << std::endl;
  TimedLogger().get("before_create").debug() << "this debug should not be visible" << std::endl;
  TimedLogger().get("before_create").warn() << "this warning should be visible in red" << std::endl;
}

void after_create_inner()
{
  TimedLogger().get("after_create_inner").info() << "this info should be visible in blue" << std::endl;
  TimedLogger().get("after_create_inner").debug() << "this debug should not be visible" << std::endl;
  TimedLogger().get("after_create_inner").warn() << "this warning should not be visible" << std::endl;
}

void after_create()
{
  auto logger = TimedLogger().get("after_create");
  logger.info() << "this info should be visible in blue" << std::endl;
  logger.debug() << "this debug should be visible in yellow" << std::endl;
  logger.warn() << "this warning should not be visible" << std::endl;
  after_create_inner();
}

void fool_level_tracking_inner()
{
  TimedLogger().get("fool_level_tracking_inner").info() << "this info should be visible in blue" << std::endl;
  TimedLogger().get("fool_level_tracking_inner").debug() << "this debug should be visible in yellow" << std::endl;
  TimedLogger().get("fool_level_tracking_inner").warn() << "this warning should not be visible" << std::endl;
}

void fool_level_tracking()
{
  TimedLogger().get("fool_level_tracking").info() << "this info should be visible in blue" << std::endl;
  TimedLogger().get("fool_level_tracking").debug() << "this debug should be visible in yellow" << std::endl;
  TimedLogger().get("fool_level_tracking").warn() << "this warning should not be visible" << std::endl;
  fool_level_tracking_inner();
}

GTEST_TEST(TimedPrefixedLogStream, all)
{
  Timer timer;
  TimedPrefixedLogStream out(timer, "prefix: ", std::cout);
  out << "sample\nline" << std::flush;
  busywait(2000);
  out << "\n" << 3 << "\n\nend" << std::endl;
} // TEST(TimedPrefixedLogStream, all)

GTEST_TEST(TimedLogger, before_create)
{
  auto logger = TimedLogger().get("main");
  auto& info = logger.info();
  info << "this info should be visible in " << TimedLogging::default_info_color() << std::endl;
  logger.debug() << "this debug should not be visible" << std::endl;
  logger.warn() << "this warning should be visible in " << TimedLogging::default_warning_color() << std::endl;
  before_create();
}

GTEST_TEST(TimedLogger, after_create)
{
  TimedLogger().create(10, 1, false, true, "blue", "yellow");
  auto logger = TimedLogger().get("main");
  logger.info() << "this info should be visible in blue" << std::endl;
  logger.debug() << "this debug should be visible in yellow" << std::endl;
  logger.warn() << "this warning should not be visible" << std::endl;
  after_create();
}

GTEST_TEST(TimedLogger, fool_level_tracking)
{
  auto logger = TimedLogger().get("");
  logger.info() << "this info should be visible in blue" << std::endl;
  logger.debug() << "this debug should be visible in yellow" << std::endl;
  logger.warn() << "this warning should not be visible" << std::endl;
  fool_level_tracking();
}

int main(int argc, char** argv)
{
#if DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS
  try {
#endif
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
#if DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS
  } catch (Dune::Exception& e) {
    std::cerr << "\nDune reported error: " << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << "\n" << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    std::abort();
  } // try
#endif // DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS
} // ... main(...)
