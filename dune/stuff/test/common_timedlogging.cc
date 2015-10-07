// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS
#define DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS 0
#endif

#include "config.h"

#include <dune/stuff/test/gtest/gtest.h>
#include <dune/stuff/common/timedlogging.hh>

#include "common.hh"

using namespace Dune;
using namespace Dune::Stuff;
using namespace Dune::Stuff::Common;

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

TEST(TimedPrefixedLogStream, all)
{
  Timer timer;
  TimedPrefixedLogStream out(timer, "prefix: ", std::cout);
  out << "sample\nline" << std::flush;
  busywait(2000);
  out << "\n" << 3 << "\n\nend" << std::endl;
} // TEST(TimedPrefixedLogStream, all)

TEST(TimedLogger, before_create)
{
  auto logger = TimedLogger().get("main");
  auto& info = logger.info();
  info << "this info should be visible in " << TimedLogging::default_info_color() << std::endl;
  logger.debug() << "this debug should not be visible" << std::endl;
  logger.warn() << "this warning should be visible in " << TimedLogging::default_warning_color() << std::endl;
  before_create();
}

TEST(TimedLogger, after_create)
{
  TimedLogger().create(10, 1, false, true, "blue", "yellow");
  auto logger = TimedLogger().get("main");
  logger.info() << "this info should be visible in blue" << std::endl;
  logger.debug() << "this debug should be visible in yellow" << std::endl;
  logger.warn() << "this warning should not be visible" << std::endl;
  after_create();
}

TEST(TimedLogger, fool_level_tracking)
{
  auto logger = TimedLogger().get("");
  logger.info() << "this info should be visible in blue" << std::endl;
  logger.debug() << "this debug should be visible in yellow" << std::endl;
  logger.warn() << "this warning should not be visible" << std::endl;
  fool_level_tracking();
}

int main(int argc, char** argv)
{
#if DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS
  try {
#endif
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
#if DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS
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
#endif // DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS
} // ... main(...)
