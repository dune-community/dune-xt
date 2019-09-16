// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014 - 2017)
//   René Fritze     (2012 - 2013, 2015 - 2016, 2018)
//   Tobias Leibner  (2014, 2016)

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/common/logging.hh>
#include <dune/xt/common/logstreams.hh>

void balh(std::ostream& out)
{
  static int c = 0;
  out << "balh " << c << "\n";
  c++;
}

void do_something_that_takes_long(std::ostream& out)
{
  out << "  there should appear five dots, but not too fast:" << std::flush;
  for (size_t i = 0; i < 5; ++i) {
    busywait(666);
    out << "." << std::flush;
  }
  out << std::endl;
} // void do_something_that_takes_long()

GTEST_TEST(LoggerTest, all)
{
  using namespace Dune::XT::Common;
  Logger().create(LOG_CONSOLE | LOG_ERROR);
  Logger().error() << "This should be in output\n";
  Logger().info() << "This should NOT be in output\n";
  Logger().debug() << "dito\n";
  Logger().flush();
  for (int i : {LOG_INFO, LOG_DEBUG, LOG_ERROR}) {
    const int id = Logger().add_stream(LOG_CONSOLE | i);
    Logger().get_stream(id) << "Create a new stream with id: " << id << std::endl;
  }
  DXTC_LOG_ERROR.suspend();
  DXTC_LOG_ERROR << "not in output\n";
  balh(DXTC_LOG_ERROR);
  DXTC_LOG_ERROR.resume();
  DXTC_LOG_ERROR << "in output\n";
  balh(DXTC_LOG_ERROR);

  // this should do nothing whatsoever
  balh(dev_null);
  Logger().flush();

  // this is the desired result:
  LogStream& err = Logger().error();
  std::cout << "begin std::cout test" << std::endl;
  do_something_that_takes_long(std::cout);
  std::cout << "end   std::cout test" << std::endl;
  std::cout << "begin Logger().error() test" << std::endl;
  do_something_that_takes_long(err);
  std::cout << "end   Logger().error() test" << std::endl;
}

GTEST_TEST(LoggerTest, file)
{
  using namespace Dune::XT::Common;
  Logger().create(LOG_INFO | LOG_CONSOLE | LOG_FILE, "test_common_logger", "", "");
  Logger().info() << "This output should be in 'test_common_logger.log'" << std::endl;
}
