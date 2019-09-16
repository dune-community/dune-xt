// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014, 2016 - 2017)
//   René Fritze     (2012 - 2013, 2015 - 2016, 2018 - 2019)
//   Tobias Leibner  (2014, 2016, 2018 - 2019)

#include <dune/xt/test/main.hxx>

#include <dune/xt/common/parallel/communicator.hh>

GTEST_TEST(Communicator, All)
{
  const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
  EXPECT_TRUE(comm == comm);
  EXPECT_FALSE(comm != comm);
}
