// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2019)

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/common/parallel/communicator.hh>

GTEST_TEST(Communicator, All)
{
  const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
  EXPECT_TRUE(comm == comm);
  EXPECT_FALSE(comm != comm);
}
