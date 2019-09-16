// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk       (2018)

#include <dune/xt/common/test/main.hxx>

#include <dune/common/unused.hh>

#include <dune/xt/common/parallel/mpi_comm_wrapper.hh>

using namespace Dune;
using namespace Dune::XT::Common;

MPIHelper::MPICommunicator implicit_convert(MPI_Comm_Wrapper mpi_comm)
{
  return mpi_comm.get();
}

GTEST_TEST(MpiCommWrapper, Global)
{
  const auto orig = MPIHelper::getCommunicator();
  const auto convert = implicit_convert(orig);
  ASSERT_EQ(orig, convert);
}

GTEST_TEST(MpiCommWrapper, Local)
{
  const auto orig = MPIHelper::getLocalCommunicator();
  const auto convert = implicit_convert(orig);
  ASSERT_EQ(orig, convert);
}
