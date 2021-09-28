// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   René Fritze    (2018 - 2020)
//   Tobias Leibner (2020)

#include <config.h>
#include "dune/xt/common/parallel/mpi_comm_wrapper.hh"


namespace Dune::XT::Common {

MPI_Comm_Wrapper::MPI_Comm_Wrapper(WrappedComm comm)
  : comm_(comm)
{}

MPI_Comm_Wrapper& MPI_Comm_Wrapper::operator=(const WrappedComm comm)
{
  this->comm_ = comm;
  return *this;
}

MPI_Comm_Wrapper::WrappedComm MPI_Comm_Wrapper::get() const
{
  return comm_;
}


} // namespace Dune::XT::Common