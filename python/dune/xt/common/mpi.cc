// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2019)
//   René Fritze     (2018 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2018, 2020)

#include "config.h"

#include <python/dune/xt/common/bindings.hh>

#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/istl/owneroverlapcopy.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/exceptions.bindings.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/fmatrix.hh>
#include <python/dune/xt/common/configuration.hh>
#include <dune/xt/common/parallel/mpi_comm_wrapper.hh>
#include <dune/xt/common/parallel/helper.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/common/string.hh>


PYBIND11_MODULE(_common_mpi, m)
{
  using namespace Dune;
  namespace py = pybind11;
  using pybind11::operator""_a;

  m.def(
      "init_mpi",
      [](const std::vector<std::string>& args) {
        int argc = XT::Common::numeric_cast<int>(args.size());
        char** argv = XT::Common::vector_to_main_args(args);
        MPIHelper::instance(argc, argv);
      },
      "args"_a = std::vector<std::string>());

  py::class_<XT::Common::MPI_Comm_Wrapper> mpi_comm_wrapper(m, "MPI_Comm_Wrapper", "MPI_Comm_Wrapper");
  mpi_comm_wrapper.def(py::init());

  XT::Common::bindings::guarded_bind([&]() {
    using Comm = CollectiveCommunication<MPIHelper::MPICommunicator>;
    py::class_<Comm> cls(m, "CollectiveCommunication", "CollectiveCommunication");
    cls.def(py::init([](XT::Common::MPI_Comm_Wrapper comm) { return std::make_unique<Comm>(comm.get()); }),
            "comm"_a = XT::Common::MPI_Comm_Wrapper());

    cls.def_property_readonly("rank", &Comm::rank);
    cls.def_property_readonly("size", &Comm::size);

    cls.def("barrier", &Comm::barrier);

    cls.def(
        "min", [](const Comm& self, double x) { return self.min(x); }, "x"_a);
    cls.def(
        "max", [](const Comm& self, double x) { return self.max(x); }, "x"_a);
    cls.def(
        "sum", [](const Comm& self, double x) { return self.sum(x); }, "x"_a);
  });

  XT::Common::bindings::guarded_bind([&]() {
    py::class_<XT::SequentialCommunication> seqcomm(m, "SequentialCommunication", "SequentialCommunication");
    seqcomm.def(py::init());
  });

#if HAVE_MPI
  XT::Common::bindings::guarded_bind([&]() {
    py::class_<OwnerOverlapCopyCommunication<unsigned long, int>>(
        m, "OwnerOverlapCopyCommunication", "OwnerOverlapCopyCommunication");
  });
#endif

  m.def("abort_all_mpi_processes", &XT::abort_all_mpi_processes);
}
