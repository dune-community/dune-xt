// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2017 - 2019)
//   Tobias Leibner  (2017, 2020)

#include <config.h>

#include "memory.hh"

#include <dune/xt/common/timings.hh>
#include <dune/xt/common/filesystem.hh>
#include <dune/xt/common/configuration.hh>

#include <sys/resource.h>

namespace Dune {
namespace XT {
namespace Common {

void mem_usage(std::string filename)
{
  auto comm = Dune::MPIHelper::getCollectiveCommunication();
  // Compute the peak memory consumption of each processes
  int who = RUSAGE_SELF;
  struct rusage usage;
  getrusage(who, &usage);
  long peakMemConsumption = usage.ru_maxrss;
  // compute the maximum and mean peak memory consumption over all processes
  long maxPeakMemConsumption = comm.max(peakMemConsumption);
  long totalPeakMemConsumption = comm.sum(peakMemConsumption);
  long meanPeakMemConsumption = totalPeakMemConsumption / comm.size();
  // write output on rank zero
  if (comm.rank() == 0) {
    std::unique_ptr<boost::filesystem::ofstream> memoryConsFile(make_ofstream(filename));
    *memoryConsFile << "global.maxPeakMemoryConsumption,global.meanPeakMemoryConsumption\n"
                    << maxPeakMemConsumption << "," << meanPeakMemConsumption << std::endl;
  }
}

void mem_usage()
{
  mem_usage(std::string(DXTC_CONFIG_GET("global.datadir", "data/")) + std::string("/memory.csv"));
}

} // namespace Common
} // namespace XT
} // namespace Dune
