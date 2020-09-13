// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   René Fritze     (2012 - 2016, 2018 - 2019)
//   Sven Kaulmann   (2013)
//   Tobias Leibner  (2014, 2018, 2020)

#include "config.h"

#include "timings.hh"

#if HAVE_LIKWID && ENABLE_PERFMON
#  include <likwid.h>
#  define DXTC_LIKWID_BEGIN_SECTION(name) LIKWID_MARKER_START(name.c_str());
#  define DXTC_LIKWID_END_SECTION(name) LIKWID_MARKER_STOP(name.c_str());
#  define DXTC_LIKWID_INIT LIKWID_MARKER_INIT
#  define DXTC_LIKWID_CLOSE LIKWID_MARKER_CLOSE
#else
#  define DXTC_LIKWID_BEGIN_SECTION(name)
#  define DXTC_LIKWID_END_SECTION(name)
#  define DXTC_LIKWID_INIT
#  define DXTC_LIKWID_CLOSE
#endif

// #include <dune/xt/common/string.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/filesystem.hh>

#include <dune/xt/common/disable_warnings.hh>
#include <boost/format.hpp>
#include <dune/xt/common/reenable_warnings.hh>

namespace Dune {
namespace XT {
namespace Common {


TimingData::TimingData(std::string _name)
  : timer_(new boost::timer::cpu_timer)
  , name(_name)
{
  timer_->start();
}

void TimingData::stop()
{
  timer_->stop();
}

TimingData::DeltaType TimingData::delta() const
{
  const auto scale = 1.0 / double(boost::timer::nanosecond_type(1e6));
  const auto elapsed = timer_->elapsed();
  const auto cast = [=](double var) { return static_cast<typename TimingData::DeltaType::value_type>(var); };
  return {{cast(elapsed.wall * scale), cast(elapsed.user * scale), cast(elapsed.system * scale)}};
}

void Timings::reset(std::string section_name)
{
  try {
    stop(section_name);
  } catch (Dune::RangeError&) {
    // ok, timer simply wasn't running
  }
  commited_deltas_[section_name] = {{0, 0, 0}};
}

void Timings::start(std::string section_name)
{
  std::lock_guard<std::mutex> lock(mutex_);

  const KnownTimersMap::iterator section = known_timers_map_.find(section_name);
  if (section != known_timers_map_.end()) {
    if (section->second.first) // timer currently running
      return;
    section->second.first = true; // set active, start with new
    section->second.second = TimingData(section_name);
  } else {
    // init new section
    known_timers_map_[section_name] = std::make_pair(true, TimingData(section_name));
  }
  DXTC_LIKWID_BEGIN_SECTION(section_name)
} // StartTiming

long Timings::stop(std::string section_name)
{
  DXTC_LIKWID_END_SECTION(section_name)
  if (known_timers_map_.find(section_name) == known_timers_map_.end())
    DUNE_THROW(Dune::RangeError, "trying to stop timer " << section_name << " that wasn't started\n");

  known_timers_map_[section_name].first = false; // marks as not running
  TimingData& timing = *(known_timers_map_[section_name].second);
  timing.stop();
  const auto dlt = timing.delta();
  if (commited_deltas_.find(section_name) == commited_deltas_.end())
    commited_deltas_[section_name] = dlt;
  else {
    for (auto i : value_range(dlt.size()))
      commited_deltas_[section_name][i] += dlt[i];
  }
  return dlt[0];
} // StopTiming

TimingData::TimeType Timings::walltime(std::string section_name) const
{
  return delta(section_name)[0];
}

TimingData::DeltaType Timings::delta(std::string section_name) const
{
  DeltaMap::const_iterator section = commited_deltas_.find(section_name);
  if (section == commited_deltas_.end()) {
    // timer might still be running
    const auto& timer_it = known_timers_map_.find(section_name);
    if (timer_it == known_timers_map_.end())
      DUNE_THROW(Dune::InvalidStateException, "no timer found: " + section_name);
    return timer_it->second.second->delta();
  }
  return section->second;
}

void Timings::stop()
{
  for (auto&& section : known_timers_map_) {
    try {
      stop(section.first);
    } catch (Dune::RangeError&) {
    }
  }
} // GetTiming

void Timings::reset()
{
  stop();
  commited_deltas_.clear();
} // Reset

void Timings::set_outputdir(std::string dir)
{
  output_dir_ = dir;
  test_create_directory(output_dir_);
}

void Timings::output_per_rank(std::string csv_base) const
{
  const auto rank = MPIHelper::getCollectiveCommunication().rank();
  boost::filesystem::path dir(output_dir_);
  boost::filesystem::path filename = dir / (boost::format("%s_p%08d.csv") % csv_base % rank).str();
  boost::filesystem::ofstream out(filename);
  output_all_measures(out, MPIHelper::getLocalCommunicator());
  std::stringstream tmp_out;
  output_all_measures(tmp_out, MPIHelper::getCommunicator());
  if (rank == 0) {
    boost::filesystem::path a_filename = dir / (boost::format("%s.csv") % csv_base).str();
    boost::filesystem::ofstream a_out(a_filename);
    a_out << tmp_out.str() << std::endl;
  }
}

void Timings::output_simple(std::ostream& out) const
{
  for (const auto& section : commited_deltas_) {
    out << csv_sep_ << section.first;
  }
  for (const auto& section : commited_deltas_) {
    out << csv_sep_ << section.second[0];
    ;
  }
  out << std::endl;
}

void Timings::output_all_measures(std::ostream& out, MPIHelper::MPICommunicator mpi_comm) const
{
  CollectiveCommunication<MPIHelper::MPICommunicator> comm(mpi_comm);
  std::stringstream stash;

  stash << "threads" << csv_sep_ << "ranks";
  for (const auto& section : commited_deltas_) {
    stash << csv_sep_ << section.first << "_avg_usr" << csv_sep_ << section.first << "_max_usr" << csv_sep_
          << section.first << "_avg_wall" << csv_sep_ << section.first << "_max_wall" << csv_sep_ << section.first
          << "_avg_sys" << csv_sep_ << section.first << "_max_sys";
  }
  const auto weight = 1 / double(comm.size());

  stash << std::endl << threadManager().max_threads() << csv_sep_ << comm.size();
  for (const auto& section : commited_deltas_) {
    const auto timings = section.second;
    auto wall = timings[0];
    auto usr = timings[1];
    auto sys = timings[2];
    const auto wall_sum = comm.sum(wall);
    const auto wall_max = comm.max(wall);
    const auto usr_sum = comm.sum(usr);
    const auto usr_max = comm.max(usr);
    const auto sys_sum = comm.sum(sys);
    const auto sys_max = comm.max(sys);
    stash << csv_sep_ << usr_sum * weight << csv_sep_ << usr_max << csv_sep_ << wall_sum * weight << csv_sep_
          << wall_max << csv_sep_ << sys_sum * weight << csv_sep_ << sys_max;
  }

  stash << std::endl;
  if (comm.rank() == 0)
    out << stash.str();
}

Timings::Timings()
  : csv_sep_(",")
{
  DXTC_LIKWID_INIT;
  reset();
  set_outputdir("profiling");
}

Timings::~Timings()
{
  DXTC_LIKWID_CLOSE;
}

OutputScopedTiming::OutputScopedTiming(const std::string& section_name, std::ostream& out)
  : ScopedTiming(section_name)
  , out_(out)
{}

OutputScopedTiming::~OutputScopedTiming()
{
  const auto duration = timings().stop(section_name_);
  out_ << "Executing " << section_name_ << " took " << duration / 1000.f << "s\n";
}


} // namespace Common
} // namespace XT
} // namespace Dune
