// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014 - 2017)
//   René Fritze     (2009 - 2016, 2018 - 2019)
//   Sven Kaulmann   (2011)
//   Tobias Leibner  (2018, 2020)

#ifndef DUNE_XT_COMMON_PROFILER_HH
#define DUNE_XT_COMMON_PROFILER_HH

#ifndef DUNE_XT_COMMON_DO_PROFILE
#  define DUNE_XT_COMMON_DO_PROFILE 0
#endif

#ifndef DUNE_XT_COMMON_DO_TIMING
#  define DUNE_XT_COMMON_DO_TIMING 0
#endif

#include <array>
#include <atomic>
#include <ctime>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/timer/timer.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/xt/common/parallel/threadstorage.hh>

namespace Dune::XT::Common {


class Timings;

//! XT::Profiler global instance
Timings& timings();

//! wraps name, start- and end time for one timing section
struct TimingData
{
private:
  std::shared_ptr<boost::timer::cpu_timer> timer_;

public:
  std::string name;

  explicit TimingData(std::string _name = "blank");

  void stop();

  using TimeType = boost::timer::nanosecond_type;
  using DeltaType = std::array<TimeType, 3>;

  /** \return array{wall,user,sys}: times elapsed since object construction in
   *milliseconds
   *  \note since typical resolutions for user+system time are 10-15ms the
   *nanosecond results are scaled accordingly
   **/
  DeltaType delta() const;
};

//! a utility class to time a limited scope of code
class ScopedTiming;

/** \brief simple inline timing class
 *  - User can set as many (even nested) named sections whose total (=system+user) time will be computed across all
 *    program instances.\n
 *  - Provides csv-conform output of process-averaged runtimes.
 **/
class Timings
{
  friend Timings& timings();

private:
  Timings();

  using KnownTimersMap = std::map<std::string, std::pair<std::atomic<bool>, UnsafePerThreadValue<TimingData>>>;
  //! section name -> seconds
  using DeltaMap = std::map<std::string, TimingData::DeltaType>;

  const TimingData& get_timing_data(std::string section_name) const;

public:
  ~Timings();

  void stop();

  //! set this to begin a named section
  void start(std::string section_name);

  //! stop named section's counter
  long stop(std::string section_name);

  //! set elapsed time back to 0 for section_name
  void reset(std::string section_name);

  //! get runtime of section in milliseconds
  TimingData::TimeType walltime(std::string section_name) const;
  //! get the full delta array
  TimingData::DeltaType delta(std::string section_name) const;

  /** creates one file local to each MPI-rank (no global averaging)
   *  one single rank-0 file with all combined/averaged measures
   ***/
  void output_per_rank(std::string csv_base) const;
  //! outputs walltime only w/o MPI-rank averaging
  void output_simple(std::ostream& out = std::cout) const;
  /** output all recorded measures
   * \note outputs average, min, max over all MPI processes associated to mpi_comm **/
  void output_all_measures(std::ostream& out = std::cout,
                           MPIHelper::MPICommunicator mpi_comm = Dune::MPIHelper::getCommunicator()) const;

  /// stops and resets all timers and data
  void reset();

  void set_outputdir(std::string dir);

private:
  DeltaMap commited_deltas_;
  //! runtime tables etc go there
  std::string output_dir_;

  KnownTimersMap known_timers_map_;
  const std::string csv_sep_;
  std::mutex mutex_;
};

//! global profiler object
DUNE_EXPORT inline Timings& timings()
{
  static Timings pf;
  return pf;
}

class ScopedTiming : public boost::noncopyable
{
protected:
  const std::string section_name_;

public:
  explicit inline ScopedTiming(const std::string& section_name)
    : section_name_(section_name)
  {
    timings().start(section_name_);
  }

  inline ~ScopedTiming()
  {
    timings().stop(section_name_);
  }
};

struct OutputScopedTiming : public ScopedTiming
{
  OutputScopedTiming(const std::string& section_name, std::ostream& out);

  ~OutputScopedTiming();

protected:
  std::ostream& out_;
};


} // namespace Dune::XT::Common

#define DXTC_TIMINGS Dune::XT::Common::timings()

#if DUNE_XT_COMMON_DO_TIMING
#  define DUNE_XT_COMMON_TIMING_SCOPE(section_name) [[maybe_unused]] Dune::XT::Common::ScopedTiming timer(section_name)
#else
#  define DUNE_XT_COMMON_TIMING_SCOPE(section_name)
#endif

#endif // DUNE_XT_COMMON_PROFILER_HH
