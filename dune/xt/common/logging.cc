// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2014, 2016 - 2017)
//   René Fritze     (2012 - 2016, 2018 - 2019)
//   Tobias Leibner  (2020)

#include "config.h"

#include <boost/format.hpp>

// #include <dune/xt/common/memory.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/filesystem.hh>

#include "logging.hh"

namespace Dune {
namespace XT {
namespace Common {

Logging::Logging()
  : streamIDs_({LOG_ERROR, LOG_DEBUG, LOG_INFO})
  , logflags_(LOG_NONE)
  , emptyLogStream_(logflags_)
{
  for (const auto id : streamIDs_)
    streammap_[id] = std::make_unique<EmptyLogStream>(logflags_);
}

void Logging::deinit()
{
  streammap_.clear();
  if ((logflags_ & LOG_FILE) != 0) {
    logfile_ << std::endl;
    logfile_.close();
  }
}

Logging::~Logging()
{
  deinit();
}

void Logging::create(int logflags, const std::string logfile, const std::string datadir, const std::string _logdir)
{
  using namespace boost::filesystem;
  const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
  boost::format log_fn("%s%s");
  if (comm.size() > 1) {
    const std::string rank = (boost::format("%08d") % comm.rank()).str();
    log_fn = boost::format("%s_p" + rank + "_%s");
  }
  logflags_ = logflags;
  path logdir = path(datadir) / _logdir;
  filename_ = logdir / (log_fn % logfile % ".log").str();
  test_create_directory(filename_.string());
  const bool file_logging = ((logflags_ & LOG_FILE) != 0);
  if (file_logging) {
    logfile_.open(filename_);
    DXT_ASSERT(logfile_.is_open());
  }

  for (const auto id : streamIDs_) {
    flagmap_[id] = logflags;
    streammap_[id] = std::make_unique<DualLogStream>(id, flagmap_[id], std::cout, logfile_);
  }
} // create

void Logging::set_prefix(std::string prefix)
{
  deinit();
  create(logflags_, prefix);
} // set_prefix

void Logging::set_stream_flags(int streamID, int flags)
{
  DXT_ASSERT(flagmap_.find(streamID) != flagmap_.end());
  // this might result in logging to diff targtes, so we flush the current targets
  flush();
  flagmap_[streamID] = flags;
}

int Logging::get_stream_flags(int streamID) const
{
  const auto it = flagmap_.find(streamID);
  if (it == flagmap_.end())
    DUNE_THROW(InvalidStateException,
               "cannot get flags for unkown Stream id, maybe you forgot to call \"Logging::create\"?");
  return it->second;
}

LogStream& Logging::get_stream(int streamId)
{
  const auto it = streammap_.find(streamId);
  if (it == streammap_.end())
    DUNE_THROW(InvalidStateException, "cannot get unknown logStream, maybe you forgot to call \"Logging::create\"?");
  else {
    DXT_ASSERT(it->second);
    return *(it->second);
  }
}

void Logging::flush()
{
  for (auto& pair : streammap_) {
    DXT_ASSERT(pair.second);
    pair.second->flush();
  }
} // flush

int Logging::add_stream(int flags)
{
  static int streamID_int = LOG_NEXT;
  streamID_int <<= 1;
  int streamID = streamID_int;
  streamIDs_.push_back(streamID);
  flagmap_[streamID] = (flags | streamID);
  streammap_[streamID] = std::make_unique<DualLogStream>(streamID, flagmap_[streamID], std::cout, logfile_);
  return streamID_int;
} // add_stream

void Logging::resume(LogStream::PriorityType prio)
{
  for (auto& pair : streammap_) {
    pair.second->resume(prio);
  }
} // resume

void Logging::suspend(LogStream::PriorityType prio)
{
  for (auto& pair : streammap_) {
    pair.second->suspend(prio);
  }
} // suspend

} // namespace Common
} // namespace XT
} // namespace Dune
