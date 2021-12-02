// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017, 2020)
//   René Fritze     (2013 - 2016, 2018 - 2020)
//   Sven Kaulmann   (2014)
//   Tobias Leibner  (2016, 2019 - 2020)

#include "config.h"

#include <boost/algorithm/string/constants.hpp>
#include <boost/format.hpp>

#include <dune/xt/common/string.hh>

#include "logstreams.hh"

namespace Dune::XT::Common {


SuspendableStrBuffer::SuspendableStrBuffer(int loglevel, int& logflags)
  : logflags_(logflags)
  , loglevel_(loglevel)
  , suspended_logflags_(logflags)
  , is_suspended_(false)
  , suspend_priority_(default_suspend_priority)
{}

void SuspendableStrBuffer::suspend(PriorityType priority)
{
  suspend_priority_ = std::max(priority, suspend_priority_);
  {
    // don't accidentally invalidate flags if already suspended
    if (!is_suspended_) {
      suspended_logflags_ = logflags_;
      logflags_ = LOG_NONE;
    }
    is_suspended_ = true;
  }
} // Suspend

void SuspendableStrBuffer::resume(PriorityType priority)
{
  if (priority >= suspend_priority_) {
    if (is_suspended_)
      logflags_ = suspended_logflags_;
    is_suspended_ = false;
    suspend_priority_ = default_suspend_priority;
  }
} // Resume

std::streamsize SuspendableStrBuffer::xsputn(const char_type* s, std::streamsize count)
{
  if (enabled()) {
    std::lock_guard<std::mutex> guard(mutex_);
    return BaseType::xsputn(s, count);
  }
  // pretend everything was written
  return std::streamsize(count);
}

SuspendableStrBuffer::int_type SuspendableStrBuffer::overflow(SuspendableStrBuffer::int_type ch)
{
  if (enabled())
    return BaseType::overflow(ch);
  // anything not equal to traits::eof is considered a success
  return traits_type::eof() + 1;
}

int SuspendableStrBuffer::pubsync()
{
  if (enabled())
    return BaseType::pubsync();
  return 0;
}

TimedPrefixedStreamBuffer::TimedPrefixedStreamBuffer(const Timer& timer, const std::string& prefix, std::ostream& out)
  : timer_(timer)
  , prefix_(prefix)
  , out_(out)
  , prefix_needed_(true)
{}

int TimedPrefixedStreamBuffer::sync()
{
  [[maybe_unused]] std::lock_guard<std::mutex> guard(mutex_);
  const std::string tmp_str = str();
  if (prefix_needed_ && !tmp_str.empty()) {
    out_ << elapsed_time_str() << prefix_;
    prefix_needed_ = false;
  }
  auto lines = tokenize(tmp_str, "\n", boost::algorithm::token_compress_off);
  DXT_ASSERT(lines.size() > 0);
  out_ << lines[0];
  for (size_t ii = 1; ii < lines.size() - 1; ++ii)
    out_ << "\n" << elapsed_time_str() << prefix_ << lines[ii];
  if (lines.size() > 1) {
    out_ << "\n";
    const auto& last = lines.back();
    if (last.empty())
      prefix_needed_ = true;
    else
      out_ << elapsed_time_str() << prefix_ << last;
  }
  out_.flush();
  str("");
  return 0;
} // ... sync(...)

std::string TimedPrefixedStreamBuffer::elapsed_time_str() const
{
  const double secs_per_week = 604800;
  const double secs_per_day = 86400;
  const double secs_per_hour = 3600;
  const double elapsed = timer_.elapsed();
  const size_t weeks(elapsed / secs_per_week);
  const size_t days((elapsed - weeks * secs_per_week) / secs_per_day);
  const size_t hours((elapsed - weeks * secs_per_week - days * secs_per_day) / 3600.0);
  const size_t minutes((elapsed - weeks * secs_per_week - days * secs_per_day - hours * secs_per_hour) / 60.0);
  const size_t seconds(elapsed - weeks * secs_per_week - days * secs_per_day - hours * secs_per_hour - minutes * 60);
  if (elapsed > secs_per_week) // more than a week
    return (boost::format("%02dw %02dd %02d:%02d:%02d|") % weeks % days % hours % minutes % seconds).str();
  if (elapsed > secs_per_day) // less than a week, more than a day
    return (boost::format("%02dd %02d:%02d:%02d|") % days % hours % minutes % seconds).str();
  if (elapsed > secs_per_hour) // less than a day, more than one hour
    return (boost::format("%02d:%02d:%02d|") % hours % minutes % seconds).str();
  else // less than one hour
    return (boost::format("%02d:%02d|") % minutes % seconds).str();
} // ... elapsed_time(...)

LogStream& LogStream::flush()
{
  this->access().pubsync();
  return *this;
}

TimedPrefixedLogStream::TimedPrefixedLogStream(const Timer& timer, const std::string& prefix, std::ostream& outstream)
  : StorageBaseType(new TimedPrefixedStreamBuffer(timer, prefix, outstream))
  , OstreamBaseType(&this->access())
{}

TimedPrefixedLogStream::~TimedPrefixedLogStream()
{
  flush();
}

int OstreamBuffer::sync()
{
  // flush buffer into stream
  std::lock_guard<std::mutex> guard(sync_mutex_);
  out_ << str();
  out_.flush();
  str("");
  return 0;
}

int EmptyBuffer::sync()
{
  str("");
  return 0;
}

int CombinedBuffer::pubsync()
{
  if (!enabled())
    return 0;
  int ret = 0;
  for (auto&& buffer_ptr : buffer_)
    ret = buffer_ptr->pubsync();
  return ret;
}

std::streamsize CombinedBuffer::xsputn(const char_type* s, std::streamsize count)
{
  std::streamsize ret = 0;
  for (auto&& buffer_ptr : buffer_)
    ret = buffer_ptr->xsputn(s, count);
  return ret;
}

CombinedBuffer::int_type CombinedBuffer::overflow(int_type ch)
{
  int_type ret = 0;
  for (auto&& buffer_ptr : buffer_)
    ret = buffer_ptr->overflow(ch);
  return ret;
}

int CombinedBuffer::sync()
{
  int ret = 0;
  for (auto&& buffer_ptr : buffer_)
    ret = buffer_ptr->sync();
  return ret;
}

DualLogStream::DualLogStream(int loglevel, int& logflags, std::ostream& outstream, std::ofstream& file)
  : LogStream(new CombinedBuffer(
      loglevel,
      logflags,
      {new OstreamBuffer(loglevel, logflags, outstream), new OstreamBuffer(loglevel, logflags, file)}))
{}

OstreamLogStream::OstreamLogStream(int loglevel, int& logflags, std::ostream& outstream)
  : LogStream(new OstreamBuffer(loglevel, logflags, outstream))
{}

EmptyLogStream::EmptyLogStream(int& logflags)
  : LogStream(new EmptyBuffer(int(LOG_NONE), logflags))
{}


} // namespace Dune::XT::Common
