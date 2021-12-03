// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2014, 2016 - 2020)
//   René Fritze     (2009 - 2016, 2018 - 2020)
//   Tobias Leibner  (2014, 2017 - 2018, 2020)

/**
 *  \file logging.hh
 *  \brief  logging
 **/
#ifndef DUNE_XT_COMMON_TIMEDLOGGING_HH
#define DUNE_XT_COMMON_TIMEDLOGGING_HH

#include <map>
#include <string>
#include <mutex>
#include <atomic>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/common/visibility.hh>

#include <dune/xt/common/color.hh>
#include <dune/xt/common/logstreams.hh>
#include <dune/xt/common/string.hh>
#include <utility>

#ifndef DUNE_XT_COMMON_TIMEDLOGGING_ENABLE_DEBUG
#  ifndef NDEBUG
#    define DUNE_XT_COMMON_TIMEDLOGGING_ENABLE_DEBUG 1
#  else
#    define DUNE_XT_COMMON_TIMEDLOGGING_ENABLE_DEBUG 0
#  endif
#endif

namespace Dune::XT::Common {


DUNE_EXPORT inline const Timer& SecondsSinceStartup()
{
  static Timer timer_;
  return timer_;
}


DUNE_EXPORT inline std::array<bool, 3>& default_logger_state()
{
  static std::array<bool, 3> state_{{false, false, true}};
  return state_;
}


/**
 * \brief A logging manager that provides info, debug and warning streams
 */
class DefaultLogger
{
  static std::string build_prefix(const std::string& prfx, const size_t cnt, const std::string& clr)
  {
    const std::string actual_color = terminal_supports_color() ? color(clr) : "";
    std::string copy_count_str;
    if (cnt > 0)
      copy_count_str += "[" + to_string(cnt) + "]";
    std::string ret;
    if (actual_color.empty())
      ret = prfx + copy_count_str + ": ";
    else
      ret = actual_color + StreamModifiers::bold + prfx + StreamModifiers::normal + copy_count_str
            + StreamModifiers::bold + ": " + StreamModifiers::normal;
    return ret;
  } // ... build_prefix(...)

public:
  std::string prefix;
  std::array<bool, 3> state;
  size_t copy_count{0};

  DefaultLogger(const std::string& prfx = "",
                const std::array<bool, 3>& initial_state = default_logger_state(),
                std::array<std::string, 3> colors = {{"blue", "darkgray", "red"}},
                bool global_timer = true);

  DefaultLogger(const DefaultLogger&);

  DefaultLogger& operator=(const DefaultLogger& other);

private:
  Timer timer_;
  std::array<std::string, 3> colors_;
  bool global_timer_;
  std::shared_ptr<std::ostream> info_;
  std::shared_ptr<std::ostream> debug_;
  std::shared_ptr<std::ostream> warn_;

public:
  bool info_enabled()
  {
    return state[0];
  }

  bool debug_enabled()
  {
    return state[1];
  }

  bool warn_enabled()
  {
    return state[2];
  }

  void enable(const std::string& prfx = "");

  void enable_like(const DefaultLogger& other);

  std::array<bool, 3> get_state_or(const std::array<bool, 3>& other_state) const;

  std::array<bool, 3> get_state_and(const std::array<bool, 3>& other_state) const;

  void state_or(const std::array<bool, 3>& other_state);

  void state_and(const std::array<bool, 3>& other_state);

  void disable();

  std::ostream& info();

  std::ostream& debug();

  std::ostream& warn();
}; // class DefaultLogger


/// \brief Convenience macro for enrivonments with a logger variable
#ifdef LOG
#  error Macro LOG already defined, open an issue at https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt/-/issues !
#else
#  define LOG(type)                                                                                                    \
    if (logger.type##_enabled())                                                                                       \
    logger.type()
#endif

/// \brief Convenience macro for classes with a logger member
#ifdef LOG_
#  error Macro LOG_ already defined, open an issue at https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt/-/issues !
#else
#  define LOG_(type)                                                                                                   \
    if (this->logger.type##_enabled())                                                                                 \
    this->logger.type()
#endif

/// \brief Convenience macro for classes with a logger member (allows to specify BaseType for diamond inheritance)
#ifdef LOG__
#  error Macro LOG_ already defined, open an issue at https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt/-/issues !
#else
#  define LOG__(base, type)                                                                                            \
    if (base ::logger.type##_enabled())                                                                                \
    base ::logger.type()
#endif


template <typename T = void>
class WithLogger
{
  using ThisType = WithLogger;

public:
  mutable DefaultLogger logger;

  WithLogger(const std::string& id, const std::array<bool, 3>& initial_state = {{true, true, true}})
    : logger(id, initial_state)
  {
    LOG_(debug) << "WithLogger(this=" << this << ")" << std::endl;
  }

  WithLogger(const ThisType& other)
    : logger(other.logger)
  {
    LOG_(debug) << "WithLogger(this=" << this << ", other=" << &other << ")" << std::endl;
  }

  WithLogger(ThisType&& source) noexcept
    : logger(std::move(source.logger))
  {
    LOG_(debug) << "WithLogger(this=" << this << ", source=" << &source << ")" << std::endl;
  }

  ~WithLogger()
  {
    LOG_(debug) << "~WithLogger(this=" << this << ")" << std::endl;
  }

  ThisType& operator=(const ThisType& other)
  {
    LOG_(debug) << "WithLogger.operator=(this=" << this << ", other=" << &other << ")" << std::endl;
  }

  ThisType& operator=(ThisType&& source) noexcept
  {
    LOG_(debug) << "WithLogger.operator=(this=" << this << ", source=" << &source << ")" << std::endl;
  }
}; // class WithLogger


/**
 * \brief A logging manager that provides info, debug and warning streams
 *
 * \note  Most likely you do not want to use this class directly but TimedLogger() instead.
 */
class TimedLogManager
{
public:
  TimedLogManager(const Timer& timer,
                  const std::string& info_prefix,
                  const std::string& debug_prefix,
                  const std::string& warning_prefix,
                  const ssize_t max_info_level,
                  const ssize_t max_debug_level,
                  const bool enable_warnings,
                  std::atomic<ssize_t>& current_level,
                  std::ostream& disabled_out = dev_null,
                  std::ostream& enabled_out = std::cout,
                  std::ostream& warn_out = std::cerr);

  ~TimedLogManager();

  std::ostream& info();

  std::ostream& debug();

  std::ostream& warn();

private:
  const Timer& timer_;
  std::atomic<ssize_t>& current_level_;
  std::shared_ptr<std::ostream> info_;
  std::shared_ptr<std::ostream> debug_;
  std::shared_ptr<std::ostream> warn_;
}; // class TimedLogManager


/**
 * \brief A logger that provides colored and prefixed streams.
 *
 * \note Most likely you do not want to use this class directly, but TimedLogger() instead.
 */
class TimedLogging
{
public:
  static constexpr ssize_t default_max_info_level = -1;
  static constexpr ssize_t default_max_debug_level = -1;
  static constexpr bool default_enable_warnings = true;
  static constexpr bool default_enable_colors = true;

  static std::string default_info_color()
  {
    return "blue";
  }
  static std::string default_debug_color()
  {
    return "darkgray";
  }
  static std::string default_warning_color()
  {
    return "red";
  }

  TimedLogging();

  /**
   * \brief sets the state
   *
   *        This method is mainly intended to be used on the global TimedLogger() instance. Before calling this method
   *        the state is set according to the defaults default_max_info_level, default_max_debug_level and
   *        default_enable_warnings.
   * \note  Calling this method more than once will throw an Exceptions::you_are_using_this_wrong, following the idea of
   *        least surprise.
   */
  void create(const ssize_t max_info_level = default_max_info_level,
              const ssize_t max_debug_level = default_max_debug_level,
              const bool enable_warnings = default_enable_warnings,
              const bool enable_colors = default_enable_colors,
              const std::string& info_color = default_info_color(),
              const std::string& debug_color = default_debug_color(),
              const std::string& warning_color = default_warning_color());

  TimedLogManager get(const std::string& id);

private:
  void update_colors();

  ssize_t max_info_level_;
  ssize_t max_debug_level_;
  bool enable_warnings_;
  bool enable_colors_;
  std::string info_prefix_;
  std::string debug_prefix_;
  std::string warning_prefix_;
  std::string info_suffix_;
  std::string debug_suffix_;
  std::string warning_suffix_;
  bool created_{false};
  std::atomic<ssize_t> current_level_;
  Timer timer_;
  std::mutex mutex_;
}; // class TimedLogging


/**
 * \brief Global instance of the timed logger.
 *
 *        This global logger instance is intended to be used in two ways:
 *        - Many classes or functions use this instance to log info, debug or warning messages. You can do so in your
 *          code by calling the TimedLogging::get() method, providing an identifier that should resemble the current
 *          scope:
\code
void user_function()
{
  TimedLogger().get("user_function").info() << "some information" << std::endl;
  for (size_t ii = 0; ii < 100; ++ii)
    TimedLogger().get("user_function").debug() << "debug output number " << ii << std::endl;
}
\endcode
            You can also hold a TimedLogManager object within the current scope or class, if wished:
\code
class UserClass
{
public:
  UserClass()
    : logger_(TimedLogger().get("UserClass"))
  {}

  void some_method()
  {
    logger_.warn() << "something is severly wrong!" << std::endl;
  }

private:
  TimedLogManager logger_;
}
\endcode
 *          Each time a new TimedLogManager is created using TimedLogging::get() the loglevel is increased, each time
 *          such a logger goes out of scope the loglevel is decreased.
 *        - You can use this instance to control the level (and style) of logging you want to have enabled in your
 *          application. You should call TimedLogging::create() as soon as possible (and only once!), until then all
 *          logging (execpt warnings) is disabled:
\code
void silent_function()
{
  auto logger = TimedLogger().get("silent_function");
  logger.info() << "This will never show!" << std::endl;
  const bool all_is_well = false;
  if (!all_is_well)
    logger.warn() << "But this warning will!" << std::endl;
}

int main()
{
  TimedLogger().create(0,   // max info level (only the first)
                       -1,  // max debug level (disabled)
                       true // warnings are enabled
                       );

  auto logger = TimedLogger().get("main");
  logger.info()  << "Welcome to my application!" << std::endl;
  logger.debug() << "This will never show!" << std::endl;

  silent_function();
}
\endcode
 *          In addition you can enable coloring of the streams (see TimedPrefixedLogStream) and give their respective
 *          colors, if wished (see the implementation of color_map() or the foreground colors of Colors for available
 *          colors):
\code
int main()
{
  TimedLogger().create(10,         // max info level
                       2,          // max debug level
                       true,       // warnings are enabled (the default)
                       true,       // colors are enabled (the default)
                       "blue",     // info color (the default)
                       "darkgrey", // debug color (the default)
                       "red"       // warning color (the default)
                      );

  auto logger = TimedLogger().get("main");
  logger.info() << "<- The 'main' prefix left of this should be blue!" << std::endl;
  logger.warn() << "<- The 'warn' prefix left of this should be red!"  << std::endl;
}
\endcode
 * \note Debug logging is only enabled if DUNE_XT_COMMON_TIMEDLOGGING_ENABLE_DEBUG is true (which is by default the case
 *       if NDEBUG is not defined) but you might still want to guard calls to logger.debug() for performance reasons.
 */
DUNE_EXPORT inline TimedLogging& TimedLogger()
{
  static TimedLogging timed_logger;
  return timed_logger;
}


/**
 * \brief Base class to generate debug output about object creation and deletion using the TimedLogging.
 * \sa    TimedLogging, TimedLogger
 * \note  Output is only generated in debug builds (or when DUNE_XT_COMMON_TIMEDLOGGING_ENABLE_DEBUG is true) and when
 *        the TimedLogger is active.
 *
 *        Simply derive from this class
\code
class Foo : EnableDebugLoggingForCtors<Foo>
{
public:
  Foo() : EnableDebugLoggingForCtors<Foo>("dune.xt.foo", "Foo") {}
};
\endcode
 *        and provide the prefix (for the logger) as well as the class name. Using your class, i.e. as in
\code
int main() {

  Foo foo;
  Foo bar(foo);
}
\endcode
 *        should yield an output like
\code
00:00|dune.xt.foo: Foo(this=0x0001)
00:00|dune.xt.foo: Foo(this=0x0002, other=0x0001)
00:00|dune.xt.foo: ~Foo(this=0x0002)
00:00|dune.xt.foo: ~Foo(this=0x0001)
\endcode
 * The same holds for the move ctor as well as move and assignment operators.
 */
template <typename T = void>
class NoOpEnableDebugLoggingForCtors
{
  using ThisType = NoOpEnableDebugLoggingForCtors<T>;

public:
  NoOpEnableDebugLoggingForCtors(const std::string&, const std::string&) {}

  NoOpEnableDebugLoggingForCtors(const ThisType& other) = default;

  NoOpEnableDebugLoggingForCtors(ThisType&& source) = default;

  ~NoOpEnableDebugLoggingForCtors() = default;

  ThisType& operator=(const ThisType& other) = default;

  ThisType& operator=(ThisType&& source) = default;

}; // class NoOpEnableDebugLoggingForCtors

template <typename T = void>
class ActiveEnableDebugLoggingForCtors
{
  using ThisType = ActiveEnableDebugLoggingForCtors<T>;

public:
  ActiveEnableDebugLoggingForCtors(const std::string& prefix, std::string class_id)
    : logger_(TimedLogger().get(prefix))
    , class_id_(std::move(class_id))
  {
    logger_.debug() << class_id_ << "(this=" << this << ")" << std::endl;
  }

  ActiveEnableDebugLoggingForCtors(const ThisType& other)
    : logger_(other.logger_)
    , class_id_(other.class_id_)
  {
    logger_.debug() << class_id_ << "(this=" << this << ", other=" << &other << ")" << std::endl;
  }

  ActiveEnableDebugLoggingForCtors(ThisType&& source) noexcept
    : logger_(std::move(source.logger_))
    , class_id_(std::move(source.class_id_))
  {
    logger_.debug() << class_id_ << "(this=" << this << ", source=" << &source << ")" << std::endl;
  }

  ~ActiveEnableDebugLoggingForCtors()
  {
    logger_.debug() << "~" << class_id_ << "(this=" << this << ")" << std::endl;
  }

  ThisType& operator=(const ThisType& other)
  {
    logger_.debug() << class_id_ << "operator=(this=" << this << ", other=" << &other << ")" << std::endl;
  }

  ThisType& operator=(ThisType&& source) noexcept
  {
    logger_.debug() << class_id_ << "operator=(this=" << this << ", source=" << &source << ")" << std::endl;
  }

protected:
  TimedLogManager logger_;

private:
  const std::string class_id_;
}; // class ActiveEnableDebugLoggingForCtors

#if DUNE_XT_COMMON_TIMEDLOGGING_ENABLE_DEBUG
template <class T>
using EnableDebugLoggingForCtors = ActiveEnableDebugLoggingForCtors<T>;
#else
template <class T>
using EnableDebugLoggingForCtors = NoOpEnableDebugLoggingForCtors<T>;
#endif

} // namespace Dune::XT::Common


#endif // DUNE_XT_COMMON_TIMED_LOGGING_HH
