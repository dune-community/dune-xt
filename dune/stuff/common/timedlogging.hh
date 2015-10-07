// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

/**
   *  \file logging.hh
   *  \brief  logging
   **/
#ifndef DUNE_STUFF_COMMON_TIMEDLOGGING_HH
#define DUNE_STUFF_COMMON_TIMEDLOGGING_HH

#include <map>
#include <string>
#include <mutex>
#include <atomic>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dune/stuff/common/logstreams.hh>
#include <dune/stuff/common/color.hh>

namespace Dune {
namespace Stuff {
namespace Common {

/**
 * \brief A logging manager that provides info, debug and warning streams
 *
 * \note  Most likely you do not want to use this class directly but TimedLogger() instead.
 */
class TimedLogManager
{
public:
  TimedLogManager(const Timer& timer, const std::string info_prefix, const std::string debug_prefix,
                  const std::string warning_prefix, const ssize_t max_info_level, const ssize_t max_debug_level,
                  const bool enable_warnings, std::atomic<ssize_t>& current_level,
                  std::ostream& disabled_out = dev_null, std::ostream& enabled_out = std::cout,
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
  static const ssize_t default_max_info_level  = -1;
  static const ssize_t default_max_debug_level = -1;
  static const bool default_enable_warnings    = true;
  static const bool default_enable_colors = true;
  static const std::string default_info_color()
  {
    return "white";
  }
  static const std::string default_debug_color()
  {
    return "darkgray";
  }
  static const std::string default_warning_color()
  {
    return "red";
  }

  TimedLogging();

  /**
   * \brief sets the state
   *
   *        This methos is mainly intended to be used on the global TimedLogger() instance. Before calling this method
   *        the state is set according to the defaults default_max_info_level, default_max_debug_level and
   *        default_enable_warnings.
   * \note  Calling this method more than once will throw an Exceptions::you_are_using_this_wrong, following the idea of
   *        least surprise.
   */
  void create(const ssize_t max_info_level = default_max_info_level,
              const ssize_t max_debug_level = default_max_debug_level,
              const bool enable_warnings = default_enable_warnings, const bool enable_colors = default_enable_colors,
              const std::string info_color    = default_info_color(),
              const std::string debug_color   = default_debug_color(),
              const std::string warning_color = default_warning_color());

  TimedLogManager get(const std::string id);

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
  bool created_;
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
  TimedLogger().create(10,          // max info level
                       2,           // max debug level
                       true,        // warnings are enabled (the default)
                       true,        // colors are enabled (the default)
                       "white",     // info color (the default)
                       "lightgrey", // debug color (the default)
                       "red"        // warning color (the default)
                      );

  auto logger = TimedLogger().get("main");
  logger.info() << "<- The 'main' prefix left of this should be white!" << std::endl;
  logger.warn() << "<- The 'warn' prefix left of this should be red!"   << std::endl;
}
\endcode
 * \note Debug logging is only enabled if NDEBUG is not defined but you might still want to guard calls to
 *       logger.debug() for performance reasons.
 */
TimedLogging& TimedLogger();

} // namespace Common
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_COMMON_TIMED_LOGGING_HH
