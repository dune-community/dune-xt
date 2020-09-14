// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2017)
//   René Fritze     (2012 - 2013, 2015 - 2016, 2018 - 2019)
//   Stefan Girke    (2012)
//   Tobias Leibner  (2014, 2020)

#ifndef DUNE_XT_COMMON_COLOR_HH
#define DUNE_XT_COMMON_COLOR_HH

#include <map>
#include <string>

namespace Dune {
namespace XT {
namespace Common {

/**
 * @brief namespace to define color constants that can be
 * used to print colored text in an output stream.
 *
 * * \todo this could go size_to libdunextcommon
 * @warning Some color codes might be unsupported by your terminal.
 */
struct Colors
{
  // foreground colors
  static const char* black;
  static const char* red;
  static const char* green;
  static const char* brown;
  static const char* blue;
  static const char* purple;
  static const char* cyan;
  static const char* lightgray;
  // light foreground colors
  static const char* darkgray;
  static const char* lightred;
  static const char* lightgreen;
  static const char* yellow;
  static const char* lightblue;
  static const char* lightpurple;
  static const char* lightcyan;
  static const char* white;

  // background colors
  static const char* bblack;
  static const char* bred;
  static const char* bgreen;
  static const char* bbrown;
  static const char* bblue;
  static const char* bpurple;
  static const char* bcyan;
  static const char* blightgray;
  // light background colors
  static const char* bdarkgray;
  static const char* blightred;
  static const char* blightgreen;
  static const char* byellow;
  static const char* blightblue;
  static const char* blightpurple;
  static const char* blightcyan;
  static const char* bwhite;
}; // struct Colors

// modifiers
struct StreamModifiers
{
  static const char* normal;
  static const char* bold;
  static const char* italic;
  static const char* underline;
  static const char* blink;
  static const char* reverse;
  static const char* enditalic;
  static const char* endunderline;
  static const char* endblink;
  static const char* endreverse;
}; // struct StreamModifiers


/**
 * @brief Chooses a color from a 256 color map for a foreground color.
 *
 * @param i The color number between 0 and 255.
 * @returns A string describing a color code.
 */
std::string color(size_t i);

std::map<std::string, std::string>& color_map();

std::string color(const std::string id);

/**
 * @brief Chooses a color from a 256 color map for a background color.
 *
 * @param i The color number between 0 and 255.
 * @returns A string describing a color code.
 */
std::string backcolor(size_t i);

// maybe you want to choose your own color
size_t template_color_chooser(size_t i);

bool terminal_supports_color();

/**
 * @brief Highlights templates depending on the "template"-level.
 *
 * @param str The string containing the template string
 * @param maxlevel The maximal template-level the string is reduced to.
 * @returns A colored template string.
 */
std::string highlight_template(std::string str, size_t maxlevel = 10000);

/**
 * @brief A simple function highlighting a whole string in a specified foreground color.
 *
 * @param str The string you want to highlight.
 * @param colornr A color number from a 256 color map between 0 and 255.
 * @returns The highlighted string.
 */
std::string highlight_string(std::string str, size_t colornr = 0);

std::string color_string(const std::string str, const std::string clr = Colors::brown);

std::string color_string_red(const std::string str);

/**
 * @brief Highlights a substring of another string in a specified color.
 *
 * @param str The string where you want to highlight substrings.
 * @param substr The sub string you want to highlight in str.
 * @param colornr A color number from a 256 color map between 0 and 255.
 * @returns The highlighted string.
 */
std::string highlight_search_string(std::string str, std::string substr, size_t colornr = 0);

} // namespace Common
} // namespace XT
} // namespace Dune

#endif // end of DUNE_XT_COMMON_COLOR_HH
