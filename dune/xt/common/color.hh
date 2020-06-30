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

#include <sstream>
#include <iostream>
#include <map>
#include <string>

#ifdef __GNUC__
#  include <cxxabi.h>
#endif
#include <string>

namespace Dune {
namespace XT {
namespace Common {

#define DS_CONST_CHAR const char*

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
  static DS_CONST_CHAR black;
  static DS_CONST_CHAR red;
  static DS_CONST_CHAR green;
  static DS_CONST_CHAR brown;
  static DS_CONST_CHAR blue;
  static DS_CONST_CHAR purple;
  static DS_CONST_CHAR cyan;
  static DS_CONST_CHAR lightgray;
  // light foreground colors
  static DS_CONST_CHAR darkgray;
  static DS_CONST_CHAR lightred;
  static DS_CONST_CHAR lightgreen;
  static DS_CONST_CHAR yellow;
  static DS_CONST_CHAR lightblue;
  static DS_CONST_CHAR lightpurple;
  static DS_CONST_CHAR lightcyan;
  static DS_CONST_CHAR white;

  // background colors
  static DS_CONST_CHAR bblack;
  static DS_CONST_CHAR bred;
  static DS_CONST_CHAR bgreen;
  static DS_CONST_CHAR bbrown;
  static DS_CONST_CHAR bblue;
  static DS_CONST_CHAR bpurple;
  static DS_CONST_CHAR bcyan;
  static DS_CONST_CHAR blightgray;
  // light background colors
  static DS_CONST_CHAR bdarkgray;
  static DS_CONST_CHAR blightred;
  static DS_CONST_CHAR blightgreen;
  static DS_CONST_CHAR byellow;
  static DS_CONST_CHAR blightblue;
  static DS_CONST_CHAR blightpurple;
  static DS_CONST_CHAR blightcyan;
  static DS_CONST_CHAR bwhite;
}; // struct Colors

// modifiers
struct StreamModifiers
{
  static DS_CONST_CHAR normal;
  static DS_CONST_CHAR bold;
  static DS_CONST_CHAR italic;
  static DS_CONST_CHAR underline;
  static DS_CONST_CHAR blink;
  static DS_CONST_CHAR reverse;
  static DS_CONST_CHAR enditalic;
  static DS_CONST_CHAR endunderline;
  static DS_CONST_CHAR endblink;
  static DS_CONST_CHAR endreverse;
}; // struct StreamModifiers

#undef DS_CONST_CHAR

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
