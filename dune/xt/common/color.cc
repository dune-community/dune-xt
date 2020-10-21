// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   René Fritze     (2013, 2015 - 2016, 2018 - 2019)
//   Tobias Leibner  (2014, 2018, 2020)

#include <config.h>

#include "color.hh"

namespace Dune::XT::Common {

// foreground colors
const char* Colors::black = "\033[30m";
const char* Colors::red = "\033[31m";
const char* Colors::green = "\033[32m";
const char* Colors::brown = "\033[33m";
const char* Colors::blue = "\033[34m";
const char* Colors::purple = "\033[35m";
const char* Colors::cyan = "\033[36m";
const char* Colors::lightgray = "\033[37m";
// light foreground colors
const char* Colors::darkgray = "\033[1;30m";
const char* Colors::lightred = "\033[1;31m";
const char* Colors::lightgreen = "\033[1;32m";
const char* Colors::yellow = "\033[1;33m";
const char* Colors::lightblue = "\033[1;34m";
const char* Colors::lightpurple = "\033[1;35m";
const char* Colors::lightcyan = "\033[1;36m";
const char* Colors::white = "\033[1;37m";
// background colors
const char* Colors::bblack = "\033[40m";
const char* Colors::bred = "\033[41m";
const char* Colors::bgreen = "\033[42m";
const char* Colors::bbrown = "\033[43m";
const char* Colors::bblue = "\033[44m";
const char* Colors::bpurple = "\033[45m";
const char* Colors::bcyan = "\033[46m";
const char* Colors::blightgray = "\033[47m";
// light background colors
const char* Colors::bdarkgray = "\033[1;40m";
const char* Colors::blightred = "\033[1;41m";
const char* Colors::blightgreen = "\033[1;42m";
const char* Colors::byellow = "\033[1;43m";
const char* Colors::blightblue = "\033[1;44m";
const char* Colors::blightpurple = "\033[1;45m";
const char* Colors::blightcyan = "\033[1;46m";
const char* Colors::bwhite = "\033[1;47m";

// modifiers
const char* StreamModifiers::normal = "\033[0m";
const char* StreamModifiers::bold = "\033[1m";
const char* StreamModifiers::italic = "\033[2m";
const char* StreamModifiers::underline = "\033[4m";
const char* StreamModifiers::blink = "\033[5m";
const char* StreamModifiers::reverse = "\033[7m";
const char* StreamModifiers::enditalic = "\033[22m";
const char* StreamModifiers::endunderline = "\033[24m";
const char* StreamModifiers::endblink = "\033[25m";
const char* StreamModifiers::endreverse = "\033[27m";

std::string color(size_t i)
{
  return "\033[38;5;" + std::to_string(i) + "m";
}

std::map<std::string, std::string>& color_map()
{
  static std::map<std::string, std::string> map_;
  static bool created_ = false;
  if (!created_) {
    map_["black"] = Colors::black;
    map_["blue"] = Colors::blue;
    map_["brown"] = Colors::brown;
    map_["cyan"] = Colors::cyan;
    map_["darkgray"] = Colors::darkgray;
    map_["green"] = Colors::green;
    map_["lightblue"] = Colors::lightblue;
    map_["lightcyan"] = Colors::lightcyan;
    map_["lightgray"] = Colors::lightgray;
    map_["lightgreen"] = Colors::lightgreen;
    map_["lightpurple"] = Colors::lightpurple;
    map_["lightred"] = Colors::lightred;
    map_["purple"] = Colors::purple;
    map_["red"] = Colors::red;
    map_["white"] = Colors::white;
    map_["yellow"] = Colors::yellow;
    created_ = true;
  }
  return map_;
} // ... color_map(...)

std::string color(const std::string id)
{
  const auto search_result = color_map().find(id);
  if (search_result != color_map().end())
    return search_result->second;
  return "";
}

std::string backcolor(size_t i)
{
  return "\033[38;5;" + std::to_string(i) + "m";
}

size_t template_color_chooser(size_t i)
{
  return i % 256;
}

/**
 * This is taken from gtest-all.cc!
 */
bool terminal_supports_color()
{
  const char* const term = getenv("TERM");
  if (term == nullptr)
    return false;
  const auto term_str = std::string(term);
  return term_str == "xterm" || term_str == "xterm-color" || term_str == "xterm-256color" || term_str == "screen"
         || term_str == "linux" || term_str == "cygwin";

} // ... terminal_supports_color(...)

std::string highlight_template(std::string str, size_t maxlevel)
{
  size_t startindex = 0;
  size_t level = 0;
  for (size_t i = 0; i < str.size(); i++) {
    if (str[i] == '<') {
      level++;
      std::string dummy = "\033[38;5;" + std::to_string(template_color_chooser(level)) + "m";
      str.insert(i, dummy);
      i += dummy.size();
      if (level == maxlevel)
        startindex = i + 1;
    } else if (str[i] == '>') {
      level--;
      std::string dummy = "\033[38;5;" + std::to_string(template_color_chooser(level)) + "m";
      str.insert(++i, dummy);
      if (level + 1 == maxlevel) {
        auto size = i - startindex - 1;
        str.erase(startindex, size);
        i = startindex + 1;
      }
      i += dummy.size();
    }
  }
  str += "\033[0m";
  return str;
} // highlight_template

std::string highlight_string(std::string str, size_t colornr)
{
  return "\033[38;5;" + std::to_string(colornr % 256) + "m" + str + "\033[0m"; //"\033[38;5;0m";
}

std::string highlight_search_string(std::string str, std::string substr, size_t colornr)
{
  std::string dummy = "\033[38;5;" + std::to_string(colornr % 256) + "m";
  std::string dummy2 = "\033[0m";
  str.insert(0, dummy2);

  size_t index = str.find(substr, 0);

  while (index != std::string::npos) {
    str.insert(index, dummy);
    str.insert(index + substr.size() + dummy.size(), dummy2);
    index = str.find(substr, index + dummy.size() + substr.size() + dummy2.size());
  }
  return str;
} // highlight_search_string

std::string color_string(const std::string str, const std::string clr)
{
  if (terminal_supports_color())
    return clr + str + StreamModifiers::normal;
  return str;
}

std::string color_string_red(const std::string str)
{
  return color_string(str, Colors::red);
}

} // namespace Dune::XT::Common
