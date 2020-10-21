// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2014, 2016 - 2017)
//   René Fritze     (2013, 2015 - 2016, 2018 - 2019)
//   Tobias Leibner  (2018 - 2020)

#include <config.h>

#include <algorithm>
#include <cctype>
#include <string>
#include <sstream>

#include "string.hh"

namespace Dune::XT::Common {


std::string to_lower(const std::string& ss)
{
  std::string out(ss);
  std::transform(out.begin(), out.end(), out.begin(), ::tolower);
  return out;
}


std::string to_upper(const std::string& ss)
{
  std::string out(ss);
  std::transform(out.begin(), out.end(), out.begin(), ::toupper);
  return out;
}


std::string to_camel_case(const std::string& ss)
{
  std::stringstream out;
  for (auto&& word : tokenize(ss, "_", boost::algorithm::token_compress_on)) {
    if (word.size() > 0) {
      out << to_upper(word.substr(0, 1));
      out << word.substr(1);
    }
  }
  return out.str();
}

void trim(std::string& s)
{
  auto wsfront = std::find_if_not(s.begin(), s.end(), [](int c) { return std::isspace(c); });
  auto wsback = std::find_if_not(s.rbegin(), s.rend(), [](int c) { return std::isspace(c); }).base();
  s = (wsback <= wsfront ? std::string() : std::string(wsfront, wsback));
}

void trim(std::vector<std::string>& v)
{
  for (auto& s : v)
    trim(s);
}


} // namespace Dune::XT::Common
