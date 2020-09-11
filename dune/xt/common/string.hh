// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2009, 2012 - 2018)
//   René Fritze     (2009 - 2016, 2018 - 2019)
//   Sven Kaulmann   (2011 - 2012)
//   Tobias Leibner  (2014, 2017 - 2020)

/**
 * \file  string.hh
 * \brief contains methods to read/write scalar types and container types from/to strings and methods to manipulate
 *        strings
 **/
#ifndef DUNE_XT_COMMON_STRING_HH
#define DUNE_XT_COMMON_STRING_HH

#include <algorithm>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>

#include <dune/xt/common/disable_warnings.hh>
#include <boost/algorithm/string/constants.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <dune/xt/common/reenable_warnings.hh>

namespace Dune {
namespace XT {
namespace Common {


// forward, std::string variant needed in internal::convert_from_string

/**
 * \brief  convenience wrapper around boost::algorithm::split to split one string into a vector of strings
 * \param  msg        the spring to be split
 * \param  separators a list of separators, duh
 * \param  mode       token_compress_off --> potentially empty strings in return,
 *                    token_compress_on --> empty tokens are discarded
 * \return all tokens in a vector, if msg contains no separators, this'll contain msg as its only element
 **/
template <class T = std::string>
inline std::vector<T>
tokenize(const std::string& msg,
         const std::string& separators,
         const boost::algorithm::token_compress_mode_type mode = boost::algorithm::token_compress_off);


} // namespace Common
} // namespace XT
} // namespace Dune

#include "string_internal.hh"

namespace Dune {
namespace XT {
namespace Common {


static constexpr const std::size_t default_to_string_precision = 6u;


/**
 * \brief Reads an object from a string.
 * \sa    internal::convert_from_string for implementations
 * \param ss   string to read from
 * \param size Determines the size of the returning container (size if T is a vector type, rows if T is a matrix type,
 *             ignored if T is a scalar type; 0 means automatic).
 * \param cols Determines the number of columns of the returning matrix if T is a matrix type (0 means automatic;
 *             ignored if T is a vector or scalar type).
 */
template <class T>
static inline T from_string(std::string ss, const size_t size = 0, const size_t cols = 0)
{
  return internal::convert_from_string<T>(ss, size, cols);
}


/**
 * \brief Converts an object to string.
 * \sa    internal::convert_to_string for implementations
 */
template <class T>
static inline std::string to_string(const T& ss, const size_t precision = default_to_string_precision)
{
  return internal::convert_to_string(ss, precision);
}


/**
 * \brief Converts each character of a string to lower case using std::tolower.
 * \note  This might not do what you expect, given your locale.
 * \sa    http://en.cppreference.com/w/cpp/string/byte/tolower
 * /sa    ::tolower
 */
std::string to_lower(const std::string& ss);


/**
 * \brief Converts each character of a string to upper case using std::toupper.
 * \note  This might not do what you expect, given your locale.
 * \sa    http://en.cppreference.com/w/cpp/string/byte/tolower
 * /sa    ::toupper
 */
std::string to_upper(const std::string& ss);


/**
 * \brief Converts this_sample_string to ThisSampleString.
 */
std::string to_camel_case(const std::string& ss);


/**
  \brief Returns a string of lengths t' whitespace (or whitespace chars).
  \param[in]  t          defines the length of the return string (after conversion to string)
  \param[in]  whitespace char, optional argument, defines entries of return string
  \return     A string of lengths t' whitespace (or whitespace chars).
  **/
template <class T>
std::string whitespaceify(const T& t, const char whitespace = ' ')
{
  const std::string s = to_string(t);
  std::string ret = "";
  for ([[maybe_unused]] auto ii : value_range(s.size()))
    ret += whitespace;
  return ret;
} // ... whitespaceify(...)


template <class T>
inline std::vector<T>
tokenize(const std::string& msg, const std::string& separators, const boost::algorithm::token_compress_mode_type mode)
{
  std::vector<std::string> strings;
  boost::algorithm::split(strings, msg, boost::algorithm::is_any_of(separators), mode);
  std::vector<T> ret(strings.size());
  size_t i = 0;
  // special case for empty strings to avoid non-default init
  std::generate(
      std::begin(ret), std::end(ret), [&]() { return strings[i++].empty() ? T() : from_string<T>(strings[i - 1]); });
  return ret;
} // ... tokenize(...)


template <>
inline std::vector<std::string>
tokenize(const std::string& msg, const std::string& separators, const boost::algorithm::token_compress_mode_type mode)
{
  std::vector<std::string> strings;
  boost::algorithm::split(strings, msg, boost::algorithm::is_any_of(separators), mode);
  return strings;
}

/**
  \brief Removes whitespace from front and back of the string.
  \param[in]  s String to be trimmed
**/
void trim(std::string& v);

/**
  \brief Removes whitespace from front and back of each string in the vector.
  \param[in]  v Vector of strings to be trimmed
**/
void trim(std::vector<std::string>& v);

//! returns string with local time in current locale's format
inline std::string stringFromTime(time_t cur_time = time(NULL))
{
  return ctime(&cur_time);
}


//! helper struct for lexical cast
// see http://stackoverflow.com/a/2079728
template <typename ElemT>
struct HexToString
{
  ElemT value;
  operator ElemT() const
  {
    return value;
  }

  friend std::istream& operator>>(std::istream& in, HexToString& out)
  {
    in >> std::hex >> out.value;
    return in;
  }
};


static inline char** vector_to_main_args(const std::vector<std::string>& args)
{
  char** argv = new char*[args.size()];
  for (auto ii : value_range(args.size())) {
    argv[ii] = new char[args[ii].length() + 1];
    strcpy(argv[ii], args[ii].c_str());
  }
  return argv;
} // ... vector_to_main_args(...)


} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_COMMON_STRING_HH
