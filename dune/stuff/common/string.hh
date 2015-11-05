// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

/**
 * \file  string.hh
 * \brief contains methods to read/write scalar types and container types from/to strings and methods to manipulate
 *        strings
 **/
#ifndef DUNE_STUFF_COMMON_STRING_HH
#define DUNE_STUFF_COMMON_STRING_HH

#include <cstring>
#include <ctime>
#include <map>
#include <cassert>
#include <algorithm>
#include <ostream>
#include <iomanip>
#include <vector>
#include <string>
#include <ctime>
#include <iostream>

#include <dune/stuff/common/disable_warnings.hh>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/common/array.hh>
#include <dune/common/bigunsignedint.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/densematrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/densevector.hh>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/debug.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/vector.hh>
#include <dune/stuff/la/container/interfaces.hh>

namespace Dune {
namespace Stuff {
namespace Common {

// forward, std::string variant needed in internal::from_string

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
tokenize(const std::string& msg, const std::string& separators,
         const boost::algorithm::token_compress_mode_type mode = boost::algorithm::token_compress_off);

} // namespace Common
} // namespace Stuff
} // namespace Dune

#include "string_internal.hh"

namespace Dune {
namespace Stuff {
namespace Common {

static constexpr const std::size_t default_to_string_precision = 6u;

/**
 * \brief Reads an object from a string.
 * \sa    internal::from_string for implementations
 * \param ss   string to read from
 * \param size Determines the size of the returning container (size if T is a vector type, rows if T is a matrix type,
 *             ignored if T is a scalar type; 0 means automatic).
 * \param cols Determines the number of columns of the returning matrix if T is a matrix type (0 means automatic;
 *             ignored if T is a vector or scalar type).
 */
template <class T>
static inline T fromString(std::string ss, const size_t size = 0, const size_t cols = 0)
{
  return internal::from_string<T>(ss, size, cols);
}

/**
 * \brief Converts an object to string.
 * \sa    internal::to_string for implementations
 */
template <class T>
static inline std::string toString(const T& ss, const size_t precision = default_to_string_precision)
{
  return internal::to_string(ss, precision);
}

/**
  \brief Returns a string of lengths t' whitespace (or whitespace chars).
  \param[in]  t          defines the length of the return string (after conversion to string)
  \param[in]  whitespace char, optional argument, defines entries of return string
  \return     A string of lengths t' whitespace (or whitespace chars).
  **/
template <class T>
std::string whitespaceify(const T& t, const char whitespace = ' ')
{
  const std::string s = toString(t);
  std::string ret = "";
  for (auto DSC_UNUSED(ii) : valueRange(s.size())) {
    ret += whitespace;
  }
  return ret;
} // ... whitespaceify(...)

template <class T>
inline std::vector<T> tokenize(const std::string& msg, const std::string& separators,
                               const boost::algorithm::token_compress_mode_type mode)
{
  std::vector<std::string> strings;
  boost::algorithm::split(strings, msg, boost::algorithm::is_any_of(separators), mode);
  std::vector<T> ret(strings.size());
  size_t i = 0;
  // special case for empty strings to avoid non-default init
  std::generate(
      std::begin(ret), std::end(ret), [&]() { return strings[i++].empty() ? T() : fromString<T>(strings[i - 1]); });
  return ret;
} // ... tokenize(...)

template <>
inline std::vector<std::string> tokenize(const std::string& msg, const std::string& separators,
                                         const boost::algorithm::token_compress_mode_type mode)
{
  std::vector<std::string> strings;
  boost::algorithm::split(strings, msg, boost::algorithm::is_any_of(separators), mode);
  return strings;
}

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

static inline char** vectorToMainArgs(const std::vector<std::string>& args)
{
  char** argv = new char*[args.size()];
  for (auto ii : valueRange(args.size())) {
    argv[ii] = new char[args[ii].length() + 1];
    strcpy(argv[ii], args[ii].c_str());
  }
  return argv;
} // ... vectorToMainArgs(...)

} // namespace Common
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_COMMON_STRING_HH
