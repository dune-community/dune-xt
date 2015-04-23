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


namespace internal {


static inline std::string trim_copy_safely(const std::string& str_in)
{
  const std::string str_out = boost::algorithm::trim_copy(str_in);
  if (str_out.find(";") != std::string::npos)
    DUNE_THROW(Exceptions::conversion_error,
               "There was an error while parsing the string below. "
                   << "The value contained a ';': '"
                   << str_out
                   << "'!\n"
                   << "This usually happens if you try to get a matrix expression with a vector type "
                   << "or if you are missing the white space after the ';' in a matrix expression!\n");
  return str_out;
} // ... trim_copy_safely(...)


template <class T>
static inline T convert_safely(const std::string ss)
{
  try {
    return boost::lexical_cast<T, std::string>(ss);
  } catch (boost::bad_lexical_cast& e) {
    DUNE_THROW(Exceptions::conversion_error,
               "Error in boost while converting the string '" << ss << "' to type '" << Typename<T>::value() << "':\n"
                                                              << e.what());
  } catch (std::exception& e) {
    DUNE_THROW(Exceptions::conversion_error,
               "Error in the stl while converting the string '" << ss << "' to type '" << Typename<T>::value() << "':\n"
                                                                << e.what());
  }
} // ... convert_safely(...)


// unspecialized variant
template <class T, bool anything = true>
struct Helper
{
  static inline T from_string(const std::string ss)
  {
    return convert_safely<T>(ss);
  }
}; // struct Helper


// variant for bool, to correctly parse true and false
template <bool anything>
struct Helper<bool, anything>
{
  static inline bool from_string(const std::string ss)
  {
    std::string ss_lower_case = ss;
    std::transform(ss_lower_case.begin(), ss_lower_case.end(), ss_lower_case.begin(), ::tolower);
    if (ss_lower_case == "true")
      return true;
    else if (ss_lower_case == "false")
      return false;
    else
      return convert_safely<bool>(ss);
  }
}; // struct Helper< bool, ... >


// variant for all basic types supported by std::sto*
#define DUNE_STUFF_COMMON_STRING_GENERATE_HELPER(tn, tns)                                                              \
  template <bool anything>                                                                                             \
  struct Helper<tn, anything>                                                                                          \
  {                                                                                                                    \
    static inline tn from_string(const std::string ss)                                                                 \
    {                                                                                                                  \
      return std::sto##tns(ss);                                                                                        \
    }                                                                                                                  \
  };

DUNE_STUFF_COMMON_STRING_GENERATE_HELPER(int, i)
DUNE_STUFF_COMMON_STRING_GENERATE_HELPER(long, l)
DUNE_STUFF_COMMON_STRING_GENERATE_HELPER(long long, ll)
DUNE_STUFF_COMMON_STRING_GENERATE_HELPER(unsigned long, ul)
DUNE_STUFF_COMMON_STRING_GENERATE_HELPER(unsigned long long, ull)
DUNE_STUFF_COMMON_STRING_GENERATE_HELPER(float, f)
DUNE_STUFF_COMMON_STRING_GENERATE_HELPER(double, d)
DUNE_STUFF_COMMON_STRING_GENERATE_HELPER(long double, ld)

#undef DUNE_STUFF_COMMON_STRING_GENERATE_HELPER


// variant for everything that is not a matrix or a vector
template <class T>
static inline typename std::enable_if<!is_vector<T>::value && !is_matrix<T>::value, T>::type
from_string(const std::string ss, const size_t UNUSED_UNLESS_DEBUG(rows) = 0,
            const size_t UNUSED_UNLESS_DEBUG(cols) = 0)
{
  assert(rows == 0);
  assert(cols == 0);
  return Helper<T>::from_string(ss);
}

template <class VectorType>
static inline typename std::enable_if<is_vector<VectorType>::value, VectorType>::type
from_string(const std::string& ss, const size_t size, const size_t UNUSED_UNLESS_DEBUG(cols) = 0)
{
  auto vector_str = ss;
  typedef typename VectorAbstraction<VectorType>::S S;
  assert(cols == 0);
  // check if this is a vector
  if (vector_str.substr(0, 1) == "[" && vector_str.substr(vector_str.size() - 1, 1) == "]") {
    vector_str = vector_str.substr(1, vector_str.size() - 2);
    // we treat this as a vector and split along ' '
    const auto tokens = tokenize<std::string>(vector_str, " ", boost::algorithm::token_compress_on);
    if (size > 0 && tokens.size() < size)
      DUNE_THROW(Exceptions::conversion_error,
                 "Vector expression (see below) has only " << tokens.size() << " elements but " << size
                                                           << " elements were requested!"
                                                           << "\n"
                                                           << "'["
                                                           << vector_str
                                                           << "]'");
    const size_t automatic_size = (size > 0) ? std::min(tokens.size(), size) : tokens.size();
    const size_t actual_size =
        VectorAbstraction<VectorType>::has_static_size ? VectorAbstraction<VectorType>::static_size : automatic_size;
    if (actual_size > automatic_size)
      DUNE_THROW(Exceptions::conversion_error,
                 "Vector expression (see below) has only " << automatic_size << " elements but " << actual_size
                                                           << " elements are required for this VectorType ("
                                                           << Typename<VectorType>::value()
                                                           << ")!"
                                                           << "\n"
                                                           << "'["
                                                           << vector_str
                                                           << "]'");
    VectorType ret = VectorAbstraction<VectorType>::create(actual_size);
    for (size_t ii = 0; ii < actual_size; ++ii)
      ret[ii] = from_string<S>(trim_copy_safely(tokens[ii]));
    return ret;
  } else {
    // we treat this as a scalar
    const auto val              = from_string<S>(trim_copy_safely(vector_str));
    const size_t automatic_size = (size == 0 ? 1 : size);
    const size_t actual_size =
        VectorAbstraction<VectorType>::has_static_size ? VectorAbstraction<VectorType>::static_size : automatic_size;
    if (actual_size > automatic_size && automatic_size != 1)
      DUNE_THROW(Exceptions::conversion_error,
                 "Vector expression (see below) has only " << automatic_size << " elements but " << actual_size
                                                           << " elements are required for this VectorType ("
                                                           << Typename<VectorType>::value()
                                                           << ")!"
                                                           << "\n"
                                                           << "'["
                                                           << vector_str
                                                           << "]'");
    VectorType ret = VectorAbstraction<VectorType>::create(actual_size);
    for (size_t ii = 0; ii < std::min(actual_size, ret.size()); ++ii)
      ret[ii] = val;
    return ret;
  }
} // ... from_string(...)

template <class MatrixType>
static inline typename std::enable_if<is_matrix<MatrixType>::value, MatrixType>::type
from_string(std::string matrix_str, const size_t rows, const size_t cols)
{
  typedef typename MatrixAbstraction<MatrixType>::S S;
  // check if this is a matrix
  if (matrix_str.substr(0, 1) == "[" && matrix_str.substr(matrix_str.size() - 1, 1) == "]") {
    matrix_str = matrix_str.substr(1, matrix_str.size() - 2);
    // we treat this as a matrix and split along ';' to obtain the rows
    const auto row_tokens = tokenize<std::string>(matrix_str, ";", boost::algorithm::token_compress_on);
    if (rows > 0 && row_tokens.size() < rows)
      DUNE_THROW(Exceptions::conversion_error,
                 "Matrix expression (see below) has only " << row_tokens.size() << " rows but " << rows
                                                           << " rows were requested!"
                                                           << "\n"
                                                           << "'["
                                                           << matrix_str
                                                           << "]'");
    const size_t automatic_rows = (rows > 0) ? std::min(row_tokens.size(), rows) : row_tokens.size();
    const size_t actual_rows =
        MatrixAbstraction<MatrixType>::has_static_size ? MatrixAbstraction<MatrixType>::static_rows : automatic_rows;
    if (actual_rows > automatic_rows)
      DUNE_THROW(Exceptions::conversion_error,
                 "Matrix expression (see below) has only " << automatic_rows << " rows but " << actual_rows
                                                           << " rows are required for this MatrixType ("
                                                           << Typename<MatrixType>::value()
                                                           << ")!"
                                                           << "\n"
                                                           << "'["
                                                           << matrix_str
                                                           << "]'");
    // compute the number of columns the matrix will have
    size_t min_cols = std::numeric_limits<size_t>::max();
    for (size_t rr = 0; rr < actual_rows; ++rr) {
      const auto row_token = boost::algorithm::trim_copy(row_tokens[rr]);
      // we treat this as a vector, so we split along ' '
      const auto column_tokens = tokenize<std::string>(row_token, " ", boost::algorithm::token_compress_on);
      min_cols                 = std::min(min_cols, column_tokens.size());
    }
    if (cols > 0 && min_cols < cols)
      DUNE_THROW(Exceptions::conversion_error,
                 "Matrix expression (see below) has only " << min_cols << " columns but " << cols
                                                           << " columns were requested!"
                                                           << "\n"
                                                           << "'["
                                                           << matrix_str
                                                           << "]'");
    const auto automatic_cols = (cols > 0) ? std::min(min_cols, cols) : min_cols;
    const size_t actual_cols =
        MatrixAbstraction<MatrixType>::has_static_size ? MatrixAbstraction<MatrixType>::static_cols : automatic_cols;
    if (actual_cols > automatic_cols)
      DUNE_THROW(Exceptions::conversion_error,
                 "Matrix expression (see below) has only " << automatic_cols << " cols but " << actual_cols
                                                           << " cols are required for this MatrixType ("
                                                           << Typename<MatrixType>::value()
                                                           << ")!"
                                                           << "\n"
                                                           << "'["
                                                           << matrix_str
                                                           << "]'");
    MatrixType ret = MatrixAbstraction<MatrixType>::create(actual_rows, actual_cols);
    // now we do the same again and build the actual matrix
    for (size_t rr = 0; rr < actual_rows; ++rr) {
      const std::string row_token = boost::algorithm::trim_copy(row_tokens[rr]);
      const auto column_tokens = tokenize<std::string>(row_token, " ", boost::algorithm::token_compress_on);
      for (size_t cc = 0; cc < actual_cols; ++cc)
        MatrixAbstraction<MatrixType>::set_entry(ret, rr, cc, from_string<S>(trim_copy_safely(column_tokens[cc])));
    }
    return ret;
  } else {
    // we treat this as a scalar
    const S val                 = from_string<S>(trim_copy_safely(matrix_str));
    const size_t automatic_rows = (rows == 0 ? 1 : rows);
    const size_t actual_rows =
        MatrixAbstraction<MatrixType>::has_static_size ? MatrixAbstraction<MatrixType>::static_rows : automatic_rows;
    if (actual_rows > automatic_rows)
      DUNE_THROW(Exceptions::conversion_error,
                 "Matrix expression (see below) has only " << automatic_rows << " rows but " << actual_rows
                                                           << " rows are required for this MatrixType ("
                                                           << Typename<MatrixType>::value()
                                                           << ")!"
                                                           << "\n"
                                                           << "'["
                                                           << matrix_str
                                                           << "]'");
    const size_t automatic_cols = (cols == 0 ? 1 : cols);
    const size_t actual_cols =
        MatrixAbstraction<MatrixType>::has_static_size ? MatrixAbstraction<MatrixType>::static_cols : automatic_cols;
    if (actual_cols > automatic_cols)
      DUNE_THROW(Exceptions::conversion_error,
                 "Matrix expression (see below) has only " << automatic_cols << " cols but " << actual_cols
                                                           << " cols are required for this MatrixType ("
                                                           << Typename<MatrixType>::value()
                                                           << ")!"
                                                           << "\n"
                                                           << "'["
                                                           << matrix_str
                                                           << "]'");
    MatrixType ret = MatrixAbstraction<MatrixType>::create(actual_rows, actual_cols);
    for (size_t rr = 0; rr < std::min(actual_rows, MatrixAbstraction<MatrixType>::rows(ret)); ++rr)
      for (size_t cc = 0; cc < std::min(actual_cols, MatrixAbstraction<MatrixType>::cols(ret)); ++cc)
        MatrixAbstraction<MatrixType>::set_entry(ret, rr, cc, val);
    return ret;
  }
} // ... from_string(...)


// variant for everything that is not a matrix, a vector or any of the types specified below
template <class T>
static inline typename std::enable_if<!is_vector<T>::value && !is_matrix<T>::value, std::string>::type
to_string(const T& ss)
{
  return std::to_string(ss);
}

template <int size>
static inline std::string to_string(const Dune::bigunsignedint<size>& ss)
{
  std::stringstream os;
  os << ss;
  return os.str();
}

inline std::string to_string(const char* ss)
{
  return std::string(ss);
}

inline std::string to_string(char ss)
{
  return std::string(1, ss);
}

inline std::string to_string(const std::string ss)
{
  return std::string(ss);
}

template <class V>
static inline typename std::enable_if<is_vector<V>::value, std::string>::type to_string(const V& vec)
{
  std::string ret = "[";
  for (auto ii : valueRange(vec.size())) {
    if (ii > 0)
      ret += " ";
    ret += to_string(vec[ii]);
  }
  ret += "]";
  return ret;
} // ... to_string(...)

template <class M>
static inline typename std::enable_if<is_matrix<M>::value, std::string>::type to_string(const M& mat)
{
  std::string ret = "[";
  for (auto rr : valueRange(MatrixAbstraction<M>::rows(mat))) {
    if (rr > 0)
      ret += "; ";
    for (auto cc : valueRange(MatrixAbstraction<M>::cols(mat))) {
      if (cc > 0)
        ret += " ";
      ret += to_string(MatrixAbstraction<M>::get_entry(mat, rr, cc));
    }
  }
  ret += "]";
  return ret;
} // ... to_string(...)


} // namespace internal


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
static inline T fromString(const std::string ss, const size_t size = 0, const size_t cols = 0)
{
  return internal::from_string<T>(ss, size, cols);
}

/**
 * \brief Converts an object to string.
 * \sa    internal::to_string for implementations
 */
template <class T>
static inline std::string toString(const T& ss)
{
  return internal::to_string(ss);
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
  for (auto DUNE_UNUSED(ii) : valueRange(s.size())) {
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
