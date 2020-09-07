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
//   Tobias Leibner  (2014 - 2015, 2017 - 2020)

#ifndef DUNE_XT_COMMON_STRING_INTERNAL_HH
#define DUNE_XT_COMMON_STRING_INTERNAL_HH

#include <cstring>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <string>

#include <dune/xt/common/disable_warnings.hh>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/algorithm/string/constants.hpp>
#include <dune/xt/common/reenable_warnings.hh>

#include <dune/xt/common/debug.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/vector.hh>
#include <dune/xt/common/matrix.hh>

namespace Dune {
namespace XT {
namespace Common {

#ifndef DUNE_XT_COMMON_STRING_HH
// only necessary for headercheck
template <class T = std::string>
inline std::vector<T>
tokenize(const std::string& msg,
         const std::string& separators,
         const boost::algorithm::token_compress_mode_type mode = boost::algorithm::token_compress_off);
#endif // DUNE_XT_COMMON_STRING_HH

namespace internal {

static inline std::string trim_copy_safely(std::string str_in)
{
  const std::string str_out = boost::algorithm::trim_copy(str_in);
  if (str_out.find(";") != std::string::npos)
    DUNE_THROW(Exceptions::conversion_error,
               "There was an error while parsing the string below. "
                   << "The value contained a ';': '" << str_out << "'!\n"
                   << "This usually happens if you try to get a matrix expression with a vector type "
                   << "or if you are missing the white space after the ';' in a matrix expression!\n");
  return str_out;
} // ... trim_copy_safely(...)

template <class T>
static inline T convert_safely(std::string ss)
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
  return T();
} // ... convert_safely(...)

template <class T>
T convert_from_string(std::string ss, const size_t rows = 0, const size_t cols = 0);

// unspecialized variant
template <class T>
struct Helper
{
  static inline T convert_from_string(std::string ss)
  {
    return convert_safely<T>(ss);
  }
}; // struct Helper

// variant for bool, to correctly parse true and false
template <>
struct Helper<bool>
{
  static inline bool convert_from_string(std::string ss)
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
#define DUNE_XT_COMMON_STRING_GENERATE_HELPER(tn, tns)                                                                 \
  template <>                                                                                                          \
  struct Helper<tn>                                                                                                    \
  {                                                                                                                    \
    static inline tn convert_from_string(std::string ss)                                                               \
    {                                                                                                                  \
      try {                                                                                                            \
        return std::sto##tns(ss);                                                                                      \
      } catch (const std::exception& ee) {                                                                             \
        DUNE_THROW(Exceptions::conversion_error,                                                                       \
                   "in stl when converting '" << ss << "' to '" << Typename<tn>::value() << "': " << ee.what());       \
      }                                                                                                                \
    }                                                                                                                  \
  };

DUNE_XT_COMMON_STRING_GENERATE_HELPER(int, i)
DUNE_XT_COMMON_STRING_GENERATE_HELPER(long, l)
DUNE_XT_COMMON_STRING_GENERATE_HELPER(long long, ll)
DUNE_XT_COMMON_STRING_GENERATE_HELPER(unsigned long, ul)
DUNE_XT_COMMON_STRING_GENERATE_HELPER(unsigned long long, ull)
DUNE_XT_COMMON_STRING_GENERATE_HELPER(float, f)
DUNE_XT_COMMON_STRING_GENERATE_HELPER(double, d)
DUNE_XT_COMMON_STRING_GENERATE_HELPER(long double, ld)

#undef DUNE_XT_COMMON_STRING_GENERATE_HELPER

// variant for unsigned int as there is no stoui
template <>
struct Helper<unsigned int>
{
  static inline unsigned int convert_from_string(std::string ss)
  {
    try {
      return XT::Common::numeric_cast<unsigned int>(std::stoul(ss));
    } catch (const XT::Common::Exceptions::external_error& ee) {
      DUNE_THROW(Exceptions::conversion_error,
                 "in XT::Common::numeric_cast when converting '" << ss << "' to '" << Typename<unsigned int>::value()
                                                                 << "': " << ee.what());
    } catch (const std::exception& ee) {
      DUNE_THROW(Exceptions::conversion_error,
                 "in stl when converting '" << ss << "' to '" << Typename<unsigned long>::value()
                                            << "': " << ee.what());
    }
  }
};

template <class ComplexType>
ComplexType complex_from_string(std::string ss, const size_t /*size*/ = 0, const size_t /*cols*/ = 0)
{
  boost::algorithm::trim(ss);
  if (ss.size() < 1)
    DUNE_THROW(Exceptions::conversion_error, "Error converting " << ss << " (too short)");
  using namespace std;
  typedef typename ComplexType::value_type T;
  T re(0), im(0);
  const auto sign_pos = ss.find("+", 1) != string::npos ? ss.find("+", 1) : ss.find("-", 1);
  auto im_pos = ss.find("i");
  if (sign_pos == string::npos && im_pos != string::npos) {
    // pure virtual case
    im = convert_from_string<T>(ss.substr(0, im_pos));
  } else {
    re = convert_from_string<T>(ss.substr(0, sign_pos));
  }
  if (sign_pos != string::npos) {
    ss = ss.substr(sign_pos);
    im_pos = ss.find("i");
    if (im_pos == string::npos)
      DUNE_THROW(Exceptions::conversion_error, "Error converting " << ss << " no imaginary unit");
    im = convert_from_string<T>(ss.substr(0, im_pos));
  }
  return ComplexType(re, im);
}

template <class VectorType>
VectorType vector_from_string(std::string vector_str, const size_t size, const size_t DXTC_DEBUG_ONLY(cols) = 0)
{
  typedef typename VectorAbstraction<VectorType>::S S;
  DXT_ASSERT(cols == 0);
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
                                                           << "'[" << vector_str << "]'");
    const size_t automatic_size = (size > 0) ? std::min(tokens.size(), size) : tokens.size();
    const size_t actual_size =
        VectorAbstraction<VectorType>::has_static_size ? VectorAbstraction<VectorType>::static_size : automatic_size;
    if (actual_size > automatic_size)
      DUNE_THROW(Exceptions::conversion_error,
                 "Vector expression (see below) has only " << automatic_size << " elements but " << actual_size
                                                           << " elements are required for this VectorType ("
                                                           << Typename<VectorType>::value() << ")!"
                                                           << "\n"
                                                           << "'[" << vector_str << "]'");
    VectorType ret = VectorAbstraction<VectorType>::create(actual_size);
    for (size_t ii = 0; ii < actual_size; ++ii)
      ret[ii] = convert_from_string<S>(trim_copy_safely(tokens[ii]));
    return ret;
  } else {
    // we treat this as a scalar
    const auto val = convert_from_string<S>(trim_copy_safely(vector_str));
    const size_t automatic_size = (size == 0 ? 1 : size);
    const size_t actual_size =
        VectorAbstraction<VectorType>::has_static_size ? VectorAbstraction<VectorType>::static_size : automatic_size;
    if (actual_size > automatic_size && automatic_size != 1)
      DUNE_THROW(Exceptions::conversion_error,
                 "Vector expression (see below) has only " << automatic_size << " elements but " << actual_size
                                                           << " elements are required for this VectorType ("
                                                           << Typename<VectorType>::value() << ")!"
                                                           << "\n"
                                                           << "'[" << vector_str << "]'");
    VectorType ret = VectorAbstraction<VectorType>::create(actual_size);
    for (size_t ii = 0; ii < std::min(actual_size, ret.size()); ++ii)
      ret[ii] = val;
    return ret;
  }
} // ... convert_from_string(...)

template <class MatrixType>
MatrixType matrix_from_string(std::string matrix_str, const size_t rows, const size_t cols)
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
                                                           << "'[" << matrix_str << "]'");
    const size_t automatic_rows = (rows > 0) ? std::min(row_tokens.size(), rows) : row_tokens.size();
    const size_t actual_rows =
        MatrixAbstraction<MatrixType>::has_static_size ? MatrixAbstraction<MatrixType>::static_rows : automatic_rows;
    if (actual_rows > automatic_rows)
      DUNE_THROW(Exceptions::conversion_error,
                 "Matrix expression (see below) has only " << automatic_rows << " rows but " << actual_rows
                                                           << " rows are required for this MatrixType ("
                                                           << Typename<MatrixType>::value() << ")!"
                                                           << "\n"
                                                           << "'[" << matrix_str << "]'");
    // compute the number of columns the matrix will have
    size_t min_cols = std::numeric_limits<size_t>::max();
    for (size_t rr = 0; rr < actual_rows; ++rr) {
      const auto row_token = boost::algorithm::trim_copy(row_tokens[rr]);
      // we treat this as a vector, so we split along ' '
      const auto column_tokens = tokenize<std::string>(row_token, " ", boost::algorithm::token_compress_on);
      min_cols = std::min(min_cols, column_tokens.size());
    }
    if (cols > 0 && min_cols < cols)
      DUNE_THROW(Exceptions::conversion_error,
                 "Matrix expression (see below) has only " << min_cols << " columns but " << cols
                                                           << " columns were requested!"
                                                           << "\n"
                                                           << "'[" << matrix_str << "]'");
    const auto automatic_cols = (cols > 0) ? std::min(min_cols, cols) : min_cols;
    const size_t actual_cols =
        MatrixAbstraction<MatrixType>::has_static_size ? MatrixAbstraction<MatrixType>::static_cols : automatic_cols;
    if (actual_cols > automatic_cols)
      DUNE_THROW(Exceptions::conversion_error,
                 "Matrix expression (see below) has only " << automatic_cols << " cols but " << actual_cols
                                                           << " cols are required for this MatrixType ("
                                                           << Typename<MatrixType>::value() << ")!"
                                                           << "\n"
                                                           << "'[" << matrix_str << "]'");
    MatrixType ret = MatrixAbstraction<MatrixType>::create(actual_rows, actual_cols);
    // now we do the same again and build the actual matrix
    for (size_t rr = 0; rr < actual_rows; ++rr) {
      const std::string row_token = boost::algorithm::trim_copy(row_tokens[rr]);
      const auto column_tokens = tokenize<std::string>(row_token, " ", boost::algorithm::token_compress_on);
      for (size_t cc = 0; cc < actual_cols; ++cc)
        MatrixAbstraction<MatrixType>::set_entry(
            ret, rr, cc, convert_from_string<S>(trim_copy_safely(column_tokens[cc])));
    }
    return ret;
  } else {
    // we treat this as a scalar
    const S val = convert_from_string<S>(trim_copy_safely(matrix_str));
    const size_t automatic_rows = (rows == 0 ? 1 : rows);
    const size_t actual_rows =
        MatrixAbstraction<MatrixType>::has_static_size ? MatrixAbstraction<MatrixType>::static_rows : automatic_rows;
    if (actual_rows > automatic_rows)
      DUNE_THROW(Exceptions::conversion_error,
                 "Matrix expression (see below) has only " << automatic_rows << " rows but " << actual_rows
                                                           << " rows are required for this MatrixType ("
                                                           << Typename<MatrixType>::value() << ")!"
                                                           << "\n"
                                                           << "'[" << matrix_str << "]'");
    const size_t automatic_cols = (cols == 0 ? 1 : cols);
    const size_t actual_cols =
        MatrixAbstraction<MatrixType>::has_static_size ? MatrixAbstraction<MatrixType>::static_cols : automatic_cols;
    if (actual_cols > automatic_cols)
      DUNE_THROW(Exceptions::conversion_error,
                 "Matrix expression (see below) has only " << automatic_cols << " cols but " << actual_cols
                                                           << " cols are required for this MatrixType ("
                                                           << Typename<MatrixType>::value() << ")!"
                                                           << "\n"
                                                           << "'[" << matrix_str << "]'");
    MatrixType ret = MatrixAbstraction<MatrixType>::create(actual_rows, actual_cols);
    for (size_t rr = 0; rr < std::min(actual_rows, MatrixAbstraction<MatrixType>::rows(ret)); ++rr)
      for (size_t cc = 0; cc < std::min(actual_cols, MatrixAbstraction<MatrixType>::cols(ret)); ++cc)
        MatrixAbstraction<MatrixType>::set_entry(ret, rr, cc, val);
    return ret;
  }
} // ... convert_from_string(...)

// main function that dispatches to specializations
template <class T>
T convert_from_string(std::string ss, const size_t rows, const size_t cols)
{
  if constexpr (is_complex<T>::value) {
    return complex_from_string<T>(ss);
  } else if constexpr (is_vector<T>::value) {
    return vector_from_string<T>(ss, rows, cols);
  } else if constexpr (is_matrix<T>::value) {
    return matrix_from_string<T>(ss, rows, cols);
  } else {
    DXT_ASSERT(rows == 0);
    DXT_ASSERT(cols == 0);
    return Helper<T>::convert_from_string(ss);
  }
}

// variant for everything that is not a matrix, a vector or any of the types specified below
template <class T>
static inline typename std::enable_if<!is_vector<T>::value && !is_matrix<T>::value, std::string>::type
convert_to_string(const T& ss, const std::size_t precision)
{
  std::ostringstream out;
  out << std::setprecision(boost::numeric_cast<int>(precision)) << ss;
  return out.str();
}

template <typename T>
static inline std::string convert_to_string(const std::complex<T>& val, const std::size_t precision)
{
  using namespace std;
  stringstream os;
  const auto im = imag(val);
  const string sign = signum(im) < 0 ? "" : "+";
  const auto real_str = convert_to_string(real(val), precision);
  const auto imag_str = convert_to_string(im, precision);
  os << real_str << sign << imag_str << "i";
  return os.str();
}

template <int size>
static inline std::string convert_to_string(const Dune::bigunsignedint<size>& ss, const std::size_t /*precision*/)
{
  std::stringstream os;
  os << ss;
  return os.str();
}

inline std::string convert_to_string(const char* ss, const std::size_t /*precision*/)
{
  return std::string(ss);
}

inline std::string convert_to_string(char ss, const std::size_t /*precision*/)
{
  return std::string(1, ss);
}

inline std::string convert_to_string(const std::string ss, const std::size_t /*precision*/)
{
  return std::string(ss);
}

// forward such that vectors of matrices can be converted
template <class M>
static inline typename std::enable_if<is_matrix<M>::value, std::string>::type
convert_to_string(const M& /*mat*/, const std::size_t /*precision*/);

template <class V>
static inline typename std::enable_if<is_vector<V>::value, std::string>::type
convert_to_string(const V& vec, const std::size_t precision)
{
  std::string ret = "[";
  for (auto ii : value_range(vec.size())) {
    if (ii > 0)
      ret += " ";
    ret += convert_to_string(vec[ii], precision);
  }
  ret += "]";
  return ret;
} // ... convert_to_string(...)

template <class M>
static inline typename std::enable_if<is_matrix<M>::value, std::string>::type
convert_to_string(const M& mat, const std::size_t precision)
{
  std::string ret = "[";
  for (auto rr : value_range(MatrixAbstraction<M>::rows(mat))) {
    if (rr > 0)
      ret += "; ";
    for (auto cc : value_range(MatrixAbstraction<M>::cols(mat))) {
      if (cc > 0)
        ret += " ";
      ret += convert_to_string(MatrixAbstraction<M>::get_entry(mat, rr, cc), precision);
    }
  }
  ret += "]";
  return ret;
} // ... convert_to_string(...)


} // namespace internal
} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_COMMON_STRING_INTERNAL_HH
