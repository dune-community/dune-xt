// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2010, 2012 - 2017)
//   René Fritze     (2009 - 2016, 2018 - 2019)
//   Sven Kaulmann   (2011 - 2012)
//   Tobias Leibner  (2014, 2018, 2020)

#ifndef DUNE_XT_COMMON_PRINT_HH
#define DUNE_XT_COMMON_PRINT_HH

#include <string>
#include <iostream>

#include <dune/xt/common/filesystem.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/string.hh>

namespace Dune {
namespace XT {
namespace Common {

template <class V>
typename std::enable_if<is_vector<V>::value, void>::type
print(const V& vec, const std::string name, std::ostream& out, const std::string prefix = "")
{
  out << prefix << name << " = ";
  if (vec.size() == 0)
    out << "[]";
  else if (vec.size() == 1)
    out << vec[0];
  else {
    out << "[" << vec[0];
    for (decltype(vec.size()) ii = 1; ii < vec.size(); ++ii)
      out << " " << vec[ii];
    out << "]";
  }
} // ... print(...)

template <class M>
typename std::enable_if<is_matrix<M>::value, void>::type
print(const M& mat, const std::string name, std::ostream& out, const std::string prefix = "")
{
  typedef MatrixAbstraction<M> MM;
  const size_t rows = MM::rows(mat);
  const size_t cols = MM::cols(mat);
  out << prefix << name << " = ";
  if (rows == 0 && cols == 0)
    out << "[]";
  else if (rows == 1) {
    out << "[" << MM::get_entry(mat, 0, 0);
    for (size_t jj = 1; jj < cols; ++jj)
      out << " " << MM::get_entry(mat, 0, jj);
    out << "]";
  } else if (cols == 1) {
    out << "[" << MM::get_entry(mat, 0, 0);
    for (size_t ii = 1; ii < rows; ++ii)
      out << " " << MM::get_entry(mat, ii, 0);
    out << "]";
  } else {
    out << "[";
    auto local_prefix = prefix + " ";
    for (size_t ii = 0; ii < rows - 1; ++ii) {
      out << MM::get_entry(mat, ii, 0);
      for (size_t jj = 1; jj < cols; ++jj)
        out << " " << MM::get_entry(mat, ii, jj);
      out << ";\n" << local_prefix;
    }
    out << MM::get_entry(mat, rows - 1, 0);
    for (size_t jj = 1; jj < cols; ++jj)
      out << " " << MM::get_entry(mat, rows - 1, jj);
    out << "]";
  }
} // ... print(...)

//! useful for visualizing sparsity patterns of matrices
template <class Matrix>
void matrix_to_gnuplot_stream(const Matrix& matrix, std::ostream& stream)
{
  unsigned long nz = 0;

  const auto cols = matrix.cols();
  for (auto row : value_range(matrix.rows())) {
    for (auto col : value_range(cols)) {
      if (matrix.find(row, col))
        stream << row << "\t" << col << "\t" << matrix(row, col) << "\n";
    }
    nz += matrix.numNonZeros(row);
    stream << "#non zeros in row " << row << " " << matrix.numNonZeros(row) << " (of " << matrix.cols() << " cols)\n";
  }
  stream << "#total non zeros " << nz << " of " << matrix.rows() * matrix.cols() << " entries\n";
} // matrix_to_gnuplot_stream

//! maps 1,2,3 to x,y,z / X,Y,Z
inline std::string dim_to_axis_name(const size_t dim, const bool capitalize = false)
{
  char c = 'x';

  c += dim;
  if (capitalize)
    c -= 32;
  return std::string() += c;
} // matrix_to_gnuplot_stream

/** Outputiterator to emulate python's str.join(iterable)
 * \see http://codereview.stackexchange.com/questions/30132/comma-formatted-stl-vectors/30181#30181
 * \example std::copy(strings.begin(), strings.end(), PrefixOutputIterator<string>(ostream, ","));
 **/
template <typename T>
class PrefixOutputIterator
{
  std::ostream& ostream;
  std::string prefix;
  bool first;

public:
  typedef std::size_t difference_type;
  typedef T value_type;
  typedef T* pointer;
  typedef T reference;
  typedef std::output_iterator_tag iterator_category;

  PrefixOutputIterator(std::ostream& o, std::string const& p = "")
    : ostream(o)
    , prefix(p)
    , first(true)
  {}

  PrefixOutputIterator& operator*()
  {
    return *this;
  }
  PrefixOutputIterator& operator++()
  {
    return *this;
  }
  PrefixOutputIterator& operator++(int)
  {
    return *this;
  }

  void operator=(T const& value)
  {
    if (first) {
      ostream << value;
      first = false;
    } else {
      ostream << prefix << value;
    }
  }
};

} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_COMMON_PRINT_HH
