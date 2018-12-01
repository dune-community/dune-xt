// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_LA_CONTAINER_IO_HH
#define DUNE_XT_LA_CONTAINER_IO_HH

#include <fstream>
#include <utility>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/string.hh>

#include "matrix-interface.hh"
#include "vector-interface.hh"

namespace Dune {
namespace XT {
namespace LA {


template <class M>
void to_file(const MatrixInterface<M>& matrix, const std::string& filename, const std::string& mode = "ascii")
{
  if (filename.empty())
    DUNE_THROW(Common::Exceptions::wrong_input_given, "'filename' must not be empty!");
  if (mode != "ascii")
    DUNE_THROW(NotImplemented, "Currently, only 'ascii' is implemented!");
  std::ofstream file(filename);
  if (!file.is_open())
    DUNE_THROW(IOError, "Could not open '" << filename << "' for writing!");
  const auto pattern = matrix.pattern();
  for (size_t ii = 0; ii < pattern.size(); ++ii)
    for (const auto& jj : pattern.inner(ii))
      file << ii << " " << jj << " " << std::scientific << std::setprecision(15) << matrix.get_entry(ii, jj)
           << std::endl;
} // ... to_file(...)


template <class V>
void to_file(const VectorInterface<V>& vector, const std::string& filename, const std::string& mode = "ascii")
{
  if (filename.empty())
    DUNE_THROW(Common::Exceptions::wrong_input_given, "'filename' must not be empty!");
  if (mode != "ascii")
    DUNE_THROW(NotImplemented, "Currently, only 'ascii' is implemented!");
  std::ofstream file(filename);
  if (!file.is_open())
    DUNE_THROW(IOError, "Could not open '" << filename << "' for writing!");
  for (size_t ii = 0; ii < vector.size(); ++ii)
    file << ii << " " << std::scientific << std::setprecision(15) << vector[ii] << std::endl;
} // ... to_file(...)


template <class M>
typename std::enable_if<is_matrix<M>::value, M>::type from_file(const std::string& filename,
                                                                const ssize_t min_rows = -1,
                                                                const ssize_t min_cols = -1,
                                                                const std::string& mode = "ascii")
{
  if (filename.empty())
    DUNE_THROW(Common::Exceptions::wrong_input_given, "Given filename must not be empty!");
  if (mode != "ascii")
    DUNE_THROW(NotImplemented, "Currently, only 'ascii' is implemented!");
  std::ifstream file(filename);
  if (!file.is_open())
    DUNE_THROW(IOError, "Could not open '" << filename << "' for reading!");
  typedef typename M::ScalarType R;
  std::vector<std::tuple<size_t, size_t, R>> values;
  std::set<size_t> rows;
  size_t max_row = 0;
  size_t max_col = 0;
  std::string line;
  while (std::getline(file, line)) {
    const auto words = Common::tokenize(line, " ");
    if (words.size() != 3)
      DUNE_THROW(IOError,
                 "Encountered illegal line (see below), Each line has to be of the form 'ii jj value'!\n\n"
                     << "The line is '" << line << "'");
    auto ii = Common::from_string<size_t>(words[0]);
    auto jj = Common::from_string<size_t>(words[1]);
    auto value = Common::from_string<R>(words[2]);
    rows.insert(ii);
    max_row = std::max(max_row, ii);
    max_col = std::max(max_col, jj);
    values.emplace_back(std::tuple<size_t, size_t, R>(ii, jj, value));
  }
  if (rows.size() == 0)
    DUNE_THROW(IOError, "Given file '" << filename << "' must not be empty!");
  const size_t matrix_rows = std::max(min_rows, ssize_t(max_row) + 1);
  const size_t matrix_cols = std::max(min_cols, ssize_t(max_col) + 1);
  SparsityPatternDefault pattern(matrix_rows);
  for (const auto& element : values)
    pattern.insert(std::get<0>(element), std::get<1>(element));
  for (size_t ii = 0; ii < matrix_rows; ++ii)
    if (rows.count(ii) == 0)
      pattern.insert(ii, 0); // <- ensure at least one entry in each row
  pattern.sort();
  M matrix(matrix_rows, matrix_cols, pattern);
  for (const auto& element : values)
    matrix.set_entry(std::get<0>(element), std::get<1>(element), std::get<2>(element));
  return matrix;
} // ... from_file(...)


template <class V>
typename std::enable_if<is_vector<V>::value, V>::type
from_file(const std::string& filename, const ssize_t min_size = -1, const std::string& mode = "ascii")
{
  if (filename.empty())
    DUNE_THROW(Common::Exceptions::wrong_input_given, "Given filename must not be empty!");
  if (mode != "ascii")
    DUNE_THROW(NotImplemented, "Currently, only 'ascii' is implemented!");
  std::ifstream file(filename);
  if (!file.is_open())
    DUNE_THROW(IOError, "Could not open '" << filename << "' for reading!");
  typedef typename V::ScalarType R;
  std::vector<std::tuple<size_t, R>> values;
  size_t max_size = 0;
  std::string line;
  while (std::getline(file, line)) {
    const auto words = Common::tokenize(line, " ");
    if (words.size() != 2)
      DUNE_THROW(IOError,
                 "Encountered illegal line (see below), Each line has to be of the form 'ii value'!\n\n"
                     << "The line is '" << line << "'");
    auto ii = Common::from_string<size_t>(words[0]);
    auto value = Common::from_string<R>(words[1]);
    max_size = std::max(max_size, ii);
    values.emplace_back(std::tuple<size_t, R>(ii, value));
  }
  if (values.size() == 0)
    DUNE_THROW(IOError, "Given file '" << filename << "' must not be empty!");
  const size_t vector_size = std::max(min_size, ssize_t(max_size) + 1);
  V vector(vector_size);
  for (const auto& element : values)
    vector[std::get<0>(element)] = std::get<1>(element);
  return vector;
} // ... from_file(...)


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_IO_HH
