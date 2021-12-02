// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   René Fritze    (2020)
//   Tobias Leibner (2019 - 2020)

#ifndef DUNE_XT_LA_CONTAINER_MATRIX_MARKET_HH
#define DUNE_XT_LA_CONTAINER_MATRIX_MARKET_HH

#include <fstream>
#include <iostream>
#include <string>

#include <boost/algorithm/string/constants.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/la/container/pattern.hh>

namespace Dune {
namespace XT {
namespace LA {
namespace internal {


template <class MatrixType>
MatrixType read_matrix_market_array_format(std::ifstream& matrix_file,
                                           std::string& curr_line,
                                           const bool field_qualifier_is_complex,
                                           const bool is_general,
                                           const bool is_symmetric,
                                           const bool is_skew_symmetric)
{
  using M = XT::Common::MatrixAbstraction<MatrixType>;
  using RealType = typename M::RealType;
  using ScalarType = typename M::ScalarType;
  // dense matrix format, first line contains 'num_rows num_cols', each following line contains one matrix entry,
  // entries are in column-major order
  // parse first line
  std::vector<std::string> matrix_dimensions;
  boost::algorithm::split(
      matrix_dimensions, curr_line, boost::algorithm::is_space(), boost::algorithm::token_compress_on);
  DUNE_THROW_IF(matrix_dimensions.size() != 2, Dune::IOError, "Could not read matrix dimensions!");
  const size_t rows = XT::Common::from_string<size_t>(matrix_dimensions[0]);
  const size_t cols = XT::Common::from_string<size_t>(matrix_dimensions[1]);
  MatrixType ret = M::create(rows, cols, 0., XT::LA::dense_pattern(rows, cols));
  DUNE_THROW_IF(!is_general && rows != cols,
                Dune::InvalidStateException,
                "You're trying to create a non-square symmetric/skew-symmetric/hermitian matrix!");
  // read entries
  const size_t expected_entry_size = field_qualifier_is_complex ? 2 : 1;
  std::vector<std::string> curr_entry(expected_entry_size);
  const size_t expected_num_entries =
      is_general ? rows * cols : (is_skew_symmetric ? (rows * (rows - 1)) / 2 : (rows * (rows - 1)) / 2 + rows);
  for (size_t ii = 0; ii < expected_num_entries; ++ii) {
    size_t curr_row, curr_col;
    if (is_general) {
      curr_row = ii % rows;
      curr_col = ii / rows;
    } else {
      // (skew-)symmetric or hermitian
      // only (strictly) lower triangular part is given
      // In the symmetric/hermitian case, first column has rows entries, second one rows - 1 entries, ...
      // In the skew-symmetric case, the diagonal is also ommited, so each column has one entry less
      size_t col_size = is_skew_symmetric ? rows - 1 : rows;
      curr_col = 0;
      curr_row = ii;
      while (curr_row >= col_size) {
        curr_row -= col_size;
        curr_col += 1;
        col_size -= 1;
      }
      curr_row += is_skew_symmetric ? curr_col + 1 : curr_col;
    }
    // get next line, skip blank lines
    do {
      std::getline(matrix_file, curr_line);
      XT::Common::trim(curr_line);
    } while (curr_line.empty() && matrix_file.good());
    DUNE_THROW_IF(!matrix_file.good(), Dune::IOError, "There were not enough entries for this matrix!");
    boost::algorithm::split(curr_entry, curr_line, boost::algorithm::is_space(), boost::algorithm::token_compress_on);
    DUNE_THROW_IF(
        !(curr_entry.size() == expected_entry_size), Dune::InvalidStateException, "Invalid entry encountered!");
    ScalarType entry = XT::Common::create_real_or_complex_number<ScalarType>(
        XT::Common::from_string<RealType>(curr_entry[0]),
        XT::Common::from_string<RealType>(field_qualifier_is_complex ? curr_entry[1] : "0."));
    M::set_entry(ret, curr_row, curr_col, entry);
  } // ii
  if (!is_general) {
    // fill in upper triangular part
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < ii; ++jj) {
        if (is_symmetric)
          M::set_entry(ret, jj, ii, M::get_entry(ret, ii, jj));
        else if (is_skew_symmetric)
          M::set_entry(ret, jj, ii, -1. * M::get_entry(ret, ii, jj));
        else // format is hermitian
          M::set_entry(ret, jj, ii, XT::Common::conj(M::get_entry(ret, ii, jj)));
      } // jj
    } // ii
  }
  return ret;
}

template <class MatrixType>
MatrixType read_matrix_market_coordinate_format(std::ifstream& matrix_file,
                                                std::string& curr_line,
                                                const bool field_qualifier_is_complex,
                                                const bool is_general,
                                                const bool is_symmetric,
                                                const bool is_skew_symmetric)
{
  using M = XT::Common::MatrixAbstraction<MatrixType>;
  using RealType = typename M::RealType;
  using ScalarType = typename M::ScalarType;
  // coordinate format, first line contains 'num_rows num_cols num_nonzeros', each following line contains one matrix
  // entry in the format 'row_index col_index entry', indices are 1-based
  std::vector<std::string> matrix_dimensions;
  boost::algorithm::split(
      matrix_dimensions, curr_line, boost::algorithm::is_space(), boost::algorithm::token_compress_on);
  DUNE_THROW_IF(matrix_dimensions.size() != 3, Dune::IOError, "Could not read matrix dimensions!");
  const size_t rows = XT::Common::from_string<size_t>(matrix_dimensions[0]);
  const size_t cols = XT::Common::from_string<size_t>(matrix_dimensions[1]);
  const size_t nnz = XT::Common::from_string<size_t>(matrix_dimensions[2]);
  // read entries
  std::vector<std::tuple<size_t, size_t, ScalarType>> entries(nnz);
  const size_t expected_entry_size = field_qualifier_is_complex ? 4 : 3;
  std::vector<std::string> curr_entry(expected_entry_size);
  for (size_t ii = 0; ii < nnz; ++ii) {
    // get next line, skip blank lines
    do {
      std::getline(matrix_file, curr_line);
      XT::Common::trim(curr_line);
    } while (curr_line.empty() && matrix_file.good());
    DUNE_THROW_IF(!matrix_file.good(), Dune::IOError, "There were not enough entries for this matrix!");
    boost::algorithm::split(curr_entry, curr_line, boost::algorithm::is_space(), boost::algorithm::token_compress_on);
    DUNE_THROW_IF(
        !(curr_entry.size() == expected_entry_size), Dune::InvalidStateException, "Invalid entry encountered!");
    ScalarType entry = XT::Common::create_real_or_complex_number<ScalarType>(
        XT::Common::from_string<RealType>(curr_entry[2]),
        XT::Common::from_string<RealType>(field_qualifier_is_complex ? curr_entry[3] : "0."));
    entries[ii] = std::make_tuple(
        XT::Common::from_string<size_t>(curr_entry[0]), XT::Common::from_string<size_t>(curr_entry[1]), entry);
  }
  SparsityPatternDefault pattern(rows);
  for (size_t ii = 0; ii < nnz; ++ii) {
    // entries in matrix market format are 1-based
    const size_t row = std::get<0>(entries[ii]);
    const size_t col = std::get<1>(entries[ii]);
    DUNE_THROW_IF(
        row == 0 || col == 0, Dune::InvalidStateException, "Indices in matrix market format have to be 1-based!");
    DUNE_THROW_IF(
        (!is_general) && row < col,
        Dune::InvalidStateException,
        "Only provide (strictly) lower triangular portion of matrix for symmetric/hermitian/skew-symmetric format!");
    DUNE_THROW_IF(
        is_skew_symmetric && row == col,
        Dune::InvalidStateException,
        "Only provide (strictly) lower triangular portion of matrix for symmetric/hermitian/skew-symmetric format!");
    pattern.insert(row - 1, col - 1);
    if (!is_general)
      pattern.insert(col - 1, row - 1);
  }
  pattern.sort();
  MatrixType ret = M::create(rows, cols, 0., pattern);
  for (size_t ii = 0; ii < nnz; ++ii) {
    const size_t row = std::get<0>(entries[ii]) - 1;
    const size_t col = std::get<1>(entries[ii]) - 1;
    const ScalarType entry = std::get<2>(entries[ii]);
    M::set_entry(ret, row, col, entry);
    if (!is_general) {
      if (is_symmetric)
        M::set_entry(ret, col, row, entry);
      else if (is_skew_symmetric)
        M::set_entry(ret, col, row, -entry);
      else // format is hermitian
        M::set_entry(ret, col, row, XT::Common::conj(entry));
    }
  } // ii
  return ret;
}


} // namespace internal


template <class MatrixType>
MatrixType read_matrix_market(const std::string& filename)
{
  using M = XT::Common::MatrixAbstraction<MatrixType>;
  using ScalarType = typename M::ScalarType;
  constexpr bool scalartype_is_complex = XT::Common::is_complex<ScalarType>::value;
  std::ifstream matrix_file(filename);
  DUNE_THROW_IF(!matrix_file.is_open(), Dune::IOError, "Opening matrix file for reading failed!");
  std::string curr_line;
  static const std::string matrix_market_prefix = "%%MatrixMarket";
  // Search for matrix market header
  // Find line that starts with %%MatrixMarket, compare returns 0 if match is found
  while (curr_line.compare(0, matrix_market_prefix.size(), matrix_market_prefix) && matrix_file.good()) {
    std::getline(matrix_file, curr_line);
    XT::Common::trim(curr_line);
  }
  DUNE_THROW_IF(!matrix_file.good(), Dune::IOError, "File is not a valid matrix market file!");
  std::vector<std::string> tokens;
  boost::algorithm::split(tokens, curr_line, boost::algorithm::is_space(), boost::algorithm::token_compress_on);
  DUNE_THROW_IF(
      !(tokens[0] == matrix_market_prefix), Dune::IOError, "There has to be whitespace after %%MatrixMarket!");
  const std::string object_str = XT::Common::to_lower(tokens[1]);
  const std::string format_str = XT::Common::to_lower(tokens[2]);
  const std::string field_qualifier = XT::Common::to_lower(tokens[3]);
  const std::string symmetry_qualifier = tokens.size() >= 5 ? XT::Common::to_lower(tokens[4]) : "";
  const bool field_qualifier_is_complex = (field_qualifier == "complex");
  const bool is_general = (symmetry_qualifier == "general");
  const bool is_symmetric = (symmetry_qualifier == "symmetric");
  const bool is_skew_symmetric = (symmetry_qualifier == "skew-symmetric");
  DUNE_THROW_IF(
      !(object_str == "matrix"), Dune::NotImplemented, "Only matrix market matrix files are supported by now!");
  DUNE_THROW_IF(!(format_str == "coordinate" || format_str == "array"),
                Dune::NotImplemented,
                "Only coordinate and array format are supported!");
  DUNE_THROW_IF(!(field_qualifier == "real" || field_qualifier == "complex"),
                Dune::NotImplemented,
                "Only real and complex matrices are supported by now!");
  DUNE_THROW_IF(field_qualifier_is_complex && !scalartype_is_complex,
                Dune::InvalidStateException,
                "You are trying to read a complex matrix into a real matrix type!");
  DUNE_THROW_IF(!is_general && !is_symmetric && !is_skew_symmetric && (symmetry_qualifier != "hermitian"),
                Dune::NotImplemented,
                "Symmetry qualifier has to be 'general', 'symmetric', 'skew-symmetric' or 'hermitian'");

  // Skip comments (comments start with %)
  std::getline(matrix_file, curr_line);
  XT::Common::trim(curr_line);
  while ((!(curr_line.compare(0, 1, "%")) || curr_line.empty()) && matrix_file.good()) {
    std::getline(matrix_file, curr_line);
    XT::Common::trim(curr_line);
  }
  DUNE_THROW_IF(!matrix_file.good(), Dune::IOError, "File only contains header and comments, no data found!");
  if (format_str == "array")
    return internal::read_matrix_market_array_format<MatrixType>(
        matrix_file, curr_line, field_qualifier_is_complex, is_general, is_symmetric, is_skew_symmetric);
  return internal::read_matrix_market_coordinate_format<MatrixType>(
      matrix_file, curr_line, field_qualifier_is_complex, is_general, is_symmetric, is_skew_symmetric);
}


template <class MatrixType>
void write_matrix_market(const MatrixType& mat, const std::string& filename, const int precision = 20)
{
  using M = XT::Common::MatrixAbstraction<MatrixType>;
  using ScalarType = typename M::ScalarType;
  constexpr bool scalartype_is_complex = XT::Common::is_complex<ScalarType>::value;
  std::ofstream matrix_file(filename);
  DUNE_THROW_IF(!matrix_file.is_open(), Dune::IOError, "Opening matrix file for writing failed!");
  matrix_file.precision(precision);
  static constexpr bool is_dense = (M::storage_layout == Common::StorageLayout::dense_column_major
                                    || M::storage_layout == Common::StorageLayout::dense_row_major);
  std::string format_str = is_dense ? "array" : "coordinate";
  std::string field_qualifier = scalartype_is_complex ? "complex" : "real";
  std::string symmetry_qualifier = "general";
  // Write header
  matrix_file << "%%MatrixMarket matrix " << format_str << " " << field_qualifier << " " << symmetry_qualifier
              << std::endl;
  // Write dimensions
  matrix_file << mat.rows() << " " << mat.cols() << (is_dense ? "" : (" " + XT::Common::to_string(mat.non_zeros())))
              << std::endl;
  // Write data
  if (is_dense) {
    // write all entries in column-major order
    for (size_t jj = 0; jj < mat.cols(); ++jj) {
      for (size_t ii = 0; ii < mat.rows(); ++ii) {
        const ScalarType& entry = M::get_entry(mat, ii, jj);
        const std::string entry_str = scalartype_is_complex ? (XT::Common::to_string(std::real(entry), precision) + " "
                                                               + XT::Common::to_string(std::imag(entry), precision))
                                                            : XT::Common::to_string(entry, precision);
        matrix_file << entry_str << std::endl;
      } // ii
    } // jj
  } else {
    // write non-zero entries in coordinate format
    const auto& pattern = mat.pattern();
    for (size_t ii = 0; ii < mat.rows(); ++ii) {
      for (const auto& jj : pattern.inner(ii)) {
        const ScalarType& entry = M::get_entry(mat, ii, jj);
        const std::string entry_str = scalartype_is_complex ? (XT::Common::to_string(std::real(entry), precision) + " "
                                                               + XT::Common::to_string(std::imag(entry), precision))
                                                            : XT::Common::to_string(entry, precision);
        matrix_file << ii + 1 << " " << jj + 1 << " " << entry_str << std::endl;
      } // ii
    } // jj
  }
} // ... write_matrix_market(...)


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_MATRIX_MARKET_HH
