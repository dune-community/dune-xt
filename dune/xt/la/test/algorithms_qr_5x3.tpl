// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014 - 2016)
//   Tobias Leibner  (2014 - 2015, 2017)

// This one has to come first (includes the config.h)!
#include <dune/xt/common/test/main.hxx>

#include <tuple>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/matrix.hh>

#include <dune/xt/la/algorithms/qr.hh>
#include <dune/xt/la/container.hh>
#include <dune/xt/la/test/container.hh>

// toggle output
// std::ostream& out = std::cout
std::ostream& out = DXTC_LOG.devnull();

using namespace Dune;
using namespace Dune::XT;
using namespace Dune::XT::LA;

{% for T_NAME, O_TYPE, R_TYPE, S_TYPE in config.testtypes %}
struct QrTest_5x3_{{T_NAME}} : public ::testing::Test
{
  typedef {{O_TYPE}} MatrixType;
  typedef {{S_TYPE}} VectorType;

  typedef XT::Common::MatrixAbstraction<MatrixType> M;
  typedef XT::Common::VectorAbstraction<VectorType> V;

  static constexpr size_t num_rows = 5;
  static constexpr size_t num_cols = 3;

  QrTest_5x3_{{T_NAME}}()
    : matrix_(create<MatrixType>({ { 1,  2,  3},
                                   { 4,  5,  6},
                                   { 7,  8,  9},
                                   {10, 11, 12},
                                   {13, 14, 15} }))
  {
  }

  void decomposes_correctly()
  {
    auto QR = matrix_;
    auto tau = V::template create<num_cols>(num_cols, 0.);
    std::vector<int> permutations(num_cols);
    XT::LA::qr(QR, tau, permutations);
    MatrixType R = QR;
    for (size_t ii = 0; ii < num_rows; ++ii)
      for (size_t jj = 0; jj < std::min(ii, num_cols); ++jj)
        M::set_entry(R, ii, jj, 0.);
    auto Q = calculate_q_from_qr(QR, tau);
    auto P = eye_matrix<typename M::MatrixTypeTemplate<M::static_cols, M::static_cols>>(num_cols, num_cols, dense_pattern(num_cols, num_cols));
    get_permutation_matrix(permutations, P);
    EXPECT_TRUE(XT::Common::FloatCmp::eq(Q * R - matrix_ * P, Common::zeros_like(matrix_), 1e-13, 1e-13));
    EXPECT_TRUE(XT::Common::FloatCmp::eq(Q * Common::transposed(Q) - eye_matrix<decltype(Q)>(num_rows, num_rows), Common::zeros_like(Q)));
    EXPECT_TRUE(XT::Common::FloatCmp::eq(Common::transposed(Q) * Q - eye_matrix<decltype(Q)>(num_rows, num_rows), Common::zeros_like(Q)));
  } // decomposes_correctly

  private:
    const MatrixType matrix_;
}; // struct QrTest

constexpr size_t QrTest_5x3_{{T_NAME}}::num_rows;
constexpr size_t QrTest_5x3_{{T_NAME}}::num_cols;

TEST_F(QrTest_5x3_{{T_NAME}}, behaves_correctly)
{
  this->decomposes_correctly();
}

{% endfor %}
