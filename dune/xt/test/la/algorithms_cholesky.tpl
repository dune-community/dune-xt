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
#include <vector>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/matrix.hh>

#include <dune/xt/la/algorithms/cholesky.hh>
#include <dune/xt/la/container.hh>
#include <dune/xt/la/test/container.hh>

// toggle output
// std::ostream& out = std::cout
std::ostream& out = DXTC_LOG.devnull();

using namespace Dune;
using namespace Dune::XT;
using namespace Dune::XT::LA;

{% for T_NAME, O_TYPE, R_TYPE, S_TYPE in config.testtypes %}
struct CholeskyTest_{{T_NAME}} : public ::testing::Test
{
  typedef {{O_TYPE}} MatrixType;
  typedef {{R_TYPE}} RhsType;
  typedef {{S_TYPE}} SolutionType;

  typedef Common::MatrixAbstraction<MatrixType> M;
  typedef Common::VectorAbstraction<RhsType> RhsV;
  typedef Common::VectorAbstraction<SolutionType> SolV;
  typedef typename M::ScalarType ScalarType;

  static void produces_correct_results()
  {
    const size_t dim = 5;
    const MatrixType matrix = create<MatrixType>({ {2, -1, 0, 0, 0}, {-1, 2, -1, 0, 0}, {0, -1, 2, -1, 0}, {0, 0, -1, 2, -1}, {0, 0, 0, -1, 2} }, tridiagonal_pattern(dim, dim));
    const MatrixType expected_L = create<MatrixType>({ { 1.414213562373095,                  0,                  0,                  0,                 0},
                                                       {-0.707106781186547,  1.224744871391589,                  0,                  0,                 0},
                                                       {                 0, -0.816496580927726,  1.154700538379251,                  0,                 0},
                                                       {                 0,                  0, -0.866025403784439,  1.118033988749895,                 0},
                                                       {                 0,                  0,                  0, -0.894427190999916, 1.095445115010332} },
                                                     diagonal_pattern(dim, dim) + diagonal_pattern(dim, dim, -1));
    auto L = matrix;
    cholesky(L);
    for (size_t ii = 0; ii < dim-1; ++ii)
      M::set_entry(L, ii, ii+1, 0.);
    EXPECT_TRUE(Common::FloatCmp::eq(L, expected_L));
    auto diag = RhsV::create<dim>(dim, 0.);
    auto subdiag = SolV::create<dim-1>(dim-1, 0.);
    typedef Common::VectorAbstraction<decltype(subdiag)> SubV;
    for (size_t ii = 0; ii < dim; ++ii)
      RhsV::set_entry(diag, ii, M::get_entry(matrix, ii, ii));
    for (size_t ii = 0; ii < dim-1; ++ii)
      SubV::set_entry(subdiag, ii, M::get_entry(matrix, ii+1, ii));
    tridiagonal_ldlt(diag, subdiag);
    std::vector<ScalarType> expected_diag{2, 1.5, 4./3., 1.25, 1.2};
    std::vector<ScalarType> expected_subdiag{-0.5, -2./3., -0.75, -0.8};
    for (size_t ii = 0; ii < dim; ++ii)
      EXPECT_TRUE(Common::FloatCmp::eq(RhsV::get_entry(diag, ii), expected_diag[ii]));
    for (size_t ii = 0; ii < dim-1; ++ii)
      EXPECT_TRUE(Common::FloatCmp::eq(SubV::get_entry(subdiag, ii), expected_subdiag[ii]));

    RhsType rhs = RhsV::create(dim, 1.);
    MatrixType matrix_rhs = M::create(dim, dim, 1.);
    XT::LA::solve_tridiagonal_ldlt_factorized(diag, subdiag, rhs);
    RhsType expected_result{2.5, 4., 4.5, 4, 2.5};
    EXPECT_TRUE(Common::FloatCmp::eq(rhs, expected_result));
    XT::LA::solve_tridiagonal_ldlt_factorized(diag, subdiag, matrix_rhs);
    for (size_t jj = 0; jj < dim; ++jj) {
      for (size_t ii = 0; ii < dim; ++ii)
        RhsV::set_entry(rhs, ii, M::get_entry(matrix_rhs, ii, jj));
      EXPECT_TRUE(Common::FloatCmp::eq(rhs, expected_result));
    }
  } // ... produces_correct_results(...)
}; // struct CholeskyTest

TEST_F(CholeskyTest_{{T_NAME}}, behaves_correctly)
{
  this->produces_correct_results();
}

{% endfor %}
