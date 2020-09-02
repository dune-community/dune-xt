// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   René Fritze    (2019)
//   Tobias Leibner (2018)

// This one has to come first (includes the config.h)!
#include <dune/xt/test/main.hxx>

#include <tuple>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/matrix.hh>

#include <dune/xt/la/algorithms/qr.hh>
#include <dune/xt/la/container.hh>
#include <dune/xt/test/la/container.hh>

// toggle output
// std::ostream& out = std::cout
std::ostream& out = DXTC_LOG.devnull();

using namespace Dune;
using namespace Dune::XT;
using namespace Dune::XT::LA;
using namespace Dune::XT::Common;

{% for T_NAME, O_TYPE, R_TYPE, S_TYPE in config.testtypes %}
struct QrTest_{{T_NAME}} : public ::testing::Test
{
  typedef {{O_TYPE}} MatrixType;
  typedef {{R_TYPE}} RhsType;
  typedef {{S_TYPE}} SolutionType;

  typedef XT::Common::MatrixAbstraction<MatrixType> M;
  typedef XT::Common::VectorAbstraction<RhsType> RhsV;
  typedef XT::Common::VectorAbstraction<SolutionType> SolV;

  static constexpr size_t dim = 6;

  QrTest_{{T_NAME}}()
    : matrix_(create<MatrixType>({ {8, 1, 6,   0,     0,     0},
                                   {3, 5, 7,   0,     0,     0},
                                   {4, 9, 2,   0,     0,     0},
                                   {0, 0, 0,   1,   0.5, 1./3.},
                                   {0, 0, 0, 0.5, 1./3.,  0.25},
                                   {0, 0, 0,   0,  0.25,   0.2} }))
    , zeros_(M::create(dim, dim , 0.))
  {
  }

  void decomposes_correctly()
  {
    auto QR = matrix_;
    SolutionType tau = SolV::create(dim, 0.);
    std::vector<int> permutations(dim);
    XT::LA::qr(QR, tau, permutations);

    MatrixType R = QR;
    for (size_t ii = 0; ii < dim; ++ii)
      for (size_t jj = 0; jj < ii; ++jj)
        M::set_entry(R, ii, jj, 0.);
    auto Q = calculate_q_from_qr(QR, tau);
    MatrixType P = eye_matrix<MatrixType>(dim, dim, dense_pattern(dim, dim));
    get_permutation_matrix(permutations, P);
    EXPECT_TRUE(XT::Common::FloatCmp::eq(Q * R - matrix_ * P, zeros_, 1e-13, 1e-13));
    EXPECT_TRUE(XT::Common::FloatCmp::eq(Q * transposed(Q) - eye_matrix<MatrixType>(dim, dim), zeros_));
    EXPECT_TRUE(XT::Common::FloatCmp::eq(transposed(Q) * Q - eye_matrix<MatrixType>(dim, dim), zeros_));
  } // decomposes_correctly

  void solves_correctly()
  {
    auto QR = matrix_;
    SolutionType tau = SolV::create(dim, 0.);
    std::vector<int> permutations(dim);
    XT::LA::qr(QR, tau, permutations);
    auto Q = calculate_q_from_qr(QR, tau);

    const RhsType rhs{-1., 1., -2., 2., -3., 3.};
    SolutionType solution = SolV::create(dim, 0.);
    apply_q_from_qr<Common::Transpose::no>(QR, tau, rhs, solution);
    SolutionType other_solution(solution);
    Q.mv(rhs, other_solution);
    EXPECT_TRUE(XT::Common::FloatCmp::eq(solution, other_solution));

    solution = SolV::create(dim, 0.);
    XT::LA::solve_qr_factorized(QR, tau, permutations, solution, rhs);
    const auto expected_result = SolutionType{-151/360., -23./180., 149./360., -24, 252, -300};
    EXPECT_TRUE(XT::Common::FloatCmp::eq(solution, expected_result, 1e-14, 1e-14));
  } // ... produces_correct_results(...)

  private:
    const MatrixType matrix_;
    const MatrixType zeros_;
}; // struct QrTest

TEST_F(QrTest_{{T_NAME}}, behaves_correctly)
{
  this->decomposes_correctly();
}

TEST_F(QrTest_{{T_NAME}}, solves_correctly)
{
  this->solves_correctly();
}


{% endfor %}
