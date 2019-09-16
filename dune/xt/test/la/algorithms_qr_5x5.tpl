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

  static const size_t dim = 5;

  QrTest_{{T_NAME}}()
    : matrix_(create<MatrixType>({ {17, 24,  1,  8, 15},  // Matlab's magic(5) matrix
                                   {23,  5,  7, 14, 16},
                                   { 4,  6, 13, 20, 22},
                                   {10, 12, 19, 21,  3},
                                   {11, 18, 25,  2,  9} }))
    , expected_QR_(create<MatrixType>({ {-34.7131099154195700, -20.0212542665697800, -23.4781614780638500, -23.4781614780638500, -20.021254266569780},
                                        {  0.1960064529966254, -25.5763441014027800, -11.7271553020608900, -15.4415250971869000, -17.170138775591260},
                                        {  0.3640119841365900,  -0.2293367407101866, -20.4021999336175400,   3.0540446377429140,  -9.488235385725829},
                                        {  0.5320175152765545,  -0.2346556131453573,   0.2645095781511115,  17.4930857093189500,   3.646133447917943},
                                        {  0.7000230464165191,  -0.3609814395253336,  -0.3881556302466347,   0.6310170951667101,  16.000462873471270} }))
    , expected_Q_(create<MatrixType>({ {-0.0288075600957838, -0.6421260383974130,  0.0101293156705570,  0.7647198987584700,  0.0441038399717481},
                                       {-0.2016529206704870, -0.7414138442334970, -0.0279815350015385, -0.6343959919157710,  0.0800023143673562},
                                       {-0.3744982812451900,  0.1367640854898990, -0.6279382905586910,  0.0707174553320426,  0.6646346116672680},
                                       {-0.5473436418198930,  0.0374762796538158, -0.4209775305089780,  0.0579516211542792, -0.7200208293062070},
                                       {-0.7201890023945970,  0.1336816209296280,  0.6539004914509510,  0.0662257729361626,  0.1774410305840080} }))
    , zeros_(M::create(dim, dim , 0.))
    , expected_tau_{1.02880756009578, 1.61555299707044, 1.638498004633, 1.43042835508717, 0.}
    , expected_permutations_{2, 0, 3, 1, 4}
  {
  }

  void decomposes_correctly()
  {
    auto QR = matrix_;
    SolutionType tau = SolV::create(dim, 0.);
    std::vector<int> permutations(dim);
    XT::LA::qr(QR, tau, permutations);
    EXPECT_TRUE(XT::Common::FloatCmp::eq(QR, expected_QR_));
    EXPECT_TRUE(XT::Common::FloatCmp::eq(tau, expected_tau_, 1e-14, 1e-14));
    EXPECT_TRUE(XT::Common::FloatCmp::eq(permutations, expected_permutations_));

    MatrixType R = QR;
    for (size_t ii = 0; ii < dim; ++ii)
      for (size_t jj = 0; jj < ii; ++jj)
        M::set_entry(R, ii, jj, 0.);
    auto Q = calculate_q_from_qr(QR, tau);
    EXPECT_TRUE(XT::Common::FloatCmp::eq(Q, expected_Q_));
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

    const RhsType rhs{-1., 1., -2., 2., 3.};
    SolutionType solution = SolV::create(dim, 0.);
    apply_q_from_qr<Common::Transpose::no>(QR, tau, rhs, solution);
    SolutionType other_solution(solution);
    Q.mv(rhs, other_solution);
    EXPECT_TRUE(XT::Common::FloatCmp::eq(solution, other_solution));

    solution = SolV::create(dim, 0.);
    XT::LA::solve_qr_factorized(QR, tau, permutations, solution, rhs);
    const auto expected_result = SolutionType{0.139358974358974, -0.041410256410256, 0.142564102564103, -0.056794871794872, -0.137564102564103};
    EXPECT_TRUE(XT::Common::FloatCmp::eq(solution, expected_result));
  } // ... produces_correct_results(...)

  private:
    const MatrixType matrix_;
    const MatrixType expected_QR_;
    const MatrixType expected_Q_;
    const MatrixType zeros_;
    const SolutionType expected_tau_;
    const std::vector<int> expected_permutations_;
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
