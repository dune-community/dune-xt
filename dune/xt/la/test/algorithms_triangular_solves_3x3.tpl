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

#include <dune/xt/la/algorithms/triangular_solves.hh>
#include <dune/xt/la/container.hh>
#include <dune/xt/la/test/container.hh>

// toggle output
// std::ostream& out = std::cout
std::ostream& out = DXTC_LOG.devnull();

using namespace Dune;
using namespace Dune::XT;
using namespace Dune::XT::LA;

{% for T_NAME, O_TYPE, R_TYPE, S_TYPE in config.testtypes %}
struct TriangularSolverTest_{{T_NAME}} : public ::testing::Test
{
  typedef {{O_TYPE}} MatrixType;
  typedef {{R_TYPE}} RhsType;
  typedef {{S_TYPE}} SolutionType;

  typedef XT::Common::MatrixAbstraction<MatrixType> M;
  typedef XT::Common::VectorAbstraction<RhsType> RhsV;
  typedef XT::Common::VectorAbstraction<SolutionType> SolV;

  static void produces_correct_results()
  {
    const size_t dim = 3;
    const MatrixType lower_triangular_matrix = create<MatrixType>({ {1, 0, 0}, {4, 5, 0}, {7, 8, 9} }, triangular_pattern(dim, dim, Common::MatrixPattern::lower_triangular));
    const MatrixType upper_triangular_matrix = create<MatrixType>({ {1, 2, 3}, {0, 5, 6}, {0, 0, 9} }, triangular_pattern(dim, dim, Common::MatrixPattern::upper_triangular));
    const RhsType rhs = RhsV::create(dim, 1.);
    SolutionType solution = SolV::create(dim, 0.);

    XT::LA::solve_lower_triangular(lower_triangular_matrix, solution, rhs);
    SolutionType expected_result{1, -3./5., -2./15.};
    EXPECT_TRUE(XT::Common::FloatCmp::eq(solution, expected_result));
    expected_result = SolutionType{2./15., 1./45., 1./9.};
    XT::LA::solve_lower_triangular_transposed(lower_triangular_matrix, solution, rhs);
    EXPECT_TRUE(XT::Common::FloatCmp::eq(solution, expected_result));
    expected_result = SolutionType{8./15., 1./15., 1./9.};
    XT::LA::solve_upper_triangular(upper_triangular_matrix, solution, rhs);
    EXPECT_TRUE(XT::Common::FloatCmp::eq(solution, expected_result));
    expected_result = SolutionType{1., -0.2, -4./45.};
    XT::LA::solve_upper_triangular_transposed(upper_triangular_matrix, solution, rhs);
    EXPECT_TRUE(XT::Common::FloatCmp::eq(solution, expected_result));
  } // ... produces_correct_results(...)
}; // struct TriangularSolverTest

TEST_F(TriangularSolverTest_{{T_NAME}}, behaves_correctly)
{
  this->produces_correct_results();
}

{% endfor %}
