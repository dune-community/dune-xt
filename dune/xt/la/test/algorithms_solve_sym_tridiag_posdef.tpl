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

#include <dune/xt/la/algorithms/solve_sym_tridiag_posdef.hh>
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
    const size_t dim = 5;
    const MatrixType matrix = create<MatrixType>({ {2, -1, 0, 0, 0}, {-1, 2, -1, 0, 0}, {0, -1, 2, -1, 0}, {0, 0, -1, 2, -1}, {0, 0, 0, -1, 2} }, tridiagonal_pattern(dim, dim));
    const RhsType rhs = RhsV::create(dim, 1.);
    SolutionType solution = SolV::create(dim, 0.);

    XT::LA::solve_sym_tridiag_posdef(matrix, solution, rhs);
    SolutionType expected_result{2.5, 4., 4.5, 4, 2.5};
    EXPECT_TRUE(XT::Common::FloatCmp::eq(solution, expected_result));
  } // ... produces_correct_results(...)
}; // struct TriangularSolverTest

TEST_F(TriangularSolverTest_{{T_NAME}}, behaves_correctly)
{
  this->produces_correct_results();
}

{% endfor %}
