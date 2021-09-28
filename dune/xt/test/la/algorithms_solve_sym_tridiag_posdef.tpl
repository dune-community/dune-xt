// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   René Fritze    (2019)
//   Tobias Leibner (2018, 2020)

// This one has to come first (includes the config.h)!
#include <dune/xt/test/main.hxx>

#include <tuple>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/matrix.hh>

#include <dune/xt/la/algorithms/solve_sym_tridiag_posdef.hh>
#include <dune/xt/la/container.hh>
#include <dune/xt/test/la/container.hh>

// toggle output
// std::ostream& out = std::cout
std::ostream& out = DXTC_LOG.devnull();

using namespace Dune;
using namespace Dune::XT;
using namespace Dune::XT::LA;

{% for T_NAME, O_TYPE, R_TYPE, S_TYPE in config.testtypes %}
struct TriangularSolverTest_{{T_NAME}} : public ::testing::Test
{
  using MatrixType = {{O_TYPE}};
  using RhsType = {{R_TYPE}};
  using SolutionType = {{S_TYPE}};

  using M = XT::Common::MatrixAbstraction<MatrixType>;
  using RhsV = XT::Common::VectorAbstraction<RhsType>;
  using SolV = XT::Common::VectorAbstraction<SolutionType>;

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
