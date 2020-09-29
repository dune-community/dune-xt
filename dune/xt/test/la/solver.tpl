// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   René Fritze     (2014 - 2017, 2019)
//   Tobias Leibner  (2014 - 2015, 2017, 2019)

// This one has to come first (includes the config.h)!
#include <dune/xt/test/main.hxx>

#include <tuple>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/la/container.hh>
#include <dune/xt/la/solver.hh>

#include <dune/xt/test/la/container.hh>

// toggle output
// std::ostream& out = std::cout;
std::ostream& out = DXTC_LOG.devnull();

using namespace Dune;
using namespace Dune::XT;
using namespace Dune::XT::LA;

{% for T_NAME, O_TYPE, R_TYPE, S_TYPE in config.testtypes %}
struct SolverTest_{{T_NAME}} : public ::testing::Test
{
  using MatrixType = {{O_TYPE}};
  using RhsType = {{R_TYPE}};
  using SolutionType = {{S_TYPE}};

  using SolverType = Solver<MatrixType>;

  static void produces_correct_results()
  {
    const size_t dim = 10;
    const MatrixType matrix = ContainerFactory<MatrixType>::create(dim);
    const RhsType rhs = ContainerFactory<RhsType>::create(dim);
    SolutionType solution = ContainerFactory<SolutionType>::create(dim);

    // dynamic test
    // Matrix is the identity matrix, solution and rhs are initially filled with ones
    const SolverType solver(matrix);
    solver.apply(rhs, solution);
    EXPECT_TRUE(XT::Common::FloatCmp::eq(solution, rhs));

    solution *= 0;
    solver.apply(rhs, solution);
    EXPECT_TRUE(XT::Common::FloatCmp::eq(solution, rhs));

    // static tests
    std::vector<std::string> types = SolverType::types();
    if (types.size() == 0)
      DUNE_THROW(Common::Exceptions::results_are_not_as_expected, "Solver has no types!");
    for (auto type : types) {
      out << "solving with type '" << type << "' and options" << std::endl;
      Common::Configuration options = SolverType::options(type);
      Common::Configuration options2 = SolverOptions<MatrixType>::options(type);
      EXPECT_TRUE(options == options2);
      options.report(out, "  ");

      // dynamic tests
      solution = ContainerFactory<SolutionType>::create(dim);
      solver.apply(rhs, solution, type);
      EXPECT_TRUE(XT::Common::FloatCmp::eq(solution, rhs));

      solution *= 0.;
      solver.apply(rhs, solution, options);
      EXPECT_TRUE(XT::Common::FloatCmp::eq(solution, rhs));
    }
  } // ... produces_correct_results(...)
}; // struct SolverTest

TEST_F(SolverTest_{{T_NAME}}, behaves_correctly)
{
  this->produces_correct_results();
}

{% endfor %}
