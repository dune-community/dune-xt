// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// This one has to come first (includes the config.h)!
#include "main.hxx"

#include <tuple>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>

#include "la_container.hh"

// toggle output
// std::ostream& out = std::cout;
std::ostream& out = DSC_LOG.devnull();

using namespace Dune::Stuff;
using namespace Dune::Stuff::LA;

struct SolverTest : public ::testing::Test
{
  typedef TESTMATRIXTYPE MatrixType;
  typedef TESTRHSVECTORTYPE RhsType;
  typedef TESTSOLUTIONVECTORTYPE SolutionType;

  typedef Solver<MatrixType> SolverType;

  static void produces_correct_results()
  {
    const size_t dim        = 10;
    const MatrixType matrix = ContainerFactory<MatrixType>::create(dim);
    const RhsType rhs       = ContainerFactory<RhsType>::create(dim);
    SolutionType solution = ContainerFactory<SolutionType>::create(dim);
    solution.scal(0);

    // dynamic test
    const SolverType solver(matrix);
    solver.apply(rhs, solution);
    EXPECT_TRUE(solution.almost_equal(rhs));
    solution.scal(0);

    // static tests
    std::vector<std::string> types = SolverType::types();
    if (types.size() == 0)
      DUNE_THROW(Exceptions::results_are_not_as_expected, "Solver has no types!");
    for (auto type : types) {
      out << "solving with type '" << type << "' and options" << std::endl;
      Common::Configuration options = SolverType::options(type);
      options.report(out, "  ");

      // dynamic tests
      solver.apply(rhs, solution, type);
      EXPECT_TRUE(solution.almost_equal(rhs));
      solution.scal(0);

      solver.apply(rhs, solution, options);
      EXPECT_TRUE(solution.almost_equal(rhs));
    }
  } // ... produces_correct_results(...)
}; // struct SolverTest

TEST_F(SolverTest, behaves_correctly)
{
  this->produces_correct_results();
}
