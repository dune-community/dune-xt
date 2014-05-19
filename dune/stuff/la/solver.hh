// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_SOLVER_HH
#define DUNE_STUFF_LA_SOLVER_HH

#include <string>
#include <vector>

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configtree.hh>

namespace Dune {
namespace Stuff {
namespace Exceptions {


class linear_solver_failed_bc_matrix_did_not_fulfill_requirements : public linear_solver_failed
{
};
class linear_solver_failed_bc_it_did_not_converge : public linear_solver_failed
{
};
class linear_solver_failed_bc_it_was_not_set_up_correctly : public linear_solver_failed
{
};
class linear_solver_failed_bc_the_solution_does_not_solve_the_system : public linear_solver_failed
{
};


} // namespace Exceptions
namespace LA {


class SolverUtils
{
public:
  static void check_given(const std::string& type, const std::vector<std::string>& opts)
  {
    if (std::find(opts.begin(), opts.end(), type) == opts.end()) {
      std::stringstream ss;
      for (auto opt : opts)
        ss << opt << " ";
      DUNE_THROW_COLORFULLY(Exceptions::configuration_error,
                            "Given type '" << type << "' not supported (see below for a list of supported ones). "
                                           << "Call options() first!\n"
                                           << ss.str());
    }
  }
};


template <class MatrixImp, class GridView = void>
class Solver
{
  static_assert(AlwaysFalse<MatrixImp>::value, "This is the unspecialized version of LA::Solver< ... >. Please include "
                                               "the correct header for your matrix implementation!");

public:
  typedef MatrixImp MatrixType;

  Solver(const MatrixType& /*matrix*/)
  {
    DUNE_THROW_COLORFULLY(NotImplemented,
                          "This is the unspecialized version of LA::Solver< ... >. "
                          "Please include the correct header for your matrix implementation '"
                              << Common::Typename<MatrixType>::value()
                              << "'!");
  }

  static std::vector<std::string> options()
  {
    DUNE_THROW_COLORFULLY(NotImplemented,
                          "This is the unspecialized version of LA::Solver< ... >. "
                          "Please include the correct header for your matrix implementation '"
                              << Common::Typename<MatrixType>::value()
                              << "'!");
  }

  static Common::ConfigTree options(const std::string& /*type*/)
  {
    DUNE_THROW_COLORFULLY(NotImplemented,
                          "This is the unspecialized version of LA::Solver< ... >. "
                          "Please include the correct header for your matrix implementation '"
                              << Common::Typename<MatrixType>::value()
                              << "'!");
  }

  /**
   *  Throws any of the above exceptions, if there was a problem. If none was thrown we beleive that a suitable solution
   *  was found (given the current options).
   */
  template <class RhsType, class SolutionType>
  void apply(const RhsType& /*rhs*/, SolutionType& /*solution*/) const
  {
    DUNE_THROW_COLORFULLY(NotImplemented,
                          "This is the unspecialized version of LA::Solver< ... >. "
                          "Please include the correct header for your matrix implementation '"
                              << Common::Typename<MatrixType>::value()
                              << "'!");
  }

  template <class RhsType, class SolutionType>
  void apply(const RhsType& /*rhs*/, SolutionType& /*solution*/, const std::string& /*type*/) const
  {
    DUNE_THROW_COLORFULLY(NotImplemented,
                          "This is the unspecialized version of LA::Solver< ... >. "
                          "Please include the correct header for your matrix implementation '"
                              << Common::Typename<MatrixType>::value()
                              << "'!");
  }

  template <class RhsType, class SolutionType>
  void apply(const RhsType& /*rhs*/, SolutionType& /*solution*/, const Common::ConfigTree& /*options*/) const
  {
    DUNE_THROW_COLORFULLY(NotImplemented,
                          "This is the unspecialized version of LA::Solver< ... >. "
                          "Please include the correct header for your matrix implementation '"
                              << Common::Typename<MatrixType>::value()
                              << "'!");
  }
}; // class Solver


} // namespace LA
} // namespace Stuff
} // namespace Dune

#include "solver/common.hh"
#include "solver/eigen.hh"
#include "solver/istl.hh"

#endif // DUNE_STUFF_LA_SOLVER_HH
