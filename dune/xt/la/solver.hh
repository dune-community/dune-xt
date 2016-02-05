// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2014, 2016)
//   Rene Milk       (2013, 2015)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_LA_SOLVER_HH
#define DUNE_XT_LA_SOLVER_HH

#include <string>
#include <vector>

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/parallel/helper.hh>

namespace Dune {
namespace XT {
namespace LA {
namespace Exceptions {

class linear_solver_failed_bc_data_did_not_fulfill_requirements : public Common::Exceptions::linear_solver_failed
{
};
class linear_solver_failed_bc_it_did_not_converge : public Common::Exceptions::linear_solver_failed
{
};
class linear_solver_failed_bc_it_was_not_set_up_correctly : public Common::Exceptions::linear_solver_failed
{
};
class linear_solver_failed_bc_the_solution_does_not_solve_the_system : public Common::Exceptions::linear_solver_failed
{
};

} // namespace Exceptions
namespace internal {

static const constexpr size_t max_size_to_print = 5;

} // namespace internal

class SolverUtils
{
public:
  static void check_given(const std::string& type, const std::vector<std::string>& opts)
  {
    if (std::find(opts.begin(), opts.end(), type) == opts.end()) {
      std::stringstream ss;
      for (auto opt : opts)
        ss << opt << " ";
      DUNE_THROW(Common::Exceptions::configuration_error,
                 "Given type '" << type << "' not supported (see below for a list of supported ones). "
                                << "Call options() first!\n"
                                << ss.str());
    }
  }
};

template <class MatrixImp, class CommunicatorType = SequentialCommunication>
class Solver
{
  static_assert(AlwaysFalse<MatrixImp>::value, "This is the unspecialized version of LA::Solver< ... >. Please include "
                                               "the correct header for your matrix implementation!");

public:
  typedef MatrixImp MatrixType;

  Solver(const MatrixType& /*matrix*/)
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::Solver< ... >. "
               "Please include the correct header for your matrix implementation '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
  }

  static std::vector<std::string> types()
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::Solver< ... >. "
               "Please include the correct header for your matrix implementation '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
  }

  static Common::Configuration options(const std::string /*type*/ = "")
  {
    DUNE_THROW(NotImplemented,
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
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::Solver< ... >. "
               "Please include the correct header for your matrix implementation '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
  }

  template <class RhsType, class SolutionType>
  void apply(const RhsType& /*rhs*/, SolutionType& /*solution*/, const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::Solver< ... >. "
               "Please include the correct header for your matrix implementation '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
  }

  template <class RhsType, class SolutionType>
  void apply(const RhsType& /*rhs*/, SolutionType& /*solution*/, const Common::Configuration& /*options*/) const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::Solver< ... >. "
               "Please include the correct header for your matrix implementation '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
  }
}; // class Solver

} // namespace LA
} // namespace XT
} // namespace Dune

#include "solver/common.hh"
#include "solver/eigen.hh"
#include "solver/istl.hh"

#endif // DUNE_XT_LA_SOLVER_HH
