// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright Holders: Felix Albrecht, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_SOLVER_HH
#define DUNE_STUFF_LA_SOLVER_HH

#include <string>
#include <vector>

#include <dune/stuff/common/type_utils.hh>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configtree.hh>

#include "container/interfaces.hh"
#include "container/dunedynamic.hh"
#include "container/eigen.hh"

namespace Dune {
namespace Stuff {
namespace LA {


class SolverUtils
{
protected:
  static void check_given(const std::string& type, const std::vector<std::string>& opts)
  {
    if (std::find(opts.begin(), opts.end(), type) == opts.end()) {
      std::stringstream ss;
      for (auto opt : opts)
        ss << opt << " ";
      DUNE_THROW_COLORFULLY(Exception::configuration_error,
                            "Given type '" << type << "' not supported (see below for a list of supported ones). "
                                           << "Call options() first!\n"
                                           << ss.str());
    }
  }
};


template <class MatrixImp>
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
   *  \return Returns 0 if the solve was successfull. Otherwise
   *          * 1 if the given matrix did not fulfill the requirements
   *          * 2 if an iterative solver did not converge
   *          * 3 if the solver was not set up correctly
   *          * 4 if the computed solution does not solve the system (only if NDEBUG is not defined)
   */
  template <class RhsType, class SolutionType>
  size_t apply(const RhsType& /*rhs*/, SolutionType& /*solution*/) const
  {
    DUNE_THROW_COLORFULLY(NotImplemented,
                          "This is the unspecialized version of LA::Solver< ... >. "
                          "Please include the correct header for your matrix implementation '"
                              << Common::Typename<MatrixType>::value()
                              << "'!");
  }

  template <class RhsType, class SolutionType>
  size_t apply(const RhsType& /*rhs*/, SolutionType& /*solution*/, const std::string& /*type*/) const
  {
    DUNE_THROW_COLORFULLY(NotImplemented,
                          "This is the unspecialized version of LA::Solver< ... >. "
                          "Please include the correct header for your matrix implementation '"
                              << Common::Typename<MatrixType>::value()
                              << "'!");
  }

  template <class RhsType, class SolutionType>
  size_t apply(const RhsType& /*rhs*/, SolutionType& /*solution*/, const Common::ConfigTree& /*options*/) const
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

#include "solver/dunedynamic.hh"
#include "solver/eigen.hh"

#endif // DUNE_STUFF_LA_SOLVER_HH
