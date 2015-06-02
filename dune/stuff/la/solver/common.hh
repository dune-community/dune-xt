// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_SOLVER_COMMON_HH
#define DUNE_STUFF_LA_SOLVER_COMMON_HH

#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cmath>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configuration.hh>

#include <dune/stuff/la/container/common.hh>

#include "../solver.hh"

namespace Dune {
namespace Stuff {
namespace LA {


template <class S, class CommunicatorType>
class Solver<CommonDenseMatrix<S>, CommunicatorType> : protected SolverUtils
{
public:
  typedef CommonDenseMatrix<S> MatrixType;
  typedef typename MatrixType::RealType R;

  Solver(const MatrixType& matrix)
    : matrix_(matrix)
  {
  }

  Solver(const MatrixType& matrix, const CommunicatorType& /*communicator*/)
    : matrix_(matrix)
  {
  }

  static std::vector<std::string> types()
  {
    return {"superlu"};
  }

  static Common::Configuration options(const std::string type = "")
  {
    const std::string tp = !type.empty() ? type : types()[0];
    SolverUtils::check_given(tp, types());
    return Common::Configuration({"type", "post_check_solves_system"}, {tp, "1e-5"});
  } // ... options(...)

  void apply(const CommonDenseVector<S>& rhs, CommonDenseVector<S>& solution) const
  {
    apply(rhs, solution, types()[0]);
  }

  void apply(const CommonDenseVector<S>& rhs, CommonDenseVector<S>& solution, const std::string& type) const
  {
    apply(rhs, solution, options(type));
  }

  void apply(const CommonDenseVector<S>& rhs, CommonDenseVector<S>& solution, const Common::Configuration& opts) const
  {
    if (!opts.has_key("type"))
      DUNE_THROW(Exceptions::configuration_error,
                 "Given options (see below) need to have at least the key 'type' set!\n\n" << opts);
    const auto type = opts.get<std::string>("type");
    SolverUtils::check_given(type, types());
    const Common::Configuration default_opts = options(type);
    // solve
    try {
      matrix_.backend().solve(solution.backend(), rhs.backend());
    } catch (FMatrixError&) {
      DUNE_THROW(Exceptions::linear_solver_failed_bc_data_did_not_fulfill_requirements,
                 "The dune-common backend reported 'FMatrixError'!\n"
                     << "Those were the given options:\n\n"
                     << opts);
    }
    // check
    const R post_check_solves_system_threshold =
        opts.get("post_check_solves_system", default_opts.get<R>("post_check_solves_system"));
    if (post_check_solves_system_threshold > 0) {
      auto tmp = rhs.copy();
      matrix_.mv(solution, tmp);
      tmp -= rhs;
      const R sup_norm = tmp.sup_norm();
      if (sup_norm > post_check_solves_system_threshold || std::isnan(sup_norm) || std::isinf(sup_norm))
        DUNE_THROW(Exceptions::linear_solver_failed_bc_the_solution_does_not_solve_the_system,
                   "The computed solution does not solve the system (although the dune-common backend "
                       << "reported no error) and you requested checking (see options below)! "
                       << "If you want to disable this check, set 'post_check_solves_system = 0' in the options."
                       << "\n\n"
                       << "  (A * x - b).sup_norm() = "
                       << tmp.sup_norm()
                       << "\n\n"
                       << "Those were the given options:\n\n"
                       << opts);
    }
  } // ... apply(...)

private:
  const MatrixType& matrix_;
}; // class Solver< CommonDenseMatrix< ... > >


} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_SOLVER_COMMON_HH
