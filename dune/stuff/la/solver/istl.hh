// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_SOLVER_ISTL_HH
#define DUNE_STUFF_LA_SOLVER_ISTL_HH

#include <type_traits>
#include <cmath>

#if HAVE_DUNE_ISTL
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#endif // HAVE_DUNE_ISTL

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/la/container/istl.hh>
#include <dune/stuff/la/solver/istl_amg.hh>

#include "../solver.hh"

namespace Dune {
namespace Stuff {
namespace LA {

#if HAVE_DUNE_ISTL


template <class S, class CommunicatorType>
class Solver<IstlRowMajorSparseMatrix<S>, CommunicatorType> : protected SolverUtils
{
public:
  typedef IstlRowMajorSparseMatrix<S> MatrixType;

  Solver(const MatrixType& matrix)
    : matrix_(matrix)
    , communicator_(new CommunicatorType())
  {
  }

  Solver(const MatrixType& matrix, const CommunicatorType& communicator)
    : matrix_(matrix)
    , communicator_(communicator)
  {
  }

  static std::vector<std::string> options()
  {
    return
    {
      "bicgstab.amg.ilu0"
#if !HAVE_MPI
          ,
          "bicgstab.ilut"
#endif
    };
  } // ... options()

  static Common::Configuration options(const std::string& type)
  {
    SolverUtils::check_given(type, options());
    Common::Configuration iterative_options({"max_iter", "precision", "verbose"}, {"10000", "1e-10", "0"});
    iterative_options.set("post_check_solves_system", "1e-5");
    if (type == "bicgstab.amg.ilu0") {
      iterative_options.set("smoother.iterations", "1");
      iterative_options.set("smoother.relaxation_factor", "1");
      iterative_options.set("smoother.verbose", "0");
      iterative_options.set("preconditioner.max_level", "100");
      iterative_options.set("preconditioner.coarse_target", "1000");
      iterative_options.set("preconditioner.min_coarse_rate", "1.2");
      iterative_options.set("preconditioner.prolong_damp", "1.6");
      iterative_options.set("preconditioner.anisotropy_dim", "2"); // <- this should be the dimDomain of the problem!
      iterative_options.set("preconditioner.isotropy_dim", "2"); // <- this as well
      iterative_options.set("preconditioner.verbose", "0");
#if !HAVE_MPI
    } else if (type == "bicgstab.ilut") {
      iterative_options.set("preconditioner.iterations", "2");
      iterative_options.set("preconditioner.relaxation_factor", "1.0");
#endif // !HAVE_MPI
    } else
      DUNE_THROW(Exceptions::internal_error,
                 "Given type '" << type << "' is not supported, although it was reported by options()!");
    iterative_options.set("type", type);
    return iterative_options;
  } // ... options(...)

  void apply(const IstlDenseVector<S>& rhs, IstlDenseVector<S>& solution) const
  {
    apply(rhs, solution, options()[0]);
  }

  void apply(const IstlDenseVector<S>& rhs, IstlDenseVector<S>& solution, const std::string& type) const
  {
    apply(rhs, solution, options(type));
  }

  /**
   *  \note does a copy of the rhs
   */
  void apply(const IstlDenseVector<S>& rhs, IstlDenseVector<S>& solution, const Common::Configuration& opts) const
  {
    try {
      if (!opts.has_key("type"))
        DUNE_THROW(Exceptions::configuration_error,
                   "Given options (see below) need to have at least the key 'type' set!\n\n" << opts);
      const auto type = opts.get<std::string>("type");
      SolverUtils::check_given(type, options());
      const Common::Configuration default_opts = options(type);
      IstlDenseVector<S> writable_rhs          = rhs.copy();
      // solve
      if (type == "bicgstab.amg.ilu0") {
        auto result = AmgApplicator<S, CommunicatorType>(matrix_, communicator_.storage_access())
                          .call(writable_rhs, solution, opts, default_opts);
        if (!result.converged)
          DUNE_THROW(Exceptions::linear_solver_failed_bc_it_did_not_converge,
                     "The dune-istl backend reported 'InverseOperatorResult.converged == false'!\n"
                         << "Those were the given options:\n\n"
                         << opts);
#if !HAVE_MPI
      } else if (type == "bicgstab.ilut") {
        typedef MatrixAdapter<typename MatrixType::BackendType,
                              typename IstlDenseVector<S>::BackendType,
                              typename IstlDenseVector<S>::BackendType> MatrixOperatorType;
        MatrixOperatorType matrix_operator(matrix_.backend());
        typedef SeqILUn<typename MatrixType::BackendType,
                        typename IstlDenseVector<S>::BackendType,
                        typename IstlDenseVector<S>::BackendType> PreconditionerType;
        PreconditionerType preconditioner(
            matrix_.backend(),
            opts.get("preconditioner.iterations", default_opts.get<int>("preconditioner.iterations")),
            opts.get("preconditioner.relaxation_factor", default_opts.get<S>("preconditioner.relaxation_factor")));
        typedef BiCGSTABSolver<typename IstlDenseVector<S>::BackendType> SolverType;
        SolverType solver(matrix_operator,
                          preconditioner,
                          opts.get("precision", default_opts.get<S>("precision")),
                          opts.get("max_iter", default_opts.get<int>("max_iter")),
                          opts.get("verbose", default_opts.get<int>("verbose")));
        InverseOperatorResult stat;
        solver.apply(solution.backend(), writable_rhs.backend(), stat);
        if (!stat.converged)
          DUNE_THROW(Exceptions::linear_solver_failed_bc_it_did_not_converge,
                     "The dune-istl backend reported 'InverseOperatorResult.converged == false'!\n"
                         << "Those were the given options:\n\n"
                         << opts);
#endif // !HAVE_MPI
      } else
        DUNE_THROW(Exceptions::internal_error,
                   "Given type '" << type << "' is not supported, although it was reported by options()!");
      // check (use writable_rhs as tmp)
      const S post_check_solves_system_threshold =
          opts.get("post_check_solves_system", default_opts.get<S>("post_check_solves_system"));
      if (post_check_solves_system_threshold > 0) {
        matrix_.mv(solution, writable_rhs);
        writable_rhs -= rhs;
        const S sup_norm = writable_rhs.sup_norm();
        if (sup_norm > post_check_solves_system_threshold || std::isnan(sup_norm) || std::isinf(sup_norm))
          DUNE_THROW(Exceptions::linear_solver_failed_bc_the_solution_does_not_solve_the_system,
                     "The computed solution does not solve the system (although the dune-istl backend "
                         << "reported no error) and you requested checking (see options below)!\n"
                         << "If you want to disable this check, set 'post_check_solves_system = 0' in the options."
                         << "\n\n"
                         << "  (A * x - b).sup_norm() = "
                         << writable_rhs.sup_norm()
                         << "\n\n"
                         << "Those were the given options:\n\n"
                         << opts);
      }
    } catch (ISTLError& e) {
      DUNE_THROW(Exceptions::linear_solver_failed, "The dune-istl backend reported: " << e.what());
    }
  } // ... apply(...)

private:
  const MatrixType& matrix_;
  const Common::ConstStorageProvider<CommunicatorType> communicator_;
}; // class Solver


#else // HAVE_DUNE_ISTL


template <class S, class CommunicatorType>
class Solver<IstlRowMajorSparseMatrix<S>, CommunicatorType>
{
  static_assert(Dune::AlwaysFalse<S>::value, "You are missing dune-istl!");
};


#endif // HAVE_DUNE_ISTL


} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_SOLVER_ISTL_HH
