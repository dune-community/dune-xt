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
#include <dune/istl/umfpack.hh>
#include <dune/istl/superlu.hh>
#endif // HAVE_DUNE_ISTL

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/la/container/istl.hh>
#include <dune/stuff/la/solver/istl_amg.hh>

#include <dune/common/version.hh>

#include "../solver.hh"

namespace Dune {
namespace Stuff {
namespace LA {

#if HAVE_DUNE_ISTL

/**
 * \not
 **/
template <class S, class CommunicatorType>
struct IstlSolverTraits
{
  typedef typename IstlDenseVector<S>::BackendType IstlVectorType;
  typedef typename IstlRowMajorSparseMatrix<S>::BackendType IstlMatrixType;
  typedef OverlappingSchwarzOperator<IstlMatrixType, IstlVectorType, IstlVectorType, CommunicatorType>
      MatrixOperatorType;
  typedef OverlappingSchwarzScalarProduct<IstlVectorType, CommunicatorType> ScalarproductType;

  static MatrixOperatorType make_operator(const IstlMatrixType& matrix, const CommunicatorType& communicator)
  {
    return MatrixOperatorType(matrix, communicator);
  }

  static ScalarproductType make_scalarproduct(const CommunicatorType& communicator)
  {
    return ScalarproductType(communicator);
  }

  template <class SequentialPreconditionerType>
  static BlockPreconditioner<IstlVectorType, IstlVectorType, CommunicatorType, SequentialPreconditionerType>
  make_preconditioner(SequentialPreconditionerType& seq_preconditioner, const CommunicatorType& communicator)
  {
    return BlockPreconditioner<IstlVectorType, IstlVectorType, CommunicatorType, SequentialPreconditionerType>(
        seq_preconditioner, communicator);
  }
};

template <class S>
struct IstlSolverTraits<S, SequentialCommunication>
{
  typedef typename IstlDenseVector<S>::BackendType IstlVectorType;
  typedef typename IstlRowMajorSparseMatrix<S>::BackendType IstlMatrixType;
  typedef MatrixAdapter<IstlMatrixType, IstlVectorType, IstlVectorType> MatrixOperatorType;
  typedef SeqScalarProduct<IstlVectorType> ScalarproductType;

  static MatrixOperatorType make_operator(const IstlMatrixType& matrix, const SequentialCommunication& /*communicator*/)
  {
    return MatrixOperatorType(matrix);
  }

  static ScalarproductType make_scalarproduct(const SequentialCommunication& /*communicator*/)
  {
    return ScalarproductType();
  }

  template <class SequentialPreconditionerType>
  static SequentialPreconditionerType make_preconditioner(const SequentialPreconditionerType& seq_preconditioner,
                                                          const SequentialCommunication& /*communicator*/)
  {
    return seq_preconditioner;
  }
};

template <class S, class CommunicatorType>
class Solver<IstlRowMajorSparseMatrix<S>, CommunicatorType> : protected SolverUtils
{
public:
  typedef IstlRowMajorSparseMatrix<S> MatrixType;
  typedef typename MatrixType::RealType R;

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

  static std::vector<std::string> types()
  {
    return
    {
#if !HAVE_MPI && HAVE_SUPERLU
      "superlu",
#endif
          "bicgstab.amg.ssor", "bicgstab.amg.ilu0", "bicgstab.ilut", "bicgstab.ssor", "bicgstab"
#if HAVE_UMFPACK
          ,
          "umfpack"
#endif
    };
  } // ... types()

  static Common::Configuration options(const std::string type = "")
  {
    const std::string tp = !type.empty() ? type : types()[0];
    SolverUtils::check_given(tp, types());
    Common::Configuration general_opts({"type", "post_check_solves_system", "verbose"}, {tp, "1e-5", "0"});
    Common::Configuration iterative_options({"max_iter", "precision"}, {"10000", "1e-10"});
    iterative_options += general_opts;
    if (tp.substr(0, 13) == "bicgstab.amg." || tp == "bicgstab") {
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
      return iterative_options;
    } else if (tp == "bicgstab.ilut" || tp == "bicgstab.ssor") {
      iterative_options.set("preconditioner.iterations", "2");
      iterative_options.set("preconditioner.relaxation_factor", "1.0");
      return iterative_options;
#if HAVE_UMFPACK
    } else if (tp == "umfpack") {
      return general_opts;
#endif
#if !HAVE_MPI && HAVE_SUPERLU
    } else if (tp == "superlu") {
      return general_opts;
#endif // !HAVE_MPI && HAVE_SUPERLU
    } else
      DUNE_THROW(Exceptions::internal_error, "Given solver type '" << tp << "' has no default options");
    return Common::Configuration();
  } // ... options(...)

  void apply(const IstlDenseVector<S>& rhs, IstlDenseVector<S>& solution) const
  {
    apply(rhs, solution, types()[0]);
  }

  void apply(const IstlDenseVector<S>& rhs, IstlDenseVector<S>& solution, const std::string& type) const
  {
    apply(rhs, solution, options(type));
  }

  int verbosity(const Common::Configuration& opts, const Common::Configuration& default_opts) const
  {
    const auto actual_value = opts.get("verbose", default_opts.get<int>("verbose"));
    return
#if HAVE_MPI
        (communicator_.storage_access().communicator().rank() == 0) ? actual_value : 0;
#else
        actual_value;
#endif
  }

  /**
   *  \note does a copy of the rhs
   */
  void apply(const IstlDenseVector<S>& rhs, IstlDenseVector<S>& solution, const Common::Configuration& opts) const
  {
    typedef IstlSolverTraits<S, CommunicatorType> Traits;
    typedef typename Traits::IstlVectorType IstlVectorType;
    typedef typename Traits::MatrixOperatorType MatrixOperatorType;
    typedef BiCGSTABSolver<IstlVectorType> BiCgSolverType;
    InverseOperatorResult solver_result;
    auto scalar_product = Traits::make_scalarproduct(communicator_.storage_access());

    try {
      if (!opts.has_key("type"))
        DUNE_THROW(Exceptions::configuration_error,
                   "Given options (see below) need to have at least the key 'type' set!\n\n" << opts);
      const auto type = opts.get<std::string>("type");
      SolverUtils::check_given(type, types());
      const Common::Configuration default_opts = options(type);
      IstlDenseVector<S> writable_rhs          = rhs.copy();

      if (type.substr(0, 13) == "bicgstab.amg.") {
        solver_result = AmgApplicator<S, CommunicatorType>(matrix_, communicator_.storage_access())
                            .call(writable_rhs, solution, opts, default_opts, type.substr(13));
      } else if (type == "bicgstab.ilut") {
        auto matrix_operator = Traits::make_operator(matrix_.backend(), communicator_.storage_access());
        typedef SeqILUn<typename MatrixType::BackendType, IstlVectorType, IstlVectorType> SequentialPreconditionerType;
        SequentialPreconditionerType seq_preconditioner(
            matrix_.backend(),
            opts.get("preconditioner.iterations", default_opts.get<int>("preconditioner.iterations")),
            opts.get("preconditioner.relaxation_factor", default_opts.get<S>("preconditioner.relaxation_factor")));
        auto preconditioner = Traits::make_preconditioner(seq_preconditioner, communicator_.storage_access());
        BiCgSolverType solver(matrix_operator,
                              scalar_product,
                              preconditioner,
                              opts.get("precision", default_opts.get<R>("precision")),
                              opts.get("max_iter", default_opts.get<int>("max_iter")),
                              verbosity(opts, default_opts));
        solver.apply(solution.backend(), writable_rhs.backend(), solver_result);
      } else if (type == "bicgstab.ssor") {
        auto matrix_operator = Traits::make_operator(matrix_.backend(), communicator_.storage_access());
        typedef SeqSSOR<typename MatrixType::BackendType, IstlVectorType, IstlVectorType> SequentialPreconditionerType;
        SequentialPreconditionerType seq_preconditioner(
            matrix_.backend(),
            opts.get("preconditioner.iterations", default_opts.get<int>("preconditioner.iterations")),
            opts.get("preconditioner.relaxation_factor", default_opts.get<S>("preconditioner.relaxation_factor")));
        auto preconditioner = Traits::make_preconditioner(seq_preconditioner, communicator_.storage_access());
        BiCgSolverType solver(matrix_operator,
                              scalar_product,
                              preconditioner,
                              opts.get("precision", default_opts.get<S>("precision")),
                              opts.get("max_iter", default_opts.get<int>("max_iter")),
                              verbosity(opts, default_opts));
        solver.apply(solution.backend(), writable_rhs.backend(), solver_result);
      } else if (type == "bicgstab") {
        auto matrix_operator = Traits::make_operator(matrix_.backend(), communicator_.storage_access());
        constexpr auto cat   = decltype(matrix_operator)::category;
        typedef IdentityPreconditioner<MatrixOperatorType, cat> SequentialPreconditioner;
        SequentialPreconditioner seq_preconditioner;
        auto preconditioner = Traits::make_preconditioner(seq_preconditioner, communicator_.storage_access());
        // define the BiCGStab as the actual solver
        BiCgSolverType solver(matrix_operator,
                              scalar_product,
                              preconditioner,
                              opts.get("precision", default_opts.get<S>("precision")),
                              opts.get("max_iter", default_opts.get<int>("max_iter")),
                              verbosity(opts, default_opts));
        solver.apply(solution.backend(), writable_rhs.backend(), solver_result);
#if HAVE_UMFPACK
      } else if (type == "umfpack") {
        UMFPack<typename MatrixType::BackendType> solver(matrix_.backend(),
                                                         opts.get("verbose", default_opts.get<int>("verbose")));
        solver.apply(solution.backend(), writable_rhs.backend(), solver_result);
#endif // HAVE_UMFPACK
#if !HAVE_MPI && HAVE_SUPERLU
      } else if (type == "superlu") {
        SuperLU<typename MatrixType::BackendType> solver(matrix_.backend(),
                                                         opts.get("verbose", default_opts.get<int>("verbose")));

        solver.apply(solution.backend(), writable_rhs.backend(), solver_result);
#endif // !HAVE_MPI && HAVE_SUPERLU
      } else
        DUNE_THROW(Exceptions::internal_error,
                   "Given type '" << type << "' is not supported, although it was reported by types()!");
      if (!solver_result.converged)
        DUNE_THROW(Exceptions::linear_solver_failed_bc_it_did_not_converge,
                   "The dune-istl backend reported 'InverseOperatorResult.converged == false'!\n"
                       << "Those were the given options:\n\n"
                       << opts);

      // check (use writable_rhs as tmp)
      const R post_check_solves_system_threshold =
          opts.get("post_check_solves_system", default_opts.get<R>("post_check_solves_system"));
      if (post_check_solves_system_threshold > 0) {
        matrix_.mv(solution, writable_rhs);
        communicator_.storage_access().copyOwnerToAll(writable_rhs.backend(), writable_rhs.backend());
        writable_rhs -= rhs;
        const R sup_norm = writable_rhs.sup_norm();
        if (sup_norm > post_check_solves_system_threshold || DSC::isnan(sup_norm) || DSC::isinf(sup_norm))
          DUNE_THROW(Exceptions::linear_solver_failed_bc_the_solution_does_not_solve_the_system,
                     "The computed solution does not solve the system (although the dune-istl backend "
                         << "reported no error) and you requested checking (see options below)!\n"
                         << "If you want to disable this check, set 'post_check_solves_system = 0' in the options."
                         << "\n\n"
                         << "  (A * x - b).sup_norm() = "
                         << sup_norm
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
