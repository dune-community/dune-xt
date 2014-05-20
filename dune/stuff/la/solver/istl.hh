// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_SOLVER_STUFF_HH
#define DUNE_STUFF_LA_SOLVER_STUFF_HH

#include "config.h"

#include <type_traits>
#include <cmath>

#if HAVE_DUNE_ISTL
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/paamg/amg.hh>
#endif // HAVE_DUNE_ISTL

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configtree.hh>
#include <dune/stuff/la/container/istl.hh>

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
  {
  }

  Solver(const MatrixType& matrix, const CommunicatorType&
#if HAVE_MPI
                                       communicator
#else
/*communicator*/
#endif
         )
    : matrix_(matrix)
#if HAVE_MPI
    , communicator_(communicator)
#endif
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

  static Common::ConfigTree options(const std::string& type)
  {
    SolverUtils::check_given(type, options());
    Common::ConfigTree iterative_options({"max_iter", "precision", "verbose"}, {"10000", "1e-10", "0"});
#if !HAVE_MPI
    iterative_options.set("post_check_solves_system", "1e-5");
#endif
    if (type == "bicgstab.amg.ilu0") {
      iterative_options.set("smoother.iterations", "1");
      iterative_options.set("smoother.relaxation_factor", "1");
      iterative_options.set("smoother.max_level", "15");
      iterative_options.set("smoother.coarse_target", "2000");
      iterative_options.set("smoother.min_coarse_rate", "1.2");
      iterative_options.set("smoother.prolong_damp", "1.6");
      iterative_options.set("smoother.anisotropy_dim", "2");
      iterative_options.set("smoother.verbose", "0");
#if !HAVE_MPI
    } else if (type == "bicgstab.ilut") {
      iterative_options.set("preconditioner.iterations", "2");
      iterative_options.set("preconditioner.relaxation_factor", "1.0");
#endif // !HAVE_MPI
    } else
      DUNE_THROW_COLORFULLY(Exceptions::internal_error,
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
  void apply(const IstlDenseVector<S>& rhs, IstlDenseVector<S>& solution, const Common::ConfigTree& opts) const
  {
    if (!opts.has_key("type"))
      DUNE_THROW_COLORFULLY(Exceptions::configuration_error,
                            "Given options (see below) need to have at least the key 'type' set!\n\n" << opts);
    const auto type = opts.get<std::string>("type");
    SolverUtils::check_given(type, options());
    const Common::ConfigTree default_opts = options(type);
    IstlDenseVector<S> writable_rhs       = rhs.copy();
    // solve
    if (type == "bicgstab.amg.ilu0") {
#if HAVE_MPI
      typedef typename MatrixType::BackendType IstlMatrixType;
      typedef typename IstlDenseVector<S>::BackendType IstlVectorType;

      typedef SeqILU0<IstlMatrixType, IstlVectorType, IstlVectorType, 1> SequentialPreconditionerType;
      typedef BlockPreconditioner<IstlVectorType, IstlVectorType, CommunicatorType, SequentialPreconditionerType>
          PreconditionerType;

      typedef OverlappingSchwarzOperator<IstlMatrixType, IstlVectorType, IstlVectorType, CommunicatorType>
          MatrixOperatorType;
      MatrixOperatorType matrix_operator(matrix_.backend(), communicator_);

      typedef Amg::AMG<MatrixOperatorType, IstlVectorType, PreconditionerType, CommunicatorType> SmootherType;
      typename Amg::SmootherTraits<PreconditionerType>::Arguments smoother_arguments;
      smoother_arguments.iterations = opts.get("smoother.iterations", default_opts.get<size_t>("smoother.iterations"));
      smoother_arguments.relaxationFactor =
          opts.get("smoother.relaxation_factor", default_opts.get<S>("smoother.relaxation_factor"));
      Amg::Parameters smoother_parameters(
          opts.get("smoother.max_level", default_opts.get<size_t>("smoother.max_level")),
          opts.get("smoother.coarse_target", default_opts.get<size_t>("smoother.coarse_target")),
          opts.get("smoother.min_coarse_rate", default_opts.get<S>("smoother.min_coarse_rate")),
          opts.get("smoother.prolong_damp", default_opts.get<S>("smoother.prolong_damp")));
      smoother_parameters.setDefaultValuesAnisotropic(
          opts.get("smoother.anisotropy_dim",
                   default_opts.get<size_t>("smoother.anisotropy_dim"))); // <- dim
      Amg::CoarsenCriterion<Amg::SymmetricCriterion<IstlMatrixType, Amg::FirstDiagonal>> smoother_criterion(
          smoother_parameters);
      SmootherType smoother(matrix_operator, smoother_criterion, smoother_arguments, communicator_);

      OverlappingSchwarzScalarProduct<IstlVectorType, CommunicatorType> scalar_product(communicator_);

      InverseOperatorResult statistics;
      BiCGSTABSolver<IstlVectorType> solver(
          matrix_operator,
          scalar_product,
          smoother,
          opts.get("precision", default_opts.get<S>("precision")),
          opts.get("max_iter", default_opts.get<size_t>("max_iter")),
          (communicator_.communicator().rank() == 0) ? opts.get("verbose", default_opts.get<int>("verbose")) : 0);
      solver.apply(solution.backend(), writable_rhs.backend(), statistics);
      if (!statistics.converged)
        DUNE_THROW_COLORFULLY(Exceptions::linear_solver_failed_bc_it_did_not_converge,
                              "The dune-istl backend reported 'InverseOperatorResult.converged == false'!\n"
                                  << "Those were the given options:\n\n"
                                  << opts);
#else // HAVE_MPI
      typedef MatrixAdapter<typename MatrixType::BackendType,
                            typename IstlDenseVector<S>::BackendType,
                            typename IstlDenseVector<S>::BackendType> MatrixOperatorType;
      MatrixOperatorType matrix_operator(matrix_.backend());
      typedef SeqILU0<typename MatrixType::BackendType,
                      typename IstlDenseVector<S>::BackendType,
                      typename IstlDenseVector<S>::BackendType> PreconditionerType;
      typedef Dune::Amg::AMG<MatrixOperatorType, typename IstlDenseVector<S>::BackendType, PreconditionerType>
          SmootherType;
      Dune::SeqScalarProduct<typename IstlDenseVector<S>::BackendType> scalar_product;
      typedef typename Dune::Amg::SmootherTraits<PreconditionerType>::Arguments SmootherArgs;
      SmootherArgs smootherArgs;
      smootherArgs.iterations = opts.get("smoother.iterations", default_opts.get<size_t>("smoother.iterations"));
      smootherArgs.relaxationFactor =
          opts.get("smoother.relaxation_factor", default_opts.get<S>("smoother.relaxation_factor"));
      Dune::Amg::Parameters params(
          opts.get("smoother.max_level", default_opts.get<size_t>("smoother.max_level")),
          opts.get("smoother.coarse_target", default_opts.get<size_t>("smoother.coarse_target")),
          opts.get("smoother.min_coarse_rate", default_opts.get<S>("smoother.min_coarse_rate")),
          opts.get("smoother.prolong_damp", default_opts.get<S>("smoother.prolong_damp")));
      params.setDefaultValuesAnisotropic(opts.get("smoother.anisotropy_dim",
                                                  default_opts.get<size_t>("smoother.anisotropy_dim"))); // <- dim
      typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<typename MatrixType::BackendType,
                                                                        Dune::Amg::FirstDiagonal>> AmgCriterion;
      AmgCriterion amg_criterion(params);
      amg_criterion.setDebugLevel(opts.get("smoother.verbose", default_opts.get<size_t>("smoother.verbose")));
      SmootherType smoother(matrix_operator, amg_criterion, smootherArgs);
      typedef BiCGSTABSolver<typename IstlDenseVector<S>::BackendType> SolverType;
      SolverType solver(matrix_operator,
                        scalar_product,
                        smoother,
                        opts.get("precision", default_opts.get<S>("precision")),
                        opts.get("max_iter", default_opts.get<size_t>("max_iter")),
                        opts.get("verbose", default_opts.get<size_t>("verbose")));
      InverseOperatorResult stat;
      solver.apply(solution.backend(), writable_rhs.backend(), stat);
      if (!stat.converged)
        DUNE_THROW_COLORFULLY(Exceptions::linear_solver_failed_bc_it_did_not_converge,
                              "The dune-istl backend reported 'InverseOperatorResult.converged == false'!\n"
                                  << "Those were the given options:\n\n"
                                  << opts);
#endif // HAVE_MPI
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
          opts.get("preconditioner.iterations", default_opts.get<size_t>("preconditioner.iterations")),
          opts.get("preconditioner.relaxation_factor", default_opts.get<S>("preconditioner.relaxation_factor")));
      typedef BiCGSTABSolver<typename IstlDenseVector<S>::BackendType> SolverType;
      SolverType solver(matrix_operator,
                        preconditioner,
                        opts.get("precision", default_opts.get<S>("precision")),
                        opts.get("max_iter", default_opts.get<size_t>("max_iter")),
                        opts.get("verbose", default_opts.get<int>("verbose")));
      InverseOperatorResult stat;
      solver.apply(solution.backend(), writable_rhs.backend(), stat);
      if (!stat.converged)
        DUNE_THROW_COLORFULLY(Exceptions::linear_solver_failed_bc_it_did_not_converge,
                              "The dune-istl backend reported 'InverseOperatorResult.converged == false'!\n"
                                  << "Those were the given options:\n\n"
                                  << opts);
#endif // !HAVE_MPI
    } else
      DUNE_THROW_COLORFULLY(Exceptions::internal_error,
                            "Given type '" << type << "' is not supported, although it was reported by options()!");
#if !HAVE_MPI
    // check (use writable_rhs as tmp)
    const S post_check_solves_system_theshhold =
        opts.get("post_check_solves_system", default_opts.get<S>("post_check_solves_system"));
    if (post_check_solves_system_theshhold > 0) {
      matrix_.mv(solution, writable_rhs);
      writable_rhs -= rhs;
      const S sup_norm = writable_rhs.sup_norm();
      if (sup_norm > post_check_solves_system_theshhold || std::isnan(sup_norm) || std::isinf(sup_norm))
        DUNE_THROW_COLORFULLY(
            Exceptions::linear_solver_failed_bc_the_solution_does_not_solve_the_system,
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
#endif // !HAVE_MPI
  } // ... apply(...)

private:
  const MatrixType& matrix_;
#if HAVE_MPI
  const CommunicatorType& communicator_;
#endif
}; // class Solver


#else // HAVE_DUNE_ISTL

template <class S>
class Solver<IstlRowMajorSparseMatrix<S>>
{
  static_assert(Dune::AlwaysFalse<S>::value, "You are missing dune-istl!");
};

#endif // HAVE_DUNE_ISTL

} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_SOLVER_STUFF_HH
