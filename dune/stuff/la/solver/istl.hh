// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
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

struct FoobarSequentialInformation : public Dune::Amg::SequentialInformation
{

  template <typename T>
  T globalSum(const T& t) const
  {
    return t;
  }

  template <class T1, class T2>
  void dot(const T1& x, const T1& y, T2& result) const
  {
    result = 0;
    for (typename T1::size_type i = 0; i < x.size(); i++)
      result += x[i] * y[i];
  }

  template <class T1>
  double norm(const T1& x) const
  {
    double result = 0;
    for (typename T1::size_type i = 0; i < x.size(); i++)
      result += x[i].two_norm2();
  }
};

template <class S, class CommunicatorType>
class Solver<IstlRowMajorSparseMatrix<S>, CommunicatorType> : protected SolverUtils
{
public:
  typedef IstlRowMajorSparseMatrix<S> MatrixType;

#if !HAVE_MPI
  Solver(const MatrixType& matrix)
    : matrix_(matrix)
  {
  }
#endif // !HAVE_MPI

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
    try {
      if (!opts.has_key("type"))
        DUNE_THROW_COLORFULLY(Exceptions::configuration_error,
                              "Given options (see below) need to have at least the key 'type' set!\n\n" << opts);
      const auto type = opts.get<std::string>("type");
      SolverUtils::check_given(type, options());
      const Common::ConfigTree default_opts = options(type);
      IstlDenseVector<S> writable_rhs       = rhs.copy();
      // solve
      if (type == "bicgstab.amg.ilu0") {
        // define the matrix operator
        typedef typename MatrixType::BackendType IstlMatrixType;
        typedef typename IstlDenseVector<S>::BackendType IstlVectorType;
#if HAVE_MPI
        typedef OverlappingSchwarzOperator<IstlMatrixType, IstlVectorType, IstlVectorType, CommunicatorType>
            MatrixOperatorType;
        MatrixOperatorType matrix_operator(matrix_.backend(), communicator_);
#else // HAVE_MPI
        typedef MatrixAdapter<IstlMatrixType, IstlVectorType, IstlVectorType> MatrixOperatorType;
        MatrixOperatorType matrix_operator(matrix_.backend());
#endif // HAVE_MPI

// define the scalar product
#if HAVE_MPI
        OverlappingSchwarzScalarProduct<IstlVectorType, CommunicatorType> scalar_product(communicator_);
#else
        Dune::SeqScalarProduct<typename IstlDenseVector<S>::BackendType> scalar_product;
#endif

// define ILU0 as the smoother for the AMG
#if HAVE_MPI
        typedef SeqILU0<IstlMatrixType, IstlVectorType, IstlVectorType, 1> SequentialSmootherType;
        typedef BlockPreconditioner<IstlVectorType, IstlVectorType, CommunicatorType, SequentialSmootherType>
            SmootherType;
#else // HAVE_MPI
        typedef SeqILU0<IstlMatrixType, IstlVectorType, IstlVectorType> SmootherType;
#endif
        typename Amg::SmootherTraits<SmootherType>::Arguments smoother_parameters;
        smoother_parameters.iterations =
            opts.get("smoother.iterations", default_opts.get<size_t>("smoother.iterations"));
        smoother_parameters.relaxationFactor =
            opts.get("smoother.relaxation_factor", default_opts.get<S>("smoother.relaxation_factor"));

        // define the AMG as the preconditioner for the BiCGStab solver
        Amg::Parameters amg_parameters(
            opts.get("preconditioner.max_level", default_opts.get<size_t>("preconditioner.max_level")),
            opts.get("preconditioner.coarse_target", default_opts.get<size_t>("preconditioner.coarse_target")),
            opts.get("preconditioner.min_coarse_rate", default_opts.get<S>("preconditioner.min_coarse_rate")),
            opts.get("preconditioner.prolong_damp", default_opts.get<S>("preconditioner.prolong_damp")));
        amg_parameters.setDefaultValuesIsotropic(
            opts.get("preconditioner.isotropy_dim", default_opts.get<size_t>("preconditioner.isotropy_dim")));
        amg_parameters.setDefaultValuesAnisotropic(
            opts.get("preconditioner.anisotropy_dim", default_opts.get<size_t>("preconditioner.anisotropy_dim")));
        amg_parameters.setDebugLevel(
            opts.get("preconditioner.verbose", default_opts.get<size_t>("preconditioner.verbose")));
        Amg::CoarsenCriterion<Amg::SymmetricCriterion<IstlMatrixType, Amg::FirstDiagonal>> amg_criterion(
            amg_parameters);
#if HAVE_MPI
        typedef Amg::AMG<MatrixOperatorType, IstlVectorType, SmootherType, CommunicatorType> PreconditionerType;
        PreconditionerType preconditioner(matrix_operator, amg_criterion, smoother_parameters, communicator_);
#else // HAVE_MPI
        typedef Amg::AMG<MatrixOperatorType, IstlVectorType, SmootherType> PreconditionerType;
        PreconditionerType preconditioner(matrix_operator, amg_criterion, smoother_parameters);
#endif // HAVE_MPI

        // define the BiCGStab as the actual solver
        BiCGSTABSolver<IstlVectorType> solver(
            matrix_operator,
            scalar_product,
            preconditioner,
            opts.get("precision", default_opts.get<S>("precision")),
            opts.get("max_iter", default_opts.get<size_t>("max_iter")),
#if HAVE_MPI
            (communicator_.communicator().rank() == 0) ? opts.get("verbose", default_opts.get<int>("verbose")) : 0
#else // HAVE_MPI
            opts.get("verbose", default_opts.get<int>("verbose"))
#endif
            );

        // call the solver an do the actual work
        InverseOperatorResult stats;
        solver.apply(solution.backend(), writable_rhs.backend(), stats);
        if (!stats.converged)
          DUNE_THROW_COLORFULLY(Exceptions::linear_solver_failed_bc_it_did_not_converge,
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
    } catch (ISTLError& e) {
      DUNE_THROW(Exceptions::linear_solver_failed, "The dune-istl backend reported: " << e.what());
    }
  } // ... apply(...)

private:
  const MatrixType& matrix_;
#if HAVE_MPI
  const CommunicatorType& communicator_;
#endif
}; // class Solver

template <class S>
class Solver<IstlRowMajorSparseMatrix<S>, FoobarSequentialInformation> : protected SolverUtils
{
public:
  typedef FoobarSequentialInformation CommunicatorType;
  typedef IstlRowMajorSparseMatrix<S> MatrixType;


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
    try {
      if (!opts.has_key("type"))
        DUNE_THROW_COLORFULLY(Exceptions::configuration_error,
                              "Given options (see below) need to have at least the key 'type' set!\n\n" << opts);
      const auto type = opts.get<std::string>("type");
      SolverUtils::check_given(type, options());
      const Common::ConfigTree default_opts = options(type);
      IstlDenseVector<S> writable_rhs       = rhs.copy();
      // solve
      if (type == "bicgstab.amg.ilu0") {
        // define the matrix operator
        typedef typename MatrixType::BackendType IstlMatrixType;
        typedef typename IstlDenseVector<S>::BackendType IstlVectorType;

        typedef MatrixAdapter<IstlMatrixType, IstlVectorType, IstlVectorType> MatrixOperatorType;
        MatrixOperatorType matrix_operator(matrix_.backend());

        // define the scalar product
        Dune::SeqScalarProduct<typename IstlDenseVector<S>::BackendType> scalar_product;

        // define ILU0 as the smoother for the AMG

        typedef SeqILU0<IstlMatrixType, IstlVectorType, IstlVectorType> SmootherType;

        typename Amg::SmootherTraits<SmootherType>::Arguments smoother_parameters;
        smoother_parameters.iterations =
            opts.get("smoother.iterations", default_opts.get<size_t>("smoother.iterations"));
        smoother_parameters.relaxationFactor =
            opts.get("smoother.relaxation_factor", default_opts.get<S>("smoother.relaxation_factor"));

        // define the AMG as the preconditioner for the BiCGStab solver
        Amg::Parameters amg_parameters(
            opts.get("preconditioner.max_level", default_opts.get<size_t>("preconditioner.max_level")),
            opts.get("preconditioner.coarse_target", default_opts.get<size_t>("preconditioner.coarse_target")),
            opts.get("preconditioner.min_coarse_rate", default_opts.get<S>("preconditioner.min_coarse_rate")),
            opts.get("preconditioner.prolong_damp", default_opts.get<S>("preconditioner.prolong_damp")));
        amg_parameters.setDefaultValuesIsotropic(
            opts.get("preconditioner.isotropy_dim", default_opts.get<size_t>("preconditioner.isotropy_dim")));
        amg_parameters.setDefaultValuesAnisotropic(
            opts.get("preconditioner.anisotropy_dim", default_opts.get<size_t>("preconditioner.anisotropy_dim")));
        amg_parameters.setDebugLevel(
            opts.get("preconditioner.verbose", default_opts.get<size_t>("preconditioner.verbose")));
        Amg::CoarsenCriterion<Amg::SymmetricCriterion<IstlMatrixType, Amg::FirstDiagonal>> amg_criterion(
            amg_parameters);

        typedef Amg::AMG<MatrixOperatorType, IstlVectorType, SmootherType> PreconditionerType;
        PreconditionerType preconditioner(matrix_operator, amg_criterion, smoother_parameters);

        // define the BiCGStab as the actual solver
        BiCGSTABSolver<IstlVectorType> solver(matrix_operator,
                                              scalar_product,
                                              preconditioner,
                                              opts.get("precision", default_opts.get<S>("precision")),
                                              opts.get("max_iter", default_opts.get<size_t>("max_iter")),
                                              opts.get("verbose", default_opts.get<int>("verbose")));

        // call the solver an do the actual work
        InverseOperatorResult stats;
        solver.apply(solution.backend(), writable_rhs.backend(), stats);
        if (!stats.converged)
          DUNE_THROW_COLORFULLY(Exceptions::linear_solver_failed_bc_it_did_not_converge,
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
    } catch (ISTLError& e) {
      DUNE_THROW(Exceptions::linear_solver_failed, "The dune-istl backend reported: " << e.what());
    }
  } // ... apply(...)

private:
  const MatrixType& matrix_;
  const CommunicatorType& communicator_;
};

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

#endif // DUNE_STUFF_LA_SOLVER_ISTL_HH
