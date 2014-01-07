// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright Holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_SOLVER_STUFF_HH
#define DUNE_STUFF_LA_SOLVER_STUFF_HH

#include <type_traits>

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


template <class S>
class Solver<IstlRowMajorSparseMatrix<S>> : protected SolverUtils
{
public:
  typedef IstlRowMajorSparseMatrix<S> MatrixType;

  Solver(const MatrixType& matrix)
    : matrix_(matrix)
  {
  }

  static std::vector<std::string> options()
  {
    return {"bicgstab.ilut", "bicgstab.amg.ilu0"};
  }

  static Common::ConfigTree options(const std::string& type)
  {
    SolverUtils::check_given(type, options());
    Common::ConfigTree iterative_options({"max_iter", "precision", "verbose"}, {"10000", "1e-10", "0"});
    if (type == "bicgstab.ilut") {
      iterative_options.add("preconditioner.iterations", "10");
      iterative_options.add("preconditioner.relaxation_factor", "1.0");
    } else if (type == "bicgstab.amg.ilu0") {
      iterative_options.add("smoother.iterations", "1");
      iterative_options.add("smoother.relaxation_factor", "1");
      iterative_options.add("smoother.max_level", "15");
      iterative_options.add("smoother.coarse_target", "2000");
      iterative_options.add("smoother.min_coarse_rate", "1.2");
      iterative_options.add("smoother.prolong_damp", "1.6");
      iterative_options.add("smoother.anisotropy_dim", "2");
      iterative_options.add("smoother.verbose", "0");
    } else
      DUNE_THROW_COLORFULLY(Exception::internal_error,
                            "Given type '" << type << "' is not supported, although it was reported by options()!");
    iterative_options.add("type", type);
    return iterative_options;
  } // ... options(...)

  size_t apply(const IstlDenseVector<S>& rhs, IstlDenseVector<S>& solution) const
  {
    return apply(rhs, solution, options()[0]);
  }

  size_t apply(const IstlDenseVector<S>& rhs, IstlDenseVector<S>& solution, const std::string& type) const
  {
    return apply(rhs, solution, options(type));
  } // ... apply(...)

  /**
   *  \note does a copy of the rhs
   */
  size_t apply(const IstlDenseVector<S>& rhs, IstlDenseVector<S>& solution, const Common::ConfigTree& opts) const
  {
    if (!opts.has_key("type"))
      DUNE_THROW_COLORFULLY(Exception::configuration_error,
                            "Given options (see below) need to have at least the key 'type' set!\n\n" << opts);
    const auto type = opts.get<std::string>("type");
    this->check_given(type, options());
    const Common::ConfigTree default_opts = options(type);
    IstlDenseVector<S> writable_rhs       = rhs.copy();
    // solve
    if (type == "bicgstab.ilut") {
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
        return 2;
    } else if (type == "bicgstab.amg.ilu0") {
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
        return 2;
    } else
      DUNE_THROW_COLORFULLY(Exception::internal_error,
                            "Given type '" << type << "' is not supported, although it was reported by options()!");
// check
#ifndef NDEBUG
    matrix_.mv(solution, writable_rhs);
    if (!rhs.almost_equal(writable_rhs))
      return 4;
#endif // NDEBUG
    return 0;
  } // ... apply(...)

private:
  const MatrixType& matrix_;
}; // class Solver


#else // HAVE_DUNE_ISTL

template <class S>
class Solver<IstlRowMajorSparseMatrix<S>>
{
  static_assert(Dune::AlwaysFalse<ScalarImp>::value, "You are missing dune-istl!");
};

#endif // HAVE_DUNE_ISTL

} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_SOLVER_STUFF_HH
