// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_SOLVER_ISTL_AMG_HH
#define DUNE_STUFF_LA_SOLVER_ISTL_AMG_HH

#include <type_traits>
#include <cmath>

#if HAVE_DUNE_ISTL
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/stuff/common/disable_warnings.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/stuff/common/reenable_warnings.hh>
#endif // HAVE_DUNE_ISTL

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configcontainer.hh>
#include <dune/stuff/common/parallel/helper.hh>
#include <dune/stuff/la/container/istl.hh>

namespace Dune {
namespace Stuff {
namespace LA {

#if HAVE_DUNE_ISTL


//! the general, parallel case
template <class S, class CommunicatorType>
class AmgApplicator
{
  typedef IstlRowMajorSparseMatrix<S> MatrixType;
  typedef typename MatrixType::BackendType IstlMatrixType;
  typedef typename IstlDenseVector<S>::BackendType IstlVectorType;

public:
  AmgApplicator(const MatrixType& matrix, const CommunicatorType& comm)
    : matrix_(matrix)
    , communicator_(comm)
  {
  }

  InverseOperatorResult call(IstlDenseVector<S>& rhs, IstlDenseVector<S>& solution, const Common::ConfigContainer& opts,
                             const Common::ConfigContainer& default_opts)
  {
    // define the matrix operator

    typedef OverlappingSchwarzOperator<IstlMatrixType, IstlVectorType, IstlVectorType, CommunicatorType>
        MatrixOperatorType;
    MatrixOperatorType matrix_operator(matrix_.backend(), communicator_);

    // define the scalar product
    OverlappingSchwarzScalarProduct<IstlVectorType, CommunicatorType> scalar_product(communicator_);

    // define ILU0 as the smoother for the AMG
    typedef SeqILU0<IstlMatrixType, IstlVectorType, IstlVectorType, 1> SequentialSmootherType;
    typedef BlockPreconditioner<IstlVectorType, IstlVectorType, CommunicatorType, SequentialSmootherType> SmootherType;
    typename Amg::SmootherTraits<SmootherType>::Arguments smoother_parameters;
    smoother_parameters.iterations = opts.get("smoother.iterations", default_opts.get<size_t>("smoother.iterations"));
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
    amg_parameters.setDebugLevel(opts.get("preconditioner.verbose", default_opts.get<int>("preconditioner.verbose")));
    Amg::CoarsenCriterion<Amg::SymmetricCriterion<IstlMatrixType, Amg::FirstDiagonal>> amg_criterion(amg_parameters);
    typedef Amg::AMG<MatrixOperatorType, IstlVectorType, SmootherType, CommunicatorType> PreconditionerType;
    PreconditionerType preconditioner(matrix_operator, amg_criterion, smoother_parameters, communicator_);

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

    InverseOperatorResult stats;
    solver.apply(solution.backend(), rhs.backend(), stats);
    return stats;
  } // ... call(...)
protected:
  const MatrixType& matrix_;
  const CommunicatorType& communicator_;
};


//! specialization for our faux type \ref SequentialCommunication
template <class S>
class AmgApplicator<S, SequentialCommunication>
{
  typedef IstlRowMajorSparseMatrix<S> MatrixType;
  typedef typename MatrixType::BackendType IstlMatrixType;
  typedef typename IstlDenseVector<S>::BackendType IstlVectorType;

public:
  AmgApplicator(const MatrixType& matrix, const SequentialCommunication& comm)
    : matrix_(matrix)
    , communicator_(comm)
  {
  }

  InverseOperatorResult call(IstlDenseVector<S>& rhs, IstlDenseVector<S>& solution, const Common::ConfigContainer& opts,
                             const Common::ConfigContainer& default_opts)
  {
    typedef MatrixAdapter<IstlMatrixType, IstlVectorType, IstlVectorType> MatrixOperatorType;
    MatrixOperatorType matrix_operator(matrix_.backend());

    // define the scalar product
    Dune::SeqScalarProduct<typename IstlDenseVector<S>::BackendType> scalar_product;

    // define ILU0 as the smoother for the AMG

    typedef SeqILU0<IstlMatrixType, IstlVectorType, IstlVectorType> SmootherType;

    typename Amg::SmootherTraits<SmootherType>::Arguments smoother_parameters;
    smoother_parameters.iterations = opts.get("smoother.iterations", default_opts.get<int>("smoother.iterations"));
    smoother_parameters.relaxationFactor =
        opts.get("smoother.relaxation_factor", default_opts.get<S>("smoother.relaxation_factor"));

    // define the AMG as the preconditioner for the BiCGStab solver
    Amg::Parameters amg_parameters(
        opts.get("preconditioner.max_level", default_opts.get<int>("preconditioner.max_level")),
        opts.get("preconditioner.coarse_target", default_opts.get<int>("preconditioner.coarse_target")),
        opts.get("preconditioner.min_coarse_rate", default_opts.get<S>("preconditioner.min_coarse_rate")),
        opts.get("preconditioner.prolong_damp", default_opts.get<S>("preconditioner.prolong_damp")));
    amg_parameters.setDefaultValuesIsotropic(
        opts.get("preconditioner.isotropy_dim", default_opts.get<size_t>("preconditioner.isotropy_dim")));
    amg_parameters.setDefaultValuesAnisotropic(
        opts.get("preconditioner.anisotropy_dim", default_opts.get<size_t>("preconditioner.anisotropy_dim")));
    amg_parameters.setDebugLevel(opts.get("preconditioner.verbose", default_opts.get<int>("preconditioner.verbose")));
    Amg::CoarsenCriterion<Amg::SymmetricCriterion<IstlMatrixType, Amg::FirstDiagonal>> amg_criterion(amg_parameters);

    typedef Amg::AMG<MatrixOperatorType, IstlVectorType, SmootherType> PreconditionerType;
    PreconditionerType preconditioner(matrix_operator, amg_criterion, smoother_parameters);

    // define the BiCGStab as the actual solver
    BiCGSTABSolver<IstlVectorType> solver(matrix_operator,
                                          scalar_product,
                                          preconditioner,
                                          opts.get("precision", default_opts.get<S>("precision")),
                                          opts.get("max_iter", default_opts.get<int>("max_iter")),
                                          opts.get("verbose", default_opts.get<int>("verbose")));
    InverseOperatorResult stats;
    solver.apply(solution.backend(), rhs.backend(), stats);
    return stats;
  } // ... call(...)

protected:
  const MatrixType& matrix_;
  const SequentialCommunication& communicator_;
};


#else // HAVE_DUNE_ISTL


template <class S, class T>
class AmgApplicator
{
  static_assert(Dune::AlwaysFalse<S, T>::value, "You are missing dune-istl!");
};


#endif // HAVE_DUNE_ISTL

} // namespace LA
} // namespace Stuff
} // namespace Dune


#endif // DUNE_STUFF_LA_SOLVER_ISTL_AMG_HH
