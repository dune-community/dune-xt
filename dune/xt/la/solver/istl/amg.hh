// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2014 - 2017)
//   Rene Milk        (2014 - 2016, 2018)
//   Tobias Leibner   (2014, 2017)

#ifndef DUNE_XT_LA_SOLVER_ISTL_AMG_HH
#define DUNE_XT_LA_SOLVER_ISTL_AMG_HH

#include <type_traits>
#include <cmath>

#if HAVE_DUNE_ISTL
#  include <dune/istl/operators.hh>
#  include <dune/istl/solvers.hh>
#  include <dune/istl/paamg/amg.hh>
#  include <dune/istl/preconditioners.hh>
#endif // HAVE_DUNE_ISTL

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/parallel/helper.hh>
#include <dune/xt/la/container/istl.hh>

namespace Dune {
namespace XT {
namespace LA {

#if HAVE_DUNE_ISTL

template <class O>
class IdentityPreconditioner : public Dune::Preconditioner<typename O::domain_type, typename O::range_type>
{
public:
  //! \brief The domain type of the preconditioner.
  typedef typename O::domain_type domain_type;
  //! \brief The range type of the preconditioner.
  typedef typename O::range_type range_type;
  //! \brief The field type of the preconditioner.
  typedef typename range_type::field_type field_type;
  typedef O InverseOperator;

  IdentityPreconditioner(const SolverCategory::Category cat)
    : category_(cat)
  {
  }

  //! Category of the preconditioner (see SolverCategory::Category)
  virtual SolverCategory::Category category() const override final
  {
    return category_;
  }

  virtual void pre(domain_type&, range_type&) override final
  {
  }

  virtual void apply(domain_type& v, const range_type& d) override final
  {
    v = d;
  }

  virtual void post(domain_type&) override final
  {
  }

private:
  SolverCategory::Category category_;
};

//! the general, parallel case
template <class S, class CommunicatorType>
class AmgApplicator
{
  typedef IstlRowMajorSparseMatrix<S> MatrixType;
  typedef typename MatrixType::RealType R;
  typedef typename MatrixType::BackendType IstlMatrixType;
  typedef typename IstlDenseVector<S>::BackendType IstlVectorType;

public:
  AmgApplicator(const MatrixType& matrix, const CommunicatorType& comm)
    : matrix_(matrix)
    , communicator_(comm)
  {}

  InverseOperatorResult call(IstlDenseVector<S>& rhs,
                             IstlDenseVector<S>& solution,
                             const Common::Configuration& opts,
                             const Common::Configuration& default_opts,
                             const std::string& smoother_type)
  {
    // define the matrix operator

    typedef OverlappingSchwarzOperator<IstlMatrixType, IstlVectorType, IstlVectorType, CommunicatorType>
        MatrixOperatorType;
    MatrixOperatorType matrix_operator(matrix_.backend(), communicator_);

    // define the scalar product
    OverlappingSchwarzScalarProduct<IstlVectorType, CommunicatorType> scalar_product(communicator_);

    // ILU0 as the smoother for the AMG
    typedef SeqILU0<IstlMatrixType, IstlVectorType, IstlVectorType, 1> SequentialSmootherType_ILU;
    typedef BlockPreconditioner<IstlVectorType, IstlVectorType, CommunicatorType, SequentialSmootherType_ILU>
        SmootherType_ILU;
    typename Amg::SmootherTraits<SmootherType_ILU>::Arguments smoother_parameters_ILU;
    smoother_parameters_ILU.iterations =
        opts.get("smoother.iterations", default_opts.get<size_t>("smoother.iterations"));
    smoother_parameters_ILU.relaxationFactor =
        opts.get("smoother.relaxation_factor", default_opts.get<S>("smoother.relaxation_factor"));
    // SSOR as the smoother for the amg
    typedef SeqSSOR<IstlMatrixType, IstlVectorType, IstlVectorType, 1> SequentialSmootherType_SSOR;
    typedef BlockPreconditioner<IstlVectorType, IstlVectorType, CommunicatorType, SequentialSmootherType_SSOR>
        SmootherType_SSOR;
    typename Amg::SmootherTraits<SmootherType_SSOR>::Arguments smoother_parameters_SSOR;
    smoother_parameters_SSOR.iterations =
        opts.get("smoother.iterations", default_opts.get<size_t>("smoother.iterations"));
    smoother_parameters_SSOR.relaxationFactor =
        opts.get("smoother.relaxation_factor", default_opts.get<S>("smoother.relaxation_factor"));

    // define the AMG as the preconditioner for the BiCGStab solver
    Amg::Parameters amg_parameters(
        opts.get("preconditioner.max_level", default_opts.get<size_t>("preconditioner.max_level")),
        opts.get("preconditioner.coarse_target", default_opts.get<size_t>("preconditioner.coarse_target")),
        opts.get("preconditioner.min_coarse_rate", default_opts.get<R>("preconditioner.min_coarse_rate")),
        opts.get("preconditioner.prolong_damp", default_opts.get<R>("preconditioner.prolong_damp")));
    amg_parameters.setDefaultValuesIsotropic(
        opts.get("preconditioner.isotropy_dim", default_opts.get<size_t>("preconditioner.isotropy_dim")));
    amg_parameters.setDefaultValuesAnisotropic(
        opts.get("preconditioner.anisotropy_dim", default_opts.get<size_t>("preconditioner.anisotropy_dim")));
    amg_parameters.setDebugLevel(opts.get("preconditioner.verbose", default_opts.get<int>("preconditioner.verbose")));
    Amg::CoarsenCriterion<Amg::UnSymmetricCriterion<IstlMatrixType, Amg::FirstDiagonal>> amg_criterion(amg_parameters);
    if (smoother_type == "ilu0") {
      typedef Amg::AMG<MatrixOperatorType, IstlVectorType, SmootherType_ILU, CommunicatorType> PreconditionerType_ILU;
      PreconditionerType_ILU preconditioner(matrix_operator, amg_criterion, smoother_parameters_ILU, communicator_);

      // define the BiCGStab as the actual solver
      BiCGSTABSolver<IstlVectorType> solver(
          matrix_operator,
          scalar_product,
          preconditioner,
          opts.get("precision", default_opts.get<S>("precision")),
          opts.get("max_iter", default_opts.get<size_t>("max_iter")),
#  if HAVE_MPI
          (communicator_.communicator().rank() == 0) ? opts.get("verbose", default_opts.get<int>("verbose")) : 0
#  else // HAVE_MPI
          opts.get("verbose", default_opts.get<int>("verbose"))
#  endif
      );

      InverseOperatorResult stats;
      solver.apply(solution.backend(), rhs.backend(), stats);
      return stats;
    } else if (smoother_type == "ssor") {
      typedef Amg::AMG<MatrixOperatorType, IstlVectorType, SmootherType_SSOR, CommunicatorType> PreconditionerType_SSOR;
      PreconditionerType_SSOR preconditioner(matrix_operator, amg_criterion, smoother_parameters_ILU, communicator_);

      // define the BiCGStab as the actual solver
      BiCGSTABSolver<IstlVectorType> solver(
          matrix_operator,
          scalar_product,
          preconditioner,
          opts.get("precision", default_opts.get<S>("precision")),
          opts.get("max_iter", default_opts.get<size_t>("max_iter")),
#  if HAVE_MPI
          (communicator_.communicator().rank() == 0) ? opts.get("verbose", default_opts.get<int>("verbose")) : 0
#  else // HAVE_MPI
          opts.get("verbose", default_opts.get<int>("verbose"))
#  endif
      );

      InverseOperatorResult stats;
      solver.apply(solution.backend(), rhs.backend(), stats);
      return stats;
    } else
      DUNE_THROW(Common::Exceptions::wrong_input_given, "Unknown smoother requested: " << smoother_type);
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
  typedef typename MatrixType::RealType R;
  typedef typename MatrixType::BackendType IstlMatrixType;
  typedef typename IstlDenseVector<S>::BackendType IstlVectorType;

public:
  AmgApplicator(const MatrixType& matrix, const SequentialCommunication& comm)
    : matrix_(matrix)
    , communicator_(comm)
  {}

  InverseOperatorResult call(IstlDenseVector<S>& rhs,
                             IstlDenseVector<S>& solution,
                             const Common::Configuration& opts,
                             const Common::Configuration& default_opts,
                             const std::string& smoother_type)
  {
    typedef MatrixAdapter<IstlMatrixType, IstlVectorType, IstlVectorType> MatrixOperatorType;
    MatrixOperatorType matrix_operator(matrix_.backend());

    // define the scalar product
    Dune::SeqScalarProduct<typename IstlDenseVector<S>::BackendType> scalar_product;

    // define the AMG as the preconditioner for the BiCGStab solver
    Amg::Parameters amg_parameters(
        opts.get("preconditioner.max_level", default_opts.get<int>("preconditioner.max_level")),
        opts.get("preconditioner.coarse_target", default_opts.get<int>("preconditioner.coarse_target")),
        opts.get("preconditioner.min_coarse_rate", default_opts.get<R>("preconditioner.min_coarse_rate")),
        opts.get("preconditioner.prolong_damp", default_opts.get<R>("preconditioner.prolong_damp")));
    amg_parameters.setDefaultValuesIsotropic(
        opts.get("preconditioner.isotropy_dim", default_opts.get<size_t>("preconditioner.isotropy_dim")));
    amg_parameters.setDefaultValuesAnisotropic(
        opts.get("preconditioner.anisotropy_dim", default_opts.get<size_t>("preconditioner.anisotropy_dim")));
    amg_parameters.setDebugLevel(opts.get("preconditioner.verbose", default_opts.get<int>("preconditioner.verbose")));
    Amg::CoarsenCriterion<Amg::UnSymmetricCriterion<IstlMatrixType, Amg::FirstDiagonal>> amg_criterion(amg_parameters);

    InverseOperatorResult stats;
    if (smoother_type == "ilu0") {
      typedef SeqILU0<IstlMatrixType, IstlVectorType, IstlVectorType> SmootherType;

      typename Amg::SmootherTraits<SmootherType>::Arguments smoother_parameters;
      smoother_parameters.iterations = opts.get("smoother.iterations", default_opts.get<int>("smoother.iterations"));
      smoother_parameters.relaxationFactor =
          opts.get("smoother.relaxation_factor", default_opts.get<S>("smoother.relaxation_factor"));
      typedef Amg::AMG<MatrixOperatorType, IstlVectorType, SmootherType> PreconditionerType;
      PreconditionerType preconditioner(matrix_operator, amg_criterion, smoother_parameters);
      // define the BiCGStab as the actual solver
      BiCGSTABSolver<IstlVectorType> solver(matrix_operator,
                                            scalar_product,
                                            preconditioner,
                                            opts.get("precision", default_opts.get<S>("precision")),
                                            opts.get("max_iter", default_opts.get<int>("max_iter")),
                                            opts.get("verbose", default_opts.get<int>("verbose")));
      solver.apply(solution.backend(), rhs.backend(), stats);
    } else if (smoother_type == "ssor") {
      typedef SeqSSOR<IstlMatrixType, IstlVectorType, IstlVectorType> SmootherType;

      typename Amg::SmootherTraits<SmootherType>::Arguments smoother_parameters;
      smoother_parameters.iterations = opts.get("smoother.iterations", default_opts.get<int>("smoother.iterations"));
      smoother_parameters.relaxationFactor =
          opts.get("smoother.relaxation_factor", default_opts.get<S>("smoother.relaxation_factor"));
      typedef Amg::AMG<MatrixOperatorType, IstlVectorType, SmootherType> PreconditionerType;
      PreconditionerType preconditioner(matrix_operator, amg_criterion, smoother_parameters);
      BiCGSTABSolver<IstlVectorType> solver(matrix_operator,
                                            scalar_product,
                                            preconditioner,
                                            opts.get("precision", default_opts.get<S>("precision")),
                                            opts.get("max_iter", default_opts.get<int>("max_iter")),
                                            opts.get("verbose", default_opts.get<int>("verbose")));
      solver.apply(solution.backend(), rhs.backend(), stats);
    } else {
      DUNE_THROW(Common::Exceptions::wrong_input_given, "Unknown smoother requested: " << smoother_type);
    }
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
  static_assert(Dune::AlwaysFalse<S>::value, "You are missing dune-istl!");
};

#endif // HAVE_DUNE_ISTL

} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_SOLVER_ISTL_AMG_HH
