// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner   (2019)

#ifndef DUNE_XT_LA_SOLVER_ISTL_SCHURCOMPLEMENT_HH
#define DUNE_XT_LA_SOLVER_ISTL_SCHURCOMPLEMENT_HH

#if HAVE_DUNE_ISTL
#  include <dune/istl/operators.hh>
#  include <dune/istl/solvers.hh>
#endif // HAVE_DUNE_ISTL

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/la/container/istl.hh>
#include <dune/xt/la/solver.hh>

namespace Dune {
namespace XT {
namespace LA {


#if HAVE_DUNE_ISTL


// For a saddle point matrix (A B1; B2^T C) this models the Schur complement (B2^T A^{-1} B1 - C)
template <class VectorType = IstlDenseVector<double>,
          class MatrixType = IstlRowMajorSparseMatrix<double>,
          class CommunicatorType = SequentialCommunication>
class SchurComplementOperator : public Dune::LinearOperator<VectorType, VectorType>
{
  using BaseType = Dune::LinearOperator<VectorType, VectorType>;

public:
  using Vector = VectorType;
  using Matrix = MatrixType;
  using Field = typename VectorType::ScalarType;
  using SolverType = Solver<Matrix, CommunicatorType>;

  // Matrix dimensions are
  // A: m x m, B1, B2: m x n, C: n x n
  SchurComplementOperator(const Matrix& _A,
                          const Matrix& _B1,
                          const Matrix& _B2,
                          const Matrix& _C,
                          const Common::Configuration& solver_opts = SolverOptions<Matrix>::options())
    : A_(_A)
    , A_inv_(make_solver(A_))
    , B1_(_B1)
    , B2_(_B2)
    , C_(_C)
    , solver_opts_(solver_opts)
    , m_vec_1_(_A.rows())
    , m_vec_2_(_A.rows())
    , n_vec_1_(_C.rows())
    , n_vec_2_(_C.rows())
  {}

  SchurComplementOperator(const SchurComplementOperator& other)
    : A_(other.A_)
    , A_inv_(make_solver(A_))
    , B1_(other.B1_)
    , B2_(other.B2_)
    , C_(other.C_)
    , solver_opts_(other.solver_opts_)
    , m_vec_1_(other.m_vec_1_)
    , m_vec_2_(other.m_vec_2_)
    , n_vec_1_(other.n_vec_1_)
    , n_vec_2_(other.n_vec_2_)
  {}

  /*! \brief apply operator to x:  \f$ y = S(x) \f$
        The input vector is consistent and the output must also be
     consistent on the interior+border partition.
   */
  virtual void apply(const Vector& x, Vector& y) const override final
  {
    // we want to calculate y = (B2^T A^{-1} B1 - C) x
    // calculate B1 x
    auto& B1x = m_vec_1_;
    B1_.mv(x, B1x);
    // calculate A^{-1} B1 x
    auto& AinvB1x = m_vec_2_;
    A_inv_.apply(B1x, AinvB1x, solver_opts_);
    // apply B2^T
    B2_.mtv(AinvB1x, y);
    // calculate Cx
    auto& Cx = n_vec_1_;
    C_.mv(x, Cx);
    y -= Cx;
  }

  virtual void applyscaleadd(Field alpha, const Vector& x, Vector& y) const override final
  {
    auto Sx = n_vec_2_;
    apply(x, Sx);
    Sx *= alpha;
    y += Sx;
  }

  //! Category of the linear operator (see SolverCategory::Category)
  virtual SolverCategory::Category category() const override final
  {
    return SolverCategory::Category::sequential;
  }

  const SolverType& A_inv() const
  {
    return A_inv_;
  }

  const Matrix& A() const
  {
    return A_;
  }

  const Matrix& B1() const
  {
    return B1_;
  }

  const Matrix& B2() const
  {
    return B2_;
  }

  const Matrix& C() const
  {
    return C_;
  }

private:
  const Matrix& A_;
  const SolverType A_inv_;
  const Matrix& B1_;
  const Matrix& B2_;
  const Matrix& C_;
  const Common::Configuration solver_opts_;
  // vectors to store intermediate results
  mutable Vector m_vec_1_;
  mutable Vector m_vec_2_;
  mutable Vector n_vec_1_;
  mutable Vector n_vec_2_;
};

#else // HAVE_DUNE_ISTL

template <class VectorType = IstlDenseVector<double>,
          class MatrixType = IstlRowMajorSparseMatrix<double>,
          class CommunicatorType = SequentialCommunication>
class SchurComplementOperator
{
  static_assert(Dune::AlwaysFalse<VectorType>::value, "You are missing dune-istl!");
};

#endif // HAVE_DUNE_ISTL

} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_SOLVER_ISTL_SCHURCOMPLEMENT_HH
