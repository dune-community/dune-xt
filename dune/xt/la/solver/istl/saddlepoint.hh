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

#ifndef DUNE_XT_LA_SOLVER_ISTL_SADDLEPOINT_HH
#define DUNE_XT_LA_SOLVER_ISTL_SADDLEPOINT_HH

#include "config.h"

#include <type_traits>
#include <cmath>

#if HAVE_DUNE_ISTL
#  include <dune/istl/operators.hh>
#  include <dune/istl/solvers.hh>
#endif // HAVE_DUNE_ISTL

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/la/container/istl.hh>
#include <dune/xt/la/solver.hh>

#include "preconditioners.hh"
#include "schurcomplement.hh"

namespace Dune {
namespace XT {
namespace LA {


#if HAVE_DUNE_ISTL


// Solver for saddle point system (A B1; B2^T C) (u; p) = (f; g) using the Schur complement, i.e., solve (B2^T A^{-1} B1
// - C) p = B2^T A^{-1} f - g first and then u = A^{-1} (F - B1 p)
template <class FieldType = double, class CommunicatorType = SequentialCommunication>
class SaddlePointSolver
{
public:
  using Vector = IstlDenseVector<FieldType>;
  using Matrix = IstlRowMajorSparseMatrix<FieldType>;
  using Solver = Solver<Matrix, CommunicatorType>;

  // Matrix and vector dimensions are
  // A: m x m, B1, B2: m x n, C: n x n, f: m, g: n
  SaddlePointSolver(const Matrix& A,
                    const Matrix& B1,
                    const Matrix& B2,
                    const Matrix& C,
                    const Common::Configuration& solver_opts = SolverOptions<Matrix>::options())
    : schur_complement_(A, B1, B2, C, solver_opts)
  {}

  static std::vector<std::string> types()
  {
    std::vector<std::string> ret{"direct"};
    ret.insert(ret.begin(), "cg_direct_schurcomplement");
    return ret;
  } // ... types()

  static XT::Common::Configuration options(const std::string type = "")
  {
    const std::string tp = !type.empty() ? type : types()[0];
    internal::SolverUtils::check_given(tp, types());
    Common::Configuration general_opts({"type", "post_check_solves_system", "verbose"}, {tp.c_str(), "1e-5", "0"});
    Common::Configuration iterative_options({"max_iter", "precision"}, {"10000", "1e-10"});
    iterative_options += general_opts;
    if (tp == "direct")
      return general_opts;
    else if (tp == "cg_direct_schurcomplement")
      return iterative_options;
    else
      return general_opts;
  } // ... options(...)

  void apply(const Vector& f, const Vector& g, Vector& u, Vector& p) const
  {
    apply(f, g, u, p, types()[0]);
  }

  void apply(const Vector& f, const Vector& g, Vector& u, Vector& p, const std::string& type) const
  {
    apply(f, g, u, p, options(type));
  }

  int verbosity(const Common::Configuration& opts, const Common::Configuration& default_opts) const
  {
    const auto actual_value = opts.get("verbose", default_opts.get<int>("verbose"));
    return actual_value;
  }

  void apply(const Vector& f, const Vector& g, Vector& u, Vector& p, const Common::Configuration& opts) const
  {
    const auto type = opts.get<std::string>("type");
    if (type == "direct") {
      // copy matrices to saddle point system matrix
      // create pattern first
      const size_t m = f.size();
      const size_t n = g.size();
      XT::LA::SparsityPatternDefault system_matrix_pattern(m + n);
      const auto& A = schur_complement_.A();
      const auto& B1 = schur_complement_.B1();
      const auto& B2 = schur_complement_.B2();
      const auto& C = schur_complement_.C();
      const auto pattern_A = A.pattern();
      const auto pattern_B1 = B1.pattern();
      const auto pattern_B2 = B2.pattern();
      const auto pattern_C = C.pattern();
      for (size_t ii = 0; ii < m; ++ii) {
        for (const auto& jj : pattern_A.inner(ii))
          system_matrix_pattern.insert(ii, jj);
        for (const auto& jj : pattern_B1.inner(ii))
          system_matrix_pattern.insert(ii, m + jj);
        for (const auto& jj : pattern_B2.inner(ii))
          system_matrix_pattern.insert(m + jj, ii);
      } // ii
      for (size_t ii = 0; ii < n; ++ii)
        for (const auto& jj : pattern_C.inner(ii))
          system_matrix_pattern.insert(m + ii, m + jj);
      system_matrix_pattern.sort();

      // now copy the matrices
      Matrix system_matrix(m + n, m + n, system_matrix_pattern);
      for (size_t ii = 0; ii < m; ++ii) {
        for (const auto& jj : pattern_A.inner(ii))
          system_matrix.set_entry(ii, jj, A.get_entry(ii, jj));
        for (const auto& jj : pattern_B1.inner(ii))
          system_matrix.set_entry(ii, m + jj, B1.get_entry(ii, jj));
        for (const auto& jj : pattern_B2.inner(ii))
          system_matrix.set_entry(m + jj, ii, B2.get_entry(ii, jj));
      } // ii
      for (size_t ii = 0; ii < n; ++ii)
        for (const auto& jj : pattern_C.inner(ii))
          system_matrix.set_entry(m + ii, m + jj, C.get_entry(ii, jj));

      // also copy the rhs
      Vector system_vector(m + n, 0.), solution_vector(m + n, 0.);
      for (size_t ii = 0; ii < m; ++ii)
        system_vector[ii] = f[ii];
      for (size_t ii = 0; ii < n; ++ii)
        system_vector[m + ii] = g[ii];

      // solve the system by a direct solver
      XT::LA::solve(system_matrix, system_vector, solution_vector);

      // copy to result vectors
      for (size_t ii = 0; ii < m; ++ii)
        u[ii] = solution_vector[ii];
      for (size_t ii = 0; ii < n; ++ii)
        p[ii] = solution_vector[m + ii];
    } else if (type == "cg_direct_schurcomplement") {
      // calculate rhs B2^T A^{-1} f - g
      auto Ainv_f = f;
      auto rhs_p = g;
      schur_complement_.A_inv().apply(f, Ainv_f);
      schur_complement_.B2().mtv(Ainv_f, rhs_p);
      rhs_p -= g;

      // Solve S p = rhs
      IdentityPreconditioner<SchurComplementOperator<FieldType, CommunicatorType>> prec(
          SolverCategory::Category::sequential);
      auto schur_complement_copy = schur_complement_;
      Dune::CGSolver<typename Vector::BackendType> outer_solver(schur_complement_copy, prec, 1e-10, 10000, 0, false);
      InverseOperatorResult res;
      outer_solver.apply(p.backend(), rhs_p.backend(), res);

      // Now solve u = A^{-1}(f - B1 p)
      auto rhs_u = f;
      rhs_u -= schur_complement_.B1() * p;
      schur_complement_.A_inv().apply(rhs_u, u);
    }
  } // ... apply(...)

private:
  const SchurComplementOperator<FieldType, CommunicatorType> schur_complement_;
};

#else // HAVE_DUNE_ISTL

template <class FieldType = double, class CommunicatorType = SequentialCommunication>
class SaddlePointSolver
{
  static_assert(Dune::AlwaysFalse<FieldType>::value, "You are missing dune-istl!");
};


#endif // HAVE_DUNE_ISTL

} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_SOLVER_ISTL_SADDLEPOINT_HH
