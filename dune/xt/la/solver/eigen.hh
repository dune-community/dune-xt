// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2012 - 2017)
//   Rene Milk        (2012 - 2013, 2015 - 2016)
//   Tobias Leibner   (2014, 2017)

#ifndef DUNE_XT_LA_SOLVER_EIGEN_HH
#define DUNE_XT_LA_SOLVER_EIGEN_HH

#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <complex>

#if HAVE_EIGEN
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
//#   if HAVE_UMFPACK
//#     include <Eigen/UmfPackSupport>
//#   endif
//#   include <Eigen/SPQRSupport>
//#   include <Eigen/CholmodSupport>
//#   if HAVE_SUPERLU
//#     include <Eigen/SuperLUSupport>
//#   endif
#endif // HAVE_EIGEN

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/la/container/eigen.hh>

#include "../solver.hh"

namespace Dune {
namespace XT {
namespace LA {

#if HAVE_EIGEN


template <class S, class CommunicatorType>
class Solver<EigenDenseMatrix<S>, CommunicatorType> : protected internal::SolverUtils
{
public:
  typedef EigenDenseMatrix<S> MatrixType;
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
    return {"lu.partialpiv",
            "qr.householder",
            "llt",
            "ldlt",
            "qr.colpivhouseholder",
            "qr.fullpivhouseholder",
            "lu.fullpiv"};
  } // ... types()

  static Common::Configuration options(const std::string type = "")
  {
    const std::string tp = !type.empty() ? type : types()[0];
    internal::SolverUtils::check_given(tp, types());
    Common::Configuration default_options({"type", "post_check_solves_system", "check_for_inf_nan"}, {tp, "1e-5", "1"});
    // * for symmetric matrices
    if (tp == "ldlt" || tp == "llt") {
      default_options.set("pre_check_symmetry", "1e-8");
    }
    return default_options;
  } // ... options(...)

  template <class T1, class T2>
  void apply(const EigenBaseVector<T1, S>& rhs, EigenBaseVector<T2, S>& solution) const
  {
    apply(rhs, solution, types()[0]);
  }

  template <class T1, class T2>
  void apply(const EigenBaseVector<T1, S>& rhs, EigenBaseVector<T2, S>& solution, const std::string& type) const
  {
    apply(rhs, solution, options(type));
  }

  template <class T1, class T2>
  void
  apply(const EigenBaseVector<T1, S>& rhs, EigenBaseVector<T2, S>& solution, const Common::Configuration& opts) const
  {
    if (!opts.has_key("type"))
      DUNE_THROW(Common::Exceptions::configuration_error,
                 "Given options (see below) need to have at least the key 'type' set!\n\n"
                     << opts);
    const auto type = opts.get<std::string>("type");
    internal::SolverUtils::check_given(type, types());
    const Common::Configuration default_opts = options(type);
    // check for inf or nan
    const bool check_for_inf_nan = opts.get("check_for_inf_nan", default_opts.get<bool>("check_for_inf_nan"));
    if (check_for_inf_nan) {
      for (size_t ii = 0; ii < matrix_.rows(); ++ii) {
        for (size_t jj = 0; jj < matrix_.cols(); ++jj) {
          const S& val = matrix_.backend()(ii, jj);
          if (Common::isnan(val) || Common::isinf(val)) {
            std::stringstream msg;
            msg << "Given matrix contains inf or nan and you requested checking (see options below)!\n"
                << "If you want to disable this check, set 'check_for_inf_nan = 0' in the options.\n\n"
                << "Those were the given options:\n\n"
                << opts;
            if (rhs.size() <= internal::max_size_to_print)
              msg << "\nThis was the given matrix:\n\n" << matrix_ << "\n";
            DUNE_THROW(Exceptions::linear_solver_failed_bc_data_did_not_fulfill_requirements, msg.str());
          }
        }
      }
      for (size_t ii = 0; ii < rhs.size(); ++ii) {
        const S val = rhs.get_entry(ii);
        if (Common::isnan(val) || Common::isinf(val)) {
          std::stringstream msg;
          msg << "Given rhs contains inf or nan and you requested checking (see options below)!\n"
              << "If you want to disable this check, set 'check_for_inf_nan = 0' in the options.\n\n"
              << "Those were the given options:\n\n"
              << opts;
          if (rhs.size() <= internal::max_size_to_print)
            msg << "\nThis was the given right hand side:\n\n" << rhs << "\n";
          DUNE_THROW(Exceptions::linear_solver_failed_bc_data_did_not_fulfill_requirements, msg.str());
        }
      }
    }
    // check for symmetry (if solver needs it)
    if (type == "ldlt" || type == "llt") {
      const R pre_check_symmetry_threshhold = opts.get("pre_check_symmetry", default_opts.get<R>("pre_check_symmetry"));
      if (pre_check_symmetry_threshhold > 0) {
        const MatrixType tmp(matrix_.backend() - matrix_.backend().adjoint());
        // serialize difference to compute L^\infty error (no copy done here)
        const R error = std::max(tmp.backend().cwiseAbs().minCoeff(), tmp.backend().cwiseAbs().maxCoeff());
        if (error > pre_check_symmetry_threshhold) {
          std::stringstream msg;
          msg << "Given matrix is not symmetric and you requested checking (see options below)!\n"
              << "If you want to disable this check, set 'pre_check_symmetry = 0' in the options.\n\n"
              << "  (A - A').sup_norm() = " << error << "\n\n"
              << "Those were the given options:\n\n"
              << opts;
          if (rhs.size() <= internal::max_size_to_print)
            msg << "\nThis was the given matrix A:\n\n" << matrix_ << "\n";
          DUNE_THROW(Exceptions::linear_solver_failed_bc_data_did_not_fulfill_requirements, msg.str());
        }
      }
    }
    // solve
    if (type == "qr.colpivhouseholder") {
      solution.backend() = matrix_.backend().colPivHouseholderQr().solve(rhs.backend());
    } else if (type == "qr.fullpivhouseholder")
      solution.backend() = matrix_.backend().fullPivHouseholderQr().solve(rhs.backend());
    else if (type == "qr.householder")
      solution.backend() = matrix_.backend().householderQr().solve(rhs.backend());
    else if (type == "lu.fullpiv")
      solution.backend() = matrix_.backend().fullPivLu().solve(rhs.backend());
    else if (type == "llt")
      solution.backend() = matrix_.backend().llt().solve(rhs.backend());
    else if (type == "ldlt")
      solution.backend() = matrix_.backend().ldlt().solve(rhs.backend());
    else if (type == "lu.partialpiv")
      solution.backend() = matrix_.backend().partialPivLu().solve(rhs.backend());
    else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type << "' is not supported, although it was reported by types()!");
    // check
    if (check_for_inf_nan)
      for (size_t ii = 0; ii < solution.size(); ++ii) {
        const S val = solution.get_entry(ii);
        if (Common::isnan(val) || Common::isinf(val)) {
          std::stringstream msg;
          msg << "The computed solution contains inf or nan and you requested checking (see options "
              << "below)!\n"
              << "If you want to disable this check, set 'check_for_inf_nan = 0' in the options.\n\n"
              << "Those were the given options:\n\n"
              << opts;
          if (rhs.size() <= internal::max_size_to_print)
            msg << "\nThis was the given matrix A:\n\n"
                << matrix_ << "\nThis was the given right hand side b:\n\n"
                << rhs << "\nThis is the computed solution:\n\n"
                << solution << "\n";
          DUNE_THROW(Exceptions::linear_solver_failed_bc_data_did_not_fulfill_requirements, msg.str());
        }
      }
    const R post_check_solves_system_threshold =
        opts.get("post_check_solves_system", default_opts.get<R>("post_check_solves_system"));
    if (post_check_solves_system_threshold > 0) {
      auto tmp = rhs.copy();
      tmp.backend() = matrix_.backend() * solution.backend() - rhs.backend();
      const R sup_norm = tmp.sup_norm();
      if (sup_norm > post_check_solves_system_threshold || Common::isnan(sup_norm) || Common::isinf(sup_norm)) {
        std::stringstream msg;
        msg << "The computed solution does not solve the system (although the eigen backend reported "
            << "'Success') and you requested checking (see options below)!\n"
            << "If you want to disable this check, set 'post_check_solves_system = 0' in the options."
            << "\n\n"
            << "  (A * x - b).sup_norm() = " << tmp.sup_norm() << "\n\n"
            << "Those were the given options:\n\n"
            << opts;
        if (rhs.size() <= internal::max_size_to_print)
          msg << "\nThis was the given matrix A:\n\n"
              << matrix_ << "\nThis was the given right hand side b:\n\n"
              << rhs << "\nThis is the computed solution:\n\n"
              << solution << "\n";
        DUNE_THROW(Exceptions::linear_solver_failed_bc_the_solution_does_not_solve_the_system, msg.str());
      }
    }
  } // ... apply(...)

private:
  const MatrixType& matrix_;
}; // class Solver


/**
 *  \note lu.sparse will copy the matrix to column major
 *  \note qr.sparse will copy the matrix to column major
 *  \note ldlt.simplicial will copy the matrix to column major
 *  \note llt.simplicial will copy the matrix to column major
 */
template <class S, class CommunicatorType>
class Solver<EigenRowMajorSparseMatrix<S>, CommunicatorType> : protected internal::SolverUtils
{
  typedef ::Eigen::SparseMatrix<S, ::Eigen::ColMajor> ColMajorBackendType;

public:
  typedef EigenRowMajorSparseMatrix<S> MatrixType;
  typedef typename MatrixType::RealType R;

private:
  typedef typename MatrixType::BackendType::Index EIGEN_size_t;

public:
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
    return {
        "bicgstab.ilut",
        "lu.sparse",
        "llt.simplicial" // <- does only work with symmetric matrices
        ,
        "ldlt.simplicial" // <- does only work with symmetric matrices
        ,
        "bicgstab.diagonal" // <- slow for complicated matrices
        ,
        "bicgstab.identity" // <- slow for complicated matrices
        ,
        "qr.sparse" // <- produces correct results, but is painfully slow
        ,
        "cg.diagonal.lower" // <- does only work with symmetric matrices, may produce correct results
        ,
        "cg.diagonal.upper" // <- does only work with symmetric matrices, may produce correct results
        ,
        "cg.identity.lower" // <- does only work with symmetric matrices, may produce correct results
        ,
        "cg.identity.upper" // <- does only work with symmetric matrices, may produce correct results
        //           , "spqr"                  // <- does not compile
        //           , "llt.cholmodsupernodal" // <- does not compile
        //#if HAVE_UMFPACK
        //           , "lu.umfpack"            // <- untested
        //#endif
        //#if HAVE_SUPERLU
        //           , "superlu"               // <- untested
        //#endif
    };
  } // ... types()

  static Common::Configuration options(const std::string type = "")
  {
    const std::string tp = !type.empty() ? type : types()[0];
    // check
    internal::SolverUtils::check_given(tp, types());
    // default config
    Common::Configuration default_options({"type", "post_check_solves_system", "check_for_inf_nan"}, {tp, "1e-5", "1"});
    Common::Configuration iterative_options({"max_iter", "precision"}, {"10000", "1e-10"});
    iterative_options += default_options;
    // direct solvers
    if (tp == "lu.sparse" || tp == "qr.sparse" || tp == "lu.umfpack" || tp == "spqr" || tp == "llt.cholmodsupernodal"
        || tp == "superlu")
      return default_options;
    // * for symmetric matrices
    if (tp == "ldlt.simplicial" || tp == "llt.simplicial") {
      default_options.set("pre_check_symmetry", "1e-8");
      return default_options;
    }
    // iterative solvers
    if (tp == "bicgstab.ilut") {
      iterative_options.set("preconditioner.fill_factor", "10");
      iterative_options.set("preconditioner.drop_tol", "1e-4");
    } else if (tp.substr(0, 3) == "cg.")
      iterative_options.set("pre_check_symmetry", "1e-8");
    return iterative_options;
  } // ... options(...)

  template <class T1, class T2>
  void apply(const EigenBaseVector<T1, S>& rhs, EigenBaseVector<T2, S>& solution) const
  {
    apply(rhs, solution, types()[0]);
  }

  template <class T1, class T2>
  void apply(const EigenBaseVector<T1, S>& rhs, EigenBaseVector<T2, S>& solution, const std::string& type) const
  {
    apply(rhs, solution, options(type));
  }

  template <class T1, class T2>
  void
  apply(const EigenBaseVector<T1, S>& rhs, EigenBaseVector<T2, S>& solution, const Common::Configuration& opts) const
  {
    if (!opts.has_key("type"))
      DUNE_THROW(Common::Exceptions::configuration_error,
                 "Given options (see below) need to have at least the key 'type' set!\n\n"
                     << opts);
    const auto type = opts.get<std::string>("type");
    internal::SolverUtils::check_given(type, types());
    const Common::Configuration default_opts = options(type);
    // check for inf or nan
    const bool check_for_inf_nan = opts.get("check_for_inf_nan", default_opts.get<bool>("check_for_inf_nan"));
    if (check_for_inf_nan) {
      // iterates over the non-zero entries of matrix_.backend() and checks them
      typedef typename MatrixType::BackendType::InnerIterator InnerIterator;
      for (EIGEN_size_t ii = 0; ii < matrix_.backend().outerSize(); ++ii) {
        for (InnerIterator it(matrix_.backend(), ii); it; ++it) {
          if (Common::isnan(std::real(it.value())) || Common::isnan(std::imag(it.value()))
              || Common::isinf(std::abs(it.value())))
            DUNE_THROW(Exceptions::linear_solver_failed_bc_data_did_not_fulfill_requirements,
                       "Given matrix contains inf or nan and you requested checking (see options below)!\n"
                           << "If you want to disable this check, set 'check_for_inf_nan = 0' in the options.\n\n"
                           << "Those were the given options:\n\n"
                           << opts);
        }
      }
      for (size_t ii = 0; ii < rhs.size(); ++ii) {
        const S val = rhs.get_entry(ii);
        if (Common::isnan(val) || Common::isinf(val))
          DUNE_THROW(Exceptions::linear_solver_failed_bc_data_did_not_fulfill_requirements,
                     "Given rhs contains inf or nan and you requested checking (see options below)!\n"
                         << "If you want to disable this check, set 'check_for_inf_nan = 0' in the options.\n\n"
                         << "Those were the given options:\n\n"
                         << opts);
      }
    }
    // check for symmetry (if solver needs it)
    if (type.substr(0, 3) == "cg." || type == "ldlt.simplicial" || type == "llt.simplicial") {
      const R pre_check_symmetry_threshhold = opts.get("pre_check_symmetry", default_opts.get<R>("pre_check_symmetry"));
      if (pre_check_symmetry_threshhold > 0) {
        ColMajorBackendType colmajor_copy(matrix_.backend());
        colmajor_copy -= matrix_.backend().adjoint();
        // iterates over non-zero entries as above
        typedef typename ColMajorBackendType::InnerIterator InnerIterator;
        for (EIGEN_size_t ii = 0; ii < colmajor_copy.outerSize(); ++ii) {
          for (InnerIterator it(colmajor_copy, ii); it; ++it) {
            if (std::max(std::abs(std::real(it.value())), std::abs(std::imag(it.value())))
                > pre_check_symmetry_threshhold)
              DUNE_THROW(Exceptions::linear_solver_failed_bc_data_did_not_fulfill_requirements,
                         "Given matrix is not symmetric/hermitian and you requested checking (see options below)!\n"
                             << "If you want to disable this check, set 'pre_check_symmetry = 0' in the options.\n\n"
                             << "Those were the given options:\n\n"
                             << opts);
          }
        }
      }
    }
    ::Eigen::ComputationInfo info;
    if (type == "cg.diagonal.lower") {
      typedef ::Eigen::
          ConjugateGradient<typename MatrixType::BackendType, ::Eigen::Lower, ::Eigen::DiagonalPreconditioner<S>>
              SolverType;
      SolverType solver(matrix_.backend());
      solver.setMaxIterations(opts.get("max_iter", default_opts.get<int>("max_iter")));
      solver.setTolerance(opts.get("precision", default_opts.get<R>("precision")));
      solution.backend() = solver.solve(rhs.backend());
      info = solver.info();
    } else if (type == "cg.diagonal.upper") {
      typedef ::Eigen::
          ConjugateGradient<typename MatrixType::BackendType, ::Eigen::Upper, ::Eigen::DiagonalPreconditioner<S>>
              SolverType;
      SolverType solver(matrix_.backend());
      solver.setMaxIterations(opts.get("max_iter", default_opts.get<int>("max_iter")));
      solver.setTolerance(opts.get("precision", default_opts.get<R>("precision")));
      solution.backend() = solver.solve(rhs.backend());
      info = solver.info();
    } else if (type == "cg.identity.lower") {
      typedef ::Eigen::
          ConjugateGradient<typename MatrixType::BackendType, ::Eigen::Lower, ::Eigen::IdentityPreconditioner>
              SolverType;
      SolverType solver(matrix_.backend());
      solver.setMaxIterations(opts.get("max_iter", default_opts.get<int>("max_iter")));
      solver.setTolerance(opts.get("precision", default_opts.get<R>("precision")));
      solution.backend() = solver.solve(rhs.backend());
      info = solver.info();
    } else if (type == "cg.identity.upper") {
      typedef ::Eigen::
          ConjugateGradient<typename MatrixType::BackendType, ::Eigen::Lower, ::Eigen::IdentityPreconditioner>
              SolverType;
      SolverType solver(matrix_.backend());
      solver.setMaxIterations(opts.get("max_iter", default_opts.get<int>("max_iter")));
      solver.setTolerance(opts.get("precision", default_opts.get<R>("precision")));
      solution.backend() = solver.solve(rhs.backend());
      info = solver.info();
    } else if (type == "bicgstab.ilut") {
      typedef ::Eigen::BiCGSTAB<typename MatrixType::BackendType, ::Eigen::IncompleteLUT<S>> SolverType;
      SolverType solver(matrix_.backend());
      solver.setMaxIterations(opts.get("max_iter", default_opts.get<int>("max_iter")));
      solver.setTolerance(opts.get("precision", default_opts.get<R>("precision")));
      solver.preconditioner().setDroptol(
          opts.get("preconditioner.drop_tol", default_opts.get<R>("preconditioner.drop_tol")));
      solver.preconditioner().setFillfactor(
          opts.get("preconditioner.fill_factor", default_opts.get<int>("preconditioner.fill_factor")));
      solution.backend() = solver.solve(rhs.backend());
      info = solver.info();
    } else if (type == "bicgstab.diagonal") {
      typedef ::Eigen::BiCGSTAB<typename MatrixType::BackendType, ::Eigen::DiagonalPreconditioner<S>> SolverType;
      SolverType solver(matrix_.backend());
      solver.setMaxIterations(opts.get("max_iter", default_opts.get<int>("max_iter")));
      solver.setTolerance(opts.get("precision", default_opts.get<R>("precision")));
      solution.backend() = solver.solve(rhs.backend());
      info = solver.info();
    } else if (type == "bicgstab.identity") {
      typedef ::Eigen::BiCGSTAB<typename MatrixType::BackendType, ::Eigen::IdentityPreconditioner> SolverType;
      SolverType solver(matrix_.backend());
      solver.setMaxIterations(opts.get("max_iter", default_opts.get<int>("max_iter")));
      solver.setTolerance(opts.get("precision", default_opts.get<R>("precision")));
      solution.backend() = solver.solve(rhs.backend());
      info = solver.info();
    } else if (type == "lu.sparse") {
      ColMajorBackendType colmajor_copy(matrix_.backend());
      colmajor_copy.makeCompressed();
      typedef ::Eigen::SparseLU<ColMajorBackendType> SolverType;
      SolverType solver;
      solver.analyzePattern(colmajor_copy);
      solver.factorize(colmajor_copy);
      solution.backend() = solver.solve(rhs.backend());
      info = solver.info();
    } else if (type == "qr.sparse") {
      ColMajorBackendType colmajor_copy(matrix_.backend());
      colmajor_copy.makeCompressed();
      typedef ::Eigen::SparseQR<ColMajorBackendType, ::Eigen::COLAMDOrdering<int>> SolverType;
      SolverType solver;
      solver.analyzePattern(colmajor_copy);
      solver.factorize(colmajor_copy);
      solution.backend() = solver.solve(rhs.backend());
      info = solver.info();
    } else if (type == "ldlt.simplicial") {
      ColMajorBackendType colmajor_copy(matrix_.backend());
      colmajor_copy.makeCompressed();
      typedef ::Eigen::SimplicialLDLT<ColMajorBackendType> SolverType;
      SolverType solver;
      solver.analyzePattern(colmajor_copy);
      solver.factorize(colmajor_copy);
      solution.backend() = solver.solve(rhs.backend());
      info = solver.info();
    } else if (type == "llt.simplicial") {
      ColMajorBackendType colmajor_copy(matrix_.backend());
      colmajor_copy.makeCompressed();
      typedef ::Eigen::SimplicialLLT<ColMajorBackendType> SolverType;
      SolverType solver;
      solver.analyzePattern(colmajor_copy);
      solver.factorize(colmajor_copy);
      solution.backend() = solver.solve(rhs.backend());
      info = solver.info();
      //#if HAVE_UMFPACK
      //    } else if (type == "lu.umfpack") {
      //      typedef ::Eigen::UmfPackLU< typename MatrixType::BackendType > SolverType;
      //      SolverType solver;
      //      solver.analyzePattern(matrix_.backend());
      //      solver.factorize(matrix_.backend());
      //      solution.backend() = solver.solve(rhs.backend());
      //      info = solver.info();
      //#endif // HAVE_UMFPACK
      //    } else if (type == "spqr") {
      //      ColMajorBackendType colmajor_copy(matrix_.backend());
      //      colmajor_copy.makeCompressed();
      //      typedef ::Eigen::SPQR< ColMajorBackendType > SolverType;
      //      SolverType solver;
      //      solver.analyzePattern(colmajor_copy);
      //      solver.factorize(colmajor_copy);
      //      solution.backend() = solver.solve(rhs.backend());
      //      if (solver.info() != ::Eigen::Success)
      //        return solver.info();
      //    } else if (type == "cholmodsupernodalllt") {
      //      typedef ::Eigen::CholmodSupernodalLLT< typename MatrixType::BackendType > SolverType;
      //      SolverType solver;
      //      solver.analyzePattern(matrix_.backend());
      //      solver.factorize(matrix_.backend());
      //      solution.backend() = solver.solve(rhs.backend());
      //      if (solver.info() != ::Eigen::Success)
      //        return solver.info();
      //#if HAVE_SUPERLU
      //    } else if (type == "superlu") {
      //      typedef ::Eigen::SuperLU< typename MatrixType::BackendType > SolverType;
      //      SolverType solver;
      //      solver.analyzePattern(matrix_.backend());
      //      solver.factorize(matrix_.backend());
      //      solution.backend() = solver.solve(rhs.backend());
      //      info = solver.info();
      //#endif // HAVE_SUPERLU
    } else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type << "' is not supported, although it was reported by types()!");
    // handle eigens info
    if (info != ::Eigen::Success) {
      if (info == ::Eigen::NumericalIssue)
        DUNE_THROW(Exceptions::linear_solver_failed_bc_data_did_not_fulfill_requirements,
                   "The eigen backend reported 'NumericalIssue'!\n"
                       << "=> see http://eigen.tuxfamily.org/dox/group__enums.html#ga51bc1ac16f26ebe51eae1abb77bd037b "
                          "for eigens explanation\n"
                       << "Those were the given options:\n\n"
                       << opts);
      else if (info == ::Eigen::NoConvergence)
        DUNE_THROW(Exceptions::linear_solver_failed_bc_it_did_not_converge,
                   "The eigen backend reported 'NoConvergence'!\n"
                       << "=> see http://eigen.tuxfamily.org/dox/group__enums.html#ga51bc1ac16f26ebe51eae1abb77bd037b "
                          "for eigens explanation\n"
                       << "Those were the given options:\n\n"
                       << opts);
      else if (info == ::Eigen::InvalidInput)
        DUNE_THROW(Exceptions::linear_solver_failed_bc_it_was_not_set_up_correctly,
                   "The eigen backend reported 'InvalidInput'!\n"
                       << "=> see http://eigen.tuxfamily.org/dox/group__enums.html#ga51bc1ac16f26ebe51eae1abb77bd037b "
                          "for eigens explanation\n"
                       << "Those were the given options:\n\n"
                       << opts);
      else
        DUNE_THROW(Common::Exceptions::internal_error,
                   "The eigen backend reported an unknown status!\n"
                       << "Please report this to the dune-xt developers!");
    }
    // check
    if (check_for_inf_nan)
      for (size_t ii = 0; ii < solution.size(); ++ii) {
        const S& val = solution.get_entry(ii);
        if (Common::isnan(val) || Common::isinf(val))
          DUNE_THROW(Exceptions::linear_solver_failed_bc_data_did_not_fulfill_requirements,
                     "The computed solution contains inf or nan and you requested checking (see options "
                         << "below)!\n"
                         << "If you want to disable this check, set 'check_for_inf_nan = 0' in the options.\n\n"
                         << "Those were the given options:\n\n"
                         << opts);
      }
    const R post_check_solves_system_threshold =
        opts.get("post_check_solves_system", default_opts.get<R>("post_check_solves_system"));
    if (post_check_solves_system_threshold > 0) {
      auto tmp = rhs.copy();
      tmp.backend() = matrix_.backend() * solution.backend() - rhs.backend();
      const R sup_norm = tmp.sup_norm();
      if (sup_norm > post_check_solves_system_threshold || Common::isnan(sup_norm) || Common::isinf(sup_norm))
        DUNE_THROW(Exceptions::linear_solver_failed_bc_the_solution_does_not_solve_the_system,
                   "The computed solution does not solve the system (although the eigen backend reported "
                       << "'Success') and you requested checking (see options below)!\n"
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
}; // class Solver


#else // HAVE_EIGEN

template <class S>
class Solver<EigenDenseMatrix<S>>
{
  static_assert(Dune::AlwaysFalse<S>::value, "You are missing Eigen!");
};

template <class S>
class Solver<EigenRowMajorSparseMatrix<S>>
{
  static_assert(Dune::AlwaysFalse<S>::value, "You are missing Eigen!");
};

#endif // HAVE_EIGEN

} // namespace LA
} // namespace XT
} // namespace Dune


// begin: this is what we need for the lib
#if DUNE_XT_WITH_PYTHON_BINDINGS && HAVE_EIGEN

extern template class Dune::XT::LA::Solver<Dune::XT::LA::EigenDenseMatrix<double>>;
// extern template void
// Dune::XT::LA::Solver<Dune::XT::LA::EigenDenseMatrix<double>>::apply(const Dune::XT::LA::EigenDenseVector<double>&,
//                                                                    Dune::XT::LA::EigenDenseVector<double>&) const;
// extern template void Dune::XT::LA::Solver<Dune::XT::LA::EigenDenseMatrix<double>>::apply(
//    const Dune::XT::LA::EigenDenseVector<double>&, Dune::XT::LA::EigenDenseVector<double>&, const std::string&) const;
// extern template void
// Dune::XT::LA::Solver<Dune::XT::LA::EigenDenseMatrix<double>>::apply(const Dune::XT::LA::EigenDenseVector<double>&,
//                                                                    Dune::XT::LA::EigenDenseVector<double>&,
//                                                                    const Dune::XT::Common::Configuration&) const;

extern template class Dune::XT::LA::Solver<Dune::XT::LA::EigenRowMajorSparseMatrix<double>>;
// extern template void Dune::XT::LA::Solver<Dune::XT::LA::EigenRowMajorSparseMatrix<double>>::apply(
//    const Dune::XT::LA::EigenDenseVector<double>&, Dune::XT::LA::EigenDenseVector<double>&) const;
// extern template void Dune::XT::LA::Solver<Dune::XT::LA::EigenRowMajorSparseMatrix<double>>::apply(
//    const Dune::XT::LA::EigenDenseVector<double>&, Dune::XT::LA::EigenDenseVector<double>&, const std::string&) const;
// extern template void Dune::XT::LA::Solver<Dune::XT::LA::EigenRowMajorSparseMatrix<double>>::apply(
//    const Dune::XT::LA::EigenDenseVector<double>&,
//    Dune::XT::LA::EigenDenseVector<double>&,
//    const Dune::XT::Common::Configuration&) const;

#endif // DUNE_XT_WITH_PYTHON_BINDINGS && HAVE_EIGEN
// end: this is what we need for the lib


#endif // DUNE_XT_LA_SOLVER_EIGEN_HH
