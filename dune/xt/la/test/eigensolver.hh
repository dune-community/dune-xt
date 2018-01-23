// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_XT_LA_TEST_EIGENSOLVER_HH
#define DUNE_XT_LA_TEST_EIGENSOLVER_HH

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/unused.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/xt/la/container/conversion.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/la/eigen-solver.hh>
#include <dune/xt/la/matrix-inverter.hh>

using namespace Dune;
using namespace Dune::XT;
using namespace Dune::XT::LA;


template <class MatrixImp, class F, class C, class R>
struct EigenSolverTest : public ::testing::Test
{
  using MatrixType = MatrixImp;
  using FieldType = Common::real_t<F>;
  using ComplexMatrixType = C;
  using RealMatrixType = R;
  using RealType = Common::real_t<FieldType>;
  using EigenValuesType = std::vector<Common::complex_t<FieldType>>;
  using RealEigenValuesType = std::vector<Common::real_t<FieldType>>;
  typedef EigenSolver<MatrixType> EigenSolverType;
  typedef EigenSolverOptions<MatrixType> EigenSolverOpts;

  using M = Common::MatrixAbstraction<MatrixType>;

  EigenSolverTest()
    : all_matrices_and_expected_eigenvalues_and_vectors_are_computed_(false)
    , broken_matrix_(M::create(M::has_static_size ? M::static_rows : 1, M::has_static_size ? M::static_cols : 1))
    , unit_matrix_(
          eye_matrix<MatrixType>(M::has_static_size ? M::static_rows : 1, M::has_static_size ? M::static_cols : 1))
  {
    M::set_entry(broken_matrix_, 0, 0, std::numeric_limits<typename M::S>::infinity());
  }

  static void exports_correct_types()
  {
    const bool MatrixType_is_correct = std::is_same<typename EigenSolverType::MatrixType, MatrixType>::value;
    EXPECT_TRUE(MatrixType_is_correct);
    const bool RealType_is_correct = std::is_same<typename EigenSolverType::RealType, RealType>::value;
    EXPECT_TRUE(RealType_is_correct);
    const bool RealMatrixType_is_correct =
        std::is_same<typename EigenSolverType::RealMatrixType, RealMatrixType>::value;
    EXPECT_TRUE(RealMatrixType_is_correct);
    const bool ComplexMatrixType_is_correct =
        std::is_same<typename EigenSolverType::ComplexMatrixType, ComplexMatrixType>::value;
    EXPECT_TRUE(ComplexMatrixType_is_correct);
  } // ... exports_correct_types()

  static void has_types_and_options()
  {
    std::vector<std::string> types = EigenSolverOpts::types();
    EXPECT_GT(types.size(), 0);
    Common::Configuration opts = EigenSolverOpts::options();
    EXPECT_TRUE(opts.has_key("type"));
    EXPECT_EQ(types[0], opts.get<std::string>("type"));
    for (const auto& tp : types) {
      Common::Configuration tp_opts = EigenSolverOpts::options(tp);
      EXPECT_TRUE(tp_opts.has_key("type"));
      EXPECT_EQ(tp, tp_opts.get<std::string>("type"));
    }
  } // ... has_types_and_options(...)

  void throws_on_broken_matrix_construction() const
  {
    try {
      EigenSolverType DXTC_UNUSED(default_solver)(broken_matrix_);
      FAIL() << "Expected LA::Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements";
    } catch (const LA::Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements& /*ee*/) {
    } catch (...) {
      FAIL() << "Expected LA::Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements";
    }
    for (const auto& tp : EigenSolverOpts::types()) {
      try {
        EigenSolverType DXTC_UNUSED(default_opts_solver)(broken_matrix_, tp);
        FAIL() << "Expected LA::Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements";
      } catch (const LA::Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements& /*ee*/) {
      } catch (...) {
        FAIL() << "Expected LA::Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements";
      }
      try {
        EigenSolverType DXTC_UNUSED(solver)(broken_matrix_, EigenSolverOpts::options(tp));
        FAIL() << "Expected LA::Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements";
      } catch (const LA::Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements& /*ee*/) {
      } catch (...) {
        FAIL() << "Expected LA::Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements";
      }
    }
  } // ... throws_on_broken_matrix_construction(...)

  void allows_broken_matrix_construction_when_checks_disabled() const
  {
    for (const auto& tp : EigenSolverOpts::types()) {
      Common::Configuration opts_with_disabled_check = EigenSolverOpts::options(tp);
      opts_with_disabled_check["check_for_inf_nan"] = "false";
      EigenSolverType DXTC_UNUSED(solver)(broken_matrix_, opts_with_disabled_check);
    }
  }

  void throws_on_inconsistent_given_options() const
  {
    for (const auto& tp : EigenSolverOpts::types()) {
      Common::Configuration inconsistent_opts = EigenSolverOpts::options(tp);
      inconsistent_opts["assert_positive_eigenvalues"] = "1e-15";
      inconsistent_opts["assert_negative_eigenvalues"] = "1e-15";
      try {
        EigenSolverType DXTC_UNUSED(solver)(unit_matrix_, inconsistent_opts);
        FAIL() << "Expected LA::Exceptions::eigen_solver_failed_bc_it_was_not_set_up_correctly";
      } catch (const LA::Exceptions::eigen_solver_failed_bc_it_was_not_set_up_correctly& /*ee*/) {
      } catch (...) {
        FAIL() << "Expected LA::Exceptions::eigen_solver_failed_bc_it_was_not_set_up_correctly";
      }
    }
  } // ... throws_on_inconsistent_given_options(...)

  void is_constructible() const
  {
    ASSERT_TRUE(all_matrices_and_expected_eigenvalues_and_vectors_are_computed_);
    EigenSolverType DXTC_UNUSED(default_solver)(matrix_);
    for (const auto& tp : EigenSolverOpts::types()) {
      EigenSolverType DXTC_UNUSED(default_opts_solver)(matrix_, tp);
      EigenSolverType DXTC_UNUSED(solver)(matrix_, EigenSolverOpts::options(tp));
    }
  } // ... is_constructible(...)

  static bool find_ev(const EigenValuesType& expected_evs,
                      const XT::Common::complex_t<FieldType>& actual_ev,
                      const double& tolerance)
  {
    for (const auto& expected_ev : expected_evs) {
      if (Common::FloatCmp::eq(actual_ev, expected_ev, {tolerance, tolerance}))
        return true;
    }
    return false;
  }

  static bool find_ev(const RealEigenValuesType& expected_evs, const RealType& actual_ev, const double& tolerance)
  {
    for (const auto& expected_ev : expected_evs) {
      if (Common::FloatCmp::eq(actual_ev, expected_ev, tolerance))
        return true;
    }
    return false;
  }

  void gives_correct_eigenvalues(const Common::Configuration& tolerances = {}) const
  {
    ASSERT_TRUE(all_matrices_and_expected_eigenvalues_and_vectors_are_computed_);
    for (const auto& tp : EigenSolverOpts::types()) {
      const double tolerance = tolerances.get(tp, 1e-15);
      EigenSolverType solver(matrix_, tp);
      try {
        const auto& eigenvalues = solver.eigenvalues();
        EXPECT_EQ(Common::get_matrix_rows(matrix_), eigenvalues.size());
        for (const auto& ev : eigenvalues) {
          EXPECT_TRUE(find_ev(expected_eigenvalues_, ev, tolerance))
              << "\n\nactual eigenvalue: " << ev << "\n\nexpected eigenvalues: " << expected_eigenvalues_
              << "\n\ntype: " << tp << "\n\ntolerance: " << tolerance;
        }
      } catch (const Dune::MathError&) {
        if (tolerance > 0) {
          FAIL() << "Dune::MathError thrown when trying to get eigenvalues!"
                 << "\n\ntype: " << tp << "\n\ntolerance: " << tolerance;
        }
      }
    }
  } // ... gives_correct_eigenvalues(...)

  void gives_correct_eigenvalues_in_correct_order(const Common::Configuration& tolerances = {}) const
  {
    ASSERT_TRUE(all_matrices_and_expected_eigenvalues_and_vectors_are_computed_);
    for (const auto& tp : EigenSolverOpts::types()) {
      const double tolerance = tolerances.get(tp, 1e-15);
      EigenSolverType solver(matrix_, tp);
      try {
        const auto eigenvalues = solver.eigenvalues();
        if (tolerance > 0) {
          EXPECT_TRUE(Common::FloatCmp::eq(eigenvalues, expected_eigenvalues_, {tolerance, tolerance}))
              << "\n\nactual eigenvalues: " << eigenvalues << "\n\nexpected eigenvalues: " << expected_eigenvalues_
              << "\n\ntype: " << tp << "\n\ntolerance: " << tolerance;
        } else {
          // negative tolerance: we expect a failure
          EXPECT_FALSE(Common::FloatCmp::eq(eigenvalues, expected_eigenvalues_, {tolerance, tolerance}))
              << "\n\nTHIS IS A GOOD THING! UPDATE THE EXPECTATIONS IN tolerances!\n\n"
              << "\n\nactual eigenvalues: " << eigenvalues << "\n\nexpected eigenvalues: " << expected_eigenvalues_
              << "\n\ntype: " << tp << "\n\ntolerance: " << tolerance;
        }
      } catch (const Dune::MathError&) {
        if (tolerance > 0) {
          FAIL() << "Dune::MathError thrown when trying to get eigenvalues!"
                 << "\n\ntype: " << tp << "\n\ntolerance: " << tolerance;
        }
      }
    }
  } // ... gives_correct_eigenvalues_in_correct_order(...)

  void gives_correct_eigenvectors_in_correct_order(const Common::Configuration& tolerances = {}) const
  {
    ASSERT_TRUE(all_matrices_and_expected_eigenvalues_and_vectors_are_computed_);
    for (const auto& tp : EigenSolverOpts::types()) {
      const double tolerance = tolerances.get(tp, 1e-15);
      try {
        EigenSolverType solver(matrix_, tp);
        const auto& actual_eigenvectors = solver.eigenvectors();
        EXPECT_EQ(Common::get_matrix_rows(matrix_), Common::get_matrix_rows(actual_eigenvectors));
        EXPECT_EQ(Common::get_matrix_rows(matrix_), Common::get_matrix_cols(actual_eigenvectors));
        if (tolerance > 0) {
          EXPECT_TRUE(Common::FloatCmp::eq(
              actual_eigenvectors, expected_eigenvectors_, {tolerance, tolerance}, {tolerance, tolerance}))
              << "\n\nactual eigenvectors: " << actual_eigenvectors
              << "\n\nexpected eigenvectors: " << expected_eigenvectors_ << "\n\ntolerance: " << tolerance
              << "\n\ntype: " << tp;
        } else {
          // negative tolerance: we expect a failure
          EXPECT_FALSE(Common::FloatCmp::eq(actual_eigenvectors, expected_eigenvectors_))
              << "\n\nTHIS IS A GOOD THING! UPDATE THE EXPECTATIONS IN tolerances!\n\n"
              << "\n\nactual eigenvectors: " << actual_eigenvectors
              << "\n\nexpected eigenvectors: " << expected_eigenvectors_ << "\n\ntype: " << tp;
        }
      } catch (const Dune::MathError&) {
        if (tolerance > 0) {
          FAIL() << "Dune::MathError thrown when trying to get eigenvalues!"
                 << "\n\ntype: " << tp << "\n\ntolerance: " << tolerance;
        }
      }
    }
  } // ... gives_correct_eigenvectors_in_correct_order(...)

  void gives_correct_eigendecomposition(const Common::Configuration& tolerances = {}) const
  {
    ASSERT_TRUE(all_matrices_and_expected_eigenvalues_and_vectors_are_computed_);
    const auto matrix_as_complex = convert_to<ComplexMatrixType>(matrix_);
    for (const auto& tp : EigenSolverOpts::types()) {
      const double tolerance = tolerances.get(tp, 1e-15);
      EigenSolverType solver(matrix_, tp);
      try {
        const auto eigenvalues = solver.eigenvalues();
        const ComplexMatrixType T = solver.eigenvectors();
        ComplexMatrixType T_inv = Common::MatrixAbstraction<ComplexMatrixType>::create(Common::get_matrix_cols(T),
                                                                                       Common::get_matrix_rows(T));
        try {
          T_inv = invert_matrix(T);
        } catch (...) {
          FAIL() << "Matrix inversion has to work for this test, else add this matrix to the matrix inversion tests!";
        }
        ComplexMatrixType lambda = Common::MatrixAbstraction<ComplexMatrixType>::create(
            Common::get_matrix_rows(T), Common::get_matrix_cols(T), 0.);
        for (size_t ii = 0; ii < Common::get_matrix_rows(matrix_); ++ii)
          Common::set_matrix_entry(lambda, ii, ii, eigenvalues[ii]);
        const auto matrix_decomposition_error = (T * (lambda * T_inv)) - matrix_as_complex;
        for (size_t ii = 0; ii < Common::get_matrix_rows(matrix_); ++ii)
          for (size_t jj = 0; jj < Common::get_matrix_cols(matrix_); ++jj) {
            if (tolerance > 0)
              EXPECT_LT(std::abs(Common::get_matrix_entry(matrix_decomposition_error, ii, jj)), tolerance)
                  << "\n\ntype = " << tp << "\n\ntolerance = " << tolerance;
            else {
              // negative tolerance: we expect a failure
              EXPECT_GT(std::abs(Common::get_matrix_entry(matrix_decomposition_error, ii, jj)), 0.)
                  << "\n\nTHIS IS A GOOD THING! UPDATE THE EXPECTATIONS IN tolerances!\n\n"
                  << "\n\ntype = " << tp;
            }
          }
      } catch (const Dune::MathError&) {
        if (tolerance > 0) {
          FAIL() << "Dune::MathError thrown in eigensolver!"
                 << "\n\ntype: " << tp << "\n\ntolerance: " << tolerance;
        }
      }
    }
  } // ... gives_correct_eigendecomposition(...)

  bool all_matrices_and_expected_eigenvalues_and_vectors_are_computed_;
  MatrixType broken_matrix_;
  MatrixType unit_matrix_;
  MatrixType matrix_;
  EigenValuesType expected_eigenvalues_;
  ComplexMatrixType expected_eigenvectors_;
}; // ... struct EigenSolverTest


template <class M, class F, class C, class R>
struct EigenSolverTestForMatricesWithRealEigenvaluesAndVectors : public EigenSolverTest<M, F, C, R>
{
  using BaseType = EigenSolverTest<M, F, C, R>;
  using typename BaseType::MatrixType;
  using typename BaseType::FieldType;
  using typename BaseType::RealType;
  using typename BaseType::ComplexMatrixType;
  using typename BaseType::RealMatrixType;
  using typename BaseType::EigenValuesType;
  using typename BaseType::RealEigenValuesType;
  using typename BaseType::EigenSolverType;
  using typename BaseType::EigenSolverOpts;
  using BaseType::find_ev;

  void gives_correct_real_eigenvalues(const Common::Configuration& tolerances = {}) const
  {
    ASSERT_TRUE(all_matrices_and_expected_eigenvalues_and_vectors_are_computed_);
    for (const auto& tp : EigenSolverOpts::types()) {
      const double tolerance = tolerances.get(tp, 1e-15);
      EigenSolverType solver(matrix_, tp);
      const auto& eigenvalues = solver.eigenvalues();
      EXPECT_EQ(Common::get_matrix_rows(matrix_), eigenvalues.size());
      for (const auto& complex_ev : eigenvalues) {
        if (tolerance > 0)
          EXPECT_TRUE(Common::FloatCmp::eq(0., complex_ev.imag(), tolerance)) << "\n  type: " << tp
                                                                              << "\n  tolerance: " << tolerance;
        else {
          // negative tolerance: we expect a failure
          EXPECT_FALSE(Common::FloatCmp::eq(0., complex_ev.imag()))
              << "\n\nTHIS IS A GOOD THING! UPDATE THE EXPECTATIONS IN tolerances!\n\n"
              << "\n  type: " << tp;
        }
        const auto real_ev = complex_ev.real();
        if (tolerance > 0)
          EXPECT_TRUE(find_ev(expected_real_eigenvalues_, real_ev, tolerance))
              << "\n\nactual eigenvalue: " << real_ev << "\n\nexpected eigenvalues: " << expected_real_eigenvalues_
              << "\n\ntype: " << tp << "\n\ntolerance: " << tolerance;
        else {
          // negative tolerance: we expect a failure
          EXPECT_FALSE(find_ev(expected_real_eigenvalues_, real_ev, 1e-15))
              << "\n\nTHIS IS A GOOD THING! UPDATE THE EXPECTATIONS IN tolerances!\n\n"
              << "\n\nactual eigenvalue: " << real_ev << "\n\nexpected eigenvalues: " << expected_real_eigenvalues_
              << "\n\ntype: " << tp;
        }
      }
      if (tolerance > 0) {
        const auto& real_eigenvalues = solver.real_eigenvalues();
        EXPECT_EQ(Common::get_matrix_rows(matrix_), real_eigenvalues.size());
        for (const auto& real_ev : real_eigenvalues) {
          EXPECT_TRUE(find_ev(expected_real_eigenvalues_, real_ev, tolerance))
              << "\n\nactual eigenvalue: " << real_ev << "\n\nexpected eigenvalues: " << expected_real_eigenvalues_
              << "\n\ntype: " << tp << "\n\ntolerance: " << tolerance;
        }
      } else {
        // negative tolerance: we expect a failure
        const auto& real_eigenvalues = solver.real_eigenvalues(); /// \todo: Add try/catch around this, too!
        EXPECT_EQ(Common::get_matrix_rows(matrix_), real_eigenvalues.size());
        for (const auto& real_ev : real_eigenvalues) {
          EXPECT_FALSE(find_ev(expected_real_eigenvalues_, real_ev, 1e-15))
              << "\n\nTHIS IS A GOOD THING! UPDATE THE EXPECTATIONS IN tolerances!\n\n"
              << "\n\nactual eigenvalue: " << real_ev << "\n\nexpected eigenvalues: " << expected_real_eigenvalues_
              << "\n\ntype: " << tp;
        }
      }
    }
  } // ... gives_correct_real_eigenvalues(...)

  void gives_correct_real_eigenvalues_in_correct_order(const Common::Configuration& tolerances = {}) const
  {
    ASSERT_TRUE(all_matrices_and_expected_eigenvalues_and_vectors_are_computed_);
    for (const auto& tp : EigenSolverOpts::types()) {
      const double tolerance = tolerances.get(tp, 1e-15);
      EigenSolverType solver(matrix_, tp);
      if (tolerance > 0) {
        const auto real_eigenvalues = solver.real_eigenvalues();
        EXPECT_TRUE(Common::FloatCmp::eq(real_eigenvalues, expected_real_eigenvalues_, tolerance))
            << "\n\nactual eigenvalues: " << real_eigenvalues
            << "\n\nexpected (real) eigenvalues: " << expected_real_eigenvalues_ << "\n\ntolerance: " << tolerance
            << "\n\ntype: " << tp;
      } else {
        // negative tolerance: we expect a failure
        const auto real_eigenvalues = solver.real_eigenvalues(); /// \todo: Add try/catch around this, too!
        EXPECT_FALSE(Common::FloatCmp::eq(real_eigenvalues, expected_real_eigenvalues_))
            << "\n\nTHIS IS A GOOD THING! UPDATE THE EXPECTATIONS IN tolerances!\n\n"
            << "\n\nactual eigenvalues: " << real_eigenvalues
            << "\n\nexpected (real) eigenvalues: " << expected_real_eigenvalues_ << "\n\ntype: " << tp;
      }
    }
  } // ... gives_correct_real_eigenvalues_in_correct_order(...)

  void gives_correct_max_eigenvalue(const Common::Configuration& tolerances = {}) const
  {
    ASSERT_TRUE(all_matrices_and_expected_eigenvalues_and_vectors_are_computed_);
    for (const auto& tp : EigenSolverOpts::types()) {
      const double tolerance = tolerances.get(tp, 1e-15);
      EigenSolverType solver(matrix_, tp);
      const auto actual_max_eigenvalues = solver.max_eigenvalues(1); /// \todo: Add try/catch around this, too!
      ASSERT_GE(1, actual_max_eigenvalues.size());
      if (tolerance > 0)
        EXPECT_TRUE(Common::FloatCmp::eq(expected_max_ev_, actual_max_eigenvalues[0], tolerance))
            << "\n\nactual max eigenvalue: " << actual_max_eigenvalues[0]
            << "\n\nexpected max eigenvalue: " << expected_max_ev_ << "\n\ntolerance: " << tolerance
            << "\n\ntype: " << tp;
      else {
        // negative tolerance: we expect a failure
        EXPECT_NE(expected_max_ev_, actual_max_eigenvalues[0])
            << "\n\nTHIS IS A GOOD THING! UPDATE THE EXPECTATIONS IN tolerances!\n\n"
            << "\n\nactual max eigenvalue: " << actual_max_eigenvalues[0]
            << "\n\nexpected max eigenvalue: " << expected_max_ev_ << "\n\ntype: " << tp;
      }
    }
  }

  void gives_correct_min_eigenvalue(const Common::Configuration& tolerances = {}) const
  {
    ASSERT_TRUE(all_matrices_and_expected_eigenvalues_and_vectors_are_computed_);
    for (const auto& tp : EigenSolverOpts::types()) {
      const double tolerance = tolerances.get(tp, 1e-15);
      EigenSolverType solver(matrix_, tp);
      const auto actual_min_eigenvalues = solver.min_eigenvalues(1); /// \todo: Add try/catch around this, too!
      ASSERT_GE(1, actual_min_eigenvalues.size());
      if (tolerance > 0)
        EXPECT_TRUE(Common::FloatCmp::eq(expected_min_ev_, actual_min_eigenvalues[0], tolerance))
            << "\n\nactual min eigenvalue: " << actual_min_eigenvalues[0]
            << "\n\nexpected min eigenvalue: " << expected_min_ev_ << "\n\ntolerance: " << tolerance
            << "\n\ntype: " << tp;
      else {
        // negative tolerance: we expect a failure
        EXPECT_NE(expected_min_ev_, actual_min_eigenvalues[0])
            << "\n\nTHIS IS A GOOD THING! UPDATE THE EXPECTATIONS IN tolerances!\n\n"
            << "\n\nactual min eigenvalue: " << actual_min_eigenvalues[0]
            << "\n\nexpected min eigenvalue: " << expected_min_ev_ << "\n\ntype: " << tp;
      }
    }
  }

  void gives_correct_real_eigenvectors_in_correct_order(const Common::Configuration& tolerances = {}) const
  {
    ASSERT_TRUE(all_matrices_and_expected_eigenvalues_and_vectors_are_computed_);
    for (const auto& tp : EigenSolverOpts::types()) {
      const double tolerance = tolerances.get(tp, 1e-15);
      EigenSolverType solver(matrix_, tp);
      const auto actual_eigenvectors = solver.eigenvectors();
      const size_t rows = Common::get_matrix_rows(actual_eigenvectors);
      const size_t cols = Common::get_matrix_cols(actual_eigenvectors);
      auto actual_eigenvectors_as_real = Common::create<RealMatrixType>(rows, cols);
      for (size_t ii = 0; ii < rows; ++ii)
        for (size_t jj = 0; jj < cols; ++jj) {
          const auto complex_entry = Common::get_matrix_entry(actual_eigenvectors, ii, jj);
          if (tolerance > 0)
            EXPECT_DOUBLE_EQ(0, complex_entry.imag());
          Common::set_matrix_entry(actual_eigenvectors_as_real, ii, jj, complex_entry.real());
        }
      if (tolerance > 0) {
        EXPECT_TRUE(
            Common::FloatCmp::eq(actual_eigenvectors_as_real, expected_real_eigenvectors_, tolerance, tolerance))
            << "\n\nactual eigenvectors: " << actual_eigenvectors_as_real
            << "\n\nexpected eigenvectors: " << expected_real_eigenvectors_ << "\n\ntolerance: " << tolerance
            << "\n\ntype: " << tp;
      } else {
        EXPECT_FALSE(
            Common::FloatCmp::eq(actual_eigenvectors_as_real, expected_real_eigenvectors_, tolerance, tolerance))
            << "\n\nTHIS IS A GOOD THING! UPDATE THE EXPECTATIONS IN tolerances!\n\n"
            << "\n\nactual eigenvectors: " << actual_eigenvectors_as_real
            << "\n\nexpected eigenvectors: " << expected_real_eigenvectors_ << "\n\ntype: " << tp;
      }
      if (tolerance > 0) {
        const auto actual_real_eigenvectors = solver.real_eigenvectors();
        EXPECT_EQ(Common::get_matrix_rows(matrix_), Common::get_matrix_rows(actual_real_eigenvectors));
        EXPECT_EQ(Common::get_matrix_rows(matrix_), Common::get_matrix_cols(actual_real_eigenvectors));
        EXPECT_TRUE(Common::FloatCmp::eq(actual_real_eigenvectors, expected_real_eigenvectors_, tolerance, tolerance))
            << "\n\nactual eigenvectors: " << actual_real_eigenvectors
            << "\n\nexpected eigenvectors: " << expected_real_eigenvectors_ << "\n\ntolerance: " << tolerance
            << "\n\ntype: " << tp;
      } else {
        // negative tolerance: we expect a failure
        bool all_is_as_expected = false;
        try {
          const auto actual_real_eigenvectors = solver.real_eigenvectors();
          EXPECT_FALSE(Common::FloatCmp::eq(actual_real_eigenvectors, expected_real_eigenvectors_))
              << "\n\nactual eigenvectors: " << actual_real_eigenvectors
              << "\n\nexpected eigenvectors: " << expected_real_eigenvectors_ << "\n\ntype: " << tp;
          all_is_as_expected = true;
        } catch (const Exceptions::eigen_solver_failed_bc_eigenvalues_are_not_real_as_requested&) {
          all_is_as_expected = true;
        } catch (const Exceptions::eigen_solver_failed_bc_eigenvectors_are_not_real_as_requested&) {
          all_is_as_expected = true;
        } catch (...) {
          all_is_as_expected = true;
          FAIL() << "Expected LA::Exceptions::eigen_solver_failed_bc_eigenvalues_are_not_real_as_requested or "
                    "LA::Exceptions::eigen_solver_failed_bc_eigenvectors_are_not_real_as_requested!";
        }
        if (!all_is_as_expected)
          FAIL() << "THIS IS A GOOD THING! CHECK THE EIGENVECTORS AND ALTER THE EXPECTATIONS IN tolerances!";
      }
    }
  } // ... gives_correct_real_eigenvectors_in_correct_order(...)

  void gives_correct_real_eigendecomposition(const Common::Configuration& tolerances = {}) const
  {
    ASSERT_TRUE(all_matrices_and_expected_eigenvalues_and_vectors_are_computed_);
    for (const auto& tp : EigenSolverOpts::types()) {
      const double tolerance = tolerances.get(tp, 1e-15);
      if (tolerance > 0) {
        EigenSolverType solver(matrix_, tp);
        const auto real_eigenvalues = solver.real_eigenvalues();
        const auto T = solver.real_eigenvectors();
        auto T_inv = Common::zeros_like(matrix_);
        try {
          T_inv = invert_matrix(T);
        } catch (...) {
          FAIL() << "Matrix inversion has to work for this test, else add this matrix to the matrix inversion tests!";
        }
        auto lambda = Common::zeros_like(matrix_);
        for (size_t ii = 0; ii < Common::get_matrix_rows(matrix_); ++ii)
          Common::set_matrix_entry(lambda, ii, ii, real_eigenvalues[ii]);
        const auto matrix_decomposition_error = (T * (lambda * T_inv)) - matrix_;
        for (size_t ii = 0; ii < Common::get_matrix_rows(matrix_); ++ii)
          for (size_t jj = 0; jj < Common::get_matrix_cols(matrix_); ++jj)
            EXPECT_LT(std::abs(Common::get_matrix_entry(matrix_decomposition_error, ii, jj)), tolerance)
                << "\n\ntype = " << tp << "\n\ntolerance = " << tolerance;
      } else {
        // negative tolerance: we expect a failure
        bool exception_occured = false;
        try {
          EigenSolverType solver(matrix_, tp);
          solver.real_eigenvalues();
          solver.real_eigenvectors();
        } catch (const Exceptions::eigen_solver_failed_bc_eigenvalues_are_not_real_as_requested&) {
          exception_occured = true;
        } catch (const Exceptions::eigen_solver_failed_bc_eigenvectors_are_not_real_as_requested&) {
          exception_occured = true;
        } catch (...) {
          exception_occured = true;
          FAIL() << "Expected LA::Exceptions::eigen_solver_failed_bc_eigenvalues_are_not_real_as_requested or "
                    "LA::Exceptions::eigen_solver_failed_bc_eigenvectors_are_not_real_as_requested!"
                 << "\n\ntype: " << tp << "\n\ntolerance: " << tolerance;
        }
        if (!exception_occured)
          FAIL() << "THIS IS A GOOD THING! CHECK THE EIGENDECOMPOSITION AND ALTER THE EXPECTATIONS IN tolerances!"
                 << "\n\ntype: " << tp << "\n\ntolerance: " << tolerance;
      }
    }
  } // ... gives_correct_real_eigendecomposition(...)

  using BaseType::all_matrices_and_expected_eigenvalues_and_vectors_are_computed_;
  using BaseType::matrix_;
  using BaseType::expected_eigenvalues_;
  using BaseType::expected_eigenvectors_;
  RealEigenValuesType expected_real_eigenvalues_;
  RealType expected_max_ev_;
  RealType expected_min_ev_;
  RealMatrixType expected_real_eigenvectors_;
}; // struct EigenSolverTestForMatricesWithRealEigenvaluesAndVectors


#endif // DUNE_XT_LA_TEST_EIGENSOLVER_HH
