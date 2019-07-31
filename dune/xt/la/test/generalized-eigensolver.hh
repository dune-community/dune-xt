// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_XT_LA_TEST_GENERALIZED_EIGENSOLVER_HH
#define DUNE_XT_LA_TEST_GENERALIZED_EIGENSOLVER_HH

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/unused.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/test/gtest/gtest.h>

#include <dune/xt/la/container.hh>
#include <dune/xt/la/container/conversion.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/la/generalized-eigen-solver.hh>

using namespace Dune;
using namespace Dune::XT;
using namespace Dune::XT::LA;


template <class MatrixImp, class F, class C, class R>
struct GeneralizedEigenSolverTest : public ::testing::Test
{
  using MatrixType = MatrixImp;
  using FieldType = Common::real_t<F>;
  using ComplexMatrixType = C;
  using RealMatrixType = R;
  using RealType = Common::real_t<FieldType>;
  using EigenValuesType = std::vector<Common::complex_t<FieldType>>;
  using RealEigenValuesType = std::vector<Common::real_t<FieldType>>;
  typedef GeneralizedEigenSolver<MatrixType> GeneralizedEigenSolverType;
  typedef GeneralizedEigenSolverOptions<MatrixType> GeneralizedEigenSolverOpts;

  using M = Common::MatrixAbstraction<MatrixType>;

  GeneralizedEigenSolverTest()
    : all_matrices_and_expected_eigenvalues_and_vectors_are_computed_(false)
    , broken_matrix_(M::create(M::has_static_size ? M::static_rows : 1, M::has_static_size ? M::static_cols : 1))
    , unit_matrix_(
          eye_matrix<MatrixType>(M::has_static_size ? M::static_rows : 1, M::has_static_size ? M::static_cols : 1))
  {
    M::set_entry(broken_matrix_, 0, 0, std::numeric_limits<typename M::S>::infinity());
  }

  void make_unit_matrix()
  {
    ASSERT_TRUE(all_matrices_and_expected_eigenvalues_and_vectors_are_computed_);
    unit_matrix_ = eye_matrix<MatrixType>(M::has_static_size ? M::static_rows : M::rows(matrix_),
                                          M::has_static_size ? M::static_cols : M::cols(matrix_));
  }

  static void exports_correct_types()
  {
    const bool MatrixType_is_correct = std::is_same<typename GeneralizedEigenSolverType::MatrixType, MatrixType>::value;
    EXPECT_TRUE(MatrixType_is_correct);
    const bool RealType_is_correct = std::is_same<typename GeneralizedEigenSolverType::RealType, RealType>::value;
    EXPECT_TRUE(RealType_is_correct);
    const bool RealMatrixType_is_correct =
        std::is_same<typename GeneralizedEigenSolverType::RealMatrixType, RealMatrixType>::value;
    EXPECT_TRUE(RealMatrixType_is_correct);
    const bool ComplexMatrixType_is_correct =
        std::is_same<typename GeneralizedEigenSolverType::ComplexMatrixType, ComplexMatrixType>::value;
    EXPECT_TRUE(ComplexMatrixType_is_correct);
  } // ... exports_correct_types()

  static void has_types_and_options()
  {
    std::vector<std::string> types = GeneralizedEigenSolverOpts::types();
    EXPECT_GT(types.size(), 0);
    Common::Configuration opts = GeneralizedEigenSolverOpts::options();
    EXPECT_TRUE(opts.has_key("type"));
    EXPECT_EQ(types[0], opts.get<std::string>("type"));
    for (const auto& tp : types) {
      Common::Configuration tp_opts = GeneralizedEigenSolverOpts::options(tp);
      EXPECT_TRUE(tp_opts.has_key("type"));
      EXPECT_EQ(tp, tp_opts.get<std::string>("type"));
    }
  } // ... has_types_and_options(...)

  void throws_on_broken_matrix_construction()
  {
    try {
      GeneralizedEigenSolverType DXTC_UNUSED(default_solver){broken_matrix_, unit_matrix_};
      FAIL() << "Expected LA::Exceptions::generalized_eigen_solver_failed_bc_data_did_not_fulfill_requirements";
    } catch (const LA::Exceptions::generalized_eigen_solver_failed_bc_data_did_not_fulfill_requirements& /*ee*/) {
    } catch (...) {
      FAIL() << "Expected LA::Exceptions::generalized_eigen_solver_failed_bc_data_did_not_fulfill_requirements";
    }
    for (const auto& tp : GeneralizedEigenSolverOpts::types()) {
      try {
        GeneralizedEigenSolverType DXTC_UNUSED(default_opts_solver)(broken_matrix_, unit_matrix_, tp);
        FAIL() << "Expected LA::Exceptions::generalized_eigen_solver_failed_bc_data_did_not_fulfill_requirements";
      } catch (const LA::Exceptions::generalized_eigen_solver_failed_bc_data_did_not_fulfill_requirements& /*ee*/) {
      } catch (...) {
        FAIL() << "Expected LA::Exceptions::generalized_eigen_solver_failed_bc_data_did_not_fulfill_requirements";
      }
      try {
        GeneralizedEigenSolverType DXTC_UNUSED(solver)(
            broken_matrix_, unit_matrix_, GeneralizedEigenSolverOpts::options(tp));
        FAIL() << "Expected LA::Exceptions::generalized_eigen_solver_failed_bc_data_did_not_fulfill_requirements";
      } catch (const LA::Exceptions::generalized_eigen_solver_failed_bc_data_did_not_fulfill_requirements& /*ee*/) {
      } catch (...) {
        FAIL() << "Expected LA::Exceptions::generalized_eigen_solver_failed_bc_data_did_not_fulfill_requirements";
      }
    }
  } // ... throws_on_broken_matrix_construction(...)

  void allows_broken_matrix_construction_when_checks_disabled()
  {
    for (const auto& tp : GeneralizedEigenSolverOpts::types()) {
      Common::Configuration opts_with_disabled_check = GeneralizedEigenSolverOpts::options(tp);
      opts_with_disabled_check["check_for_inf_nan"] = "false";
      GeneralizedEigenSolverType DXTC_UNUSED(solver)(broken_matrix_, unit_matrix_, opts_with_disabled_check);
    }
  }

  void throws_on_inconsistent_given_options()
  {
    for (const auto& tp : GeneralizedEigenSolverOpts::types()) {
      Common::Configuration inconsistent_opts = GeneralizedEigenSolverOpts::options(tp);
      inconsistent_opts["assert_positive_eigenvalues"] = "1e-15";
      inconsistent_opts["assert_negative_eigenvalues"] = "1e-15";
      try {
        GeneralizedEigenSolverType DXTC_UNUSED(solver)(unit_matrix_, unit_matrix_, inconsistent_opts);
        FAIL() << "Expected LA::Exceptions::generalized_eigen_solver_failed_bc_it_was_not_set_up_correctly";
      } catch (const LA::Exceptions::generalized_eigen_solver_failed_bc_it_was_not_set_up_correctly& /*ee*/) {
      } catch (...) {
        FAIL() << "Expected LA::Exceptions::generalized_eigen_solver_failed_bc_it_was_not_set_up_correctly";
      }
    }
  } // ... throws_on_inconsistent_given_options(...)

  void is_constructible()
  {
    ASSERT_TRUE(all_matrices_and_expected_eigenvalues_and_vectors_are_computed_);
    this->make_unit_matrix();
    GeneralizedEigenSolverType DXTC_UNUSED(default_solver){matrix_, unit_matrix_};
    for (const auto& tp : GeneralizedEigenSolverOpts::types()) {
      GeneralizedEigenSolverType DXTC_UNUSED(default_opts_solver)(matrix_, unit_matrix_, tp);
      GeneralizedEigenSolverType DXTC_UNUSED(solver)(matrix_, unit_matrix_, GeneralizedEigenSolverOpts::options(tp));
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

  void gives_correct_eigenvalues(const Common::Configuration& tolerances = {})
  {
    ASSERT_TRUE(all_matrices_and_expected_eigenvalues_and_vectors_are_computed_);
    this->make_unit_matrix();
    for (const auto& tp : GeneralizedEigenSolverOpts::types()) {
      const double tolerance = tolerances.get(tp, 1e-15);
      GeneralizedEigenSolverType solver(matrix_, unit_matrix_, tp);
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

  bool all_matrices_and_expected_eigenvalues_and_vectors_are_computed_;
  MatrixType broken_matrix_;
  MatrixType unit_matrix_;
  MatrixType matrix_;
  EigenValuesType expected_eigenvalues_;
}; // ... struct GeneralizedEigenSolverTest


template <class M, class F, class C, class R>
struct GeneralizedEigenSolverTestForMatricesWithRealEigenvaluesAndVectors
  : public GeneralizedEigenSolverTest<M, F, C, R>
{
  using BaseType = GeneralizedEigenSolverTest<M, F, C, R>;
  using BaseType::find_ev;
  using typename BaseType::ComplexMatrixType;
  using typename BaseType::EigenValuesType;
  using typename BaseType::FieldType;
  using typename BaseType::GeneralizedEigenSolverOpts;
  using typename BaseType::GeneralizedEigenSolverType;
  using typename BaseType::MatrixType;
  using typename BaseType::RealEigenValuesType;
  using typename BaseType::RealMatrixType;
  using typename BaseType::RealType;

  void gives_correct_real_eigenvalues(const Common::Configuration& tolerances = {})
  {
    ASSERT_TRUE(all_matrices_and_expected_eigenvalues_and_vectors_are_computed_);
    this->make_unit_matrix();
    for (const auto& tp : GeneralizedEigenSolverOpts::types()) {
      const double tolerance = tolerances.get(tp, 1e-15);
      GeneralizedEigenSolverType solver(matrix_, unit_matrix_, tp);
      const auto& eigenvalues = solver.eigenvalues();
      EXPECT_EQ(Common::get_matrix_rows(matrix_), eigenvalues.size());
      for (const auto& complex_ev : eigenvalues) {
        if (tolerance > 0)
          EXPECT_TRUE(Common::FloatCmp::eq(0., complex_ev.imag(), tolerance))
              << "\n  type: " << tp << "\n  tolerance: " << tolerance;
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

  void gives_correct_max_eigenvalue(const Common::Configuration& tolerances = {})
  {
    ASSERT_TRUE(all_matrices_and_expected_eigenvalues_and_vectors_are_computed_);
    this->make_unit_matrix();
    for (const auto& tp : GeneralizedEigenSolverOpts::types()) {
      const double tolerance = tolerances.get(tp, 1e-15);
      GeneralizedEigenSolverType solver(matrix_, unit_matrix_, tp);
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

  void gives_correct_min_eigenvalue(const Common::Configuration& tolerances = {})
  {
    ASSERT_TRUE(all_matrices_and_expected_eigenvalues_and_vectors_are_computed_);
    this->make_unit_matrix();
    for (const auto& tp : GeneralizedEigenSolverOpts::types()) {
      const double tolerance = tolerances.get(tp, 1e-15);
      GeneralizedEigenSolverType solver(matrix_, unit_matrix_, tp);
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

  using BaseType::all_matrices_and_expected_eigenvalues_and_vectors_are_computed_;
  using BaseType::expected_eigenvalues_;
  using BaseType::matrix_;
  using BaseType::unit_matrix_;
  RealEigenValuesType expected_real_eigenvalues_;
  RealType expected_max_ev_;
  RealType expected_min_ev_;
}; // struct GeneralizedEigenSolverTestForMatricesWithRealEigenvaluesAndVectors


#endif // DUNE_XT_LA_TEST_GENERALIZED_EIGENSOLVER_HH
