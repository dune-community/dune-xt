// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014 - 2016)
//   Tobias Leibner  (2014 - 2015, 2017)

#include <dune/xt/common/test/main.hxx> // <- has to come first (includes the config.h)!

#include <dune/xt/la/test/eigensolver.hh>

constexpr int TESTMATRIXSIZE = 5;

{% for T_NAME, TESTMATRIXTYPE, TESTFIELDTYPE, TESTCOMPLEXMATRIXTYPE, TESTREALMATRIXTYPE in config.testtypes %}
struct EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}
: public EigenSolverTestForMatricesWithRealEigenvaluesAndVectors<{{TESTMATRIXTYPE}},
        {{TESTFIELDTYPE}},
        {{TESTCOMPLEXMATRIXTYPE}},
        {{TESTREALMATRIXTYPE}}>
{
  using BaseType = EigenSolverTestForMatricesWithRealEigenvaluesAndVectors;
  using typename BaseType::MatrixType;
  using typename BaseType::ComplexMatrixType;
  using typename BaseType::RealMatrixType;
  using typename BaseType::EigenValuesType;
  using typename BaseType::RealEigenValuesType;

  EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}()
  {
    matrix_ = XT::LA::eye_matrix<MatrixType>(TESTMATRIXSIZE, TESTMATRIXSIZE);
    for (size_t ii = 0; ii < TESTMATRIXSIZE; ++ii)
      XT::Common::set_matrix_entry(matrix_, ii, ii, ii + 1);
    expected_eigenvalues_ = XT::Common::create<EigenValuesType>(TESTMATRIXSIZE);
    for (size_t ii = 0; ii < TESTMATRIXSIZE; ++ii)
      expected_eigenvalues_[ii] = ii + 1;
    expected_eigenvectors_ = XT::LA::eye_matrix<ComplexMatrixType>(TESTMATRIXSIZE, TESTMATRIXSIZE);
    expected_real_eigenvalues_ = XT::Common::create<RealEigenValuesType>(TESTMATRIXSIZE);
    for (size_t ii = 0; ii < TESTMATRIXSIZE; ++ii)
      expected_real_eigenvalues_[ii] = ii + 1;
    expected_max_ev_ = TESTMATRIXSIZE;
    expected_min_ev_ = 1;
    expected_real_eigenvectors_ = XT::LA::eye_matrix<RealMatrixType>(TESTMATRIXSIZE, TESTMATRIXSIZE);
    all_matrices_and_expected_eigenvalues_and_vectors_are_computed_ = true;
  }

  using BaseType::all_matrices_and_expected_eigenvalues_and_vectors_are_computed_;
  using BaseType::matrix_;
  using BaseType::expected_eigenvalues_;
  using BaseType::expected_eigenvectors_;
  using BaseType::expected_real_eigenvalues_;
  using BaseType::expected_max_ev_;
  using BaseType::expected_min_ev_;
  using BaseType::expected_real_eigenvectors_;
}; // struct EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}


TEST_F(EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}, exports_correct_types)
{
  exports_correct_types();
}

TEST_F(EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}, has_types_and_options)
{
  has_types_and_options();
}

TEST_F(EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}, throws_on_broken_matrix_construction)
{
  throws_on_broken_matrix_construction();
}

TEST_F(EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}, allows_broken_matrix_construction_when_checks_disabled)
{
  allows_broken_matrix_construction_when_checks_disabled();
}

TEST_F(EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}, throws_on_inconsistent_given_options)
{
  throws_on_inconsistent_given_options();
}

TEST_F(EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}, is_constructible)
{
  is_constructible();
}

TEST_F(EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}, gives_correct_eigenvalues)
{
  gives_correct_eigenvalues();
}

TEST_F(EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}, gives_correct_eigenvalues_in_correct_order)
{
  gives_correct_eigenvalues_in_correct_order();
}

TEST_F(EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}, gives_correct_real_eigenvalues)
{
  gives_correct_real_eigenvalues();
}

TEST_F(EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}, gives_correct_real_eigenvalues_in_correct_order)
{
  gives_correct_real_eigenvalues_in_correct_order();
}

TEST_F(EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}, gives_correct_max_eigenvalue)
{
  gives_correct_max_eigenvalue();
}

TEST_F(EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}, gives_correct_min_eigenvalue)
{
  gives_correct_min_eigenvalue();
}

TEST_F(EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}, gives_correct_eigenvectors_in_correct_order)
{
  gives_correct_eigenvectors_in_correct_order();
}

TEST_F(EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}, gives_correct_real_eigenvectors_in_correct_order)
{
  gives_correct_real_eigenvectors_in_correct_order();
}

TEST_F(EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}, gives_correct_eigendecomposition)
{
  gives_correct_eigendecomposition();
}

TEST_F(EigenSolverForRowWiseScaledUniMatrix_{{T_NAME}}, gives_correct_real_eigendecomposition)
{
  gives_correct_real_eigendecomposition();
}
{% endfor %}
