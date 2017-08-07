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

// This one has to come first (includes the config.h)!
#include <dune/xt/common/test/main.hxx>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/la/eigen-solver.hh>

// toggle output
// std::ostream& out = std::cout;
std::ostream& out = DXTC_LOG.devnull();

using namespace Dune;
using namespace Dune::XT;
using namespace Dune::XT::LA;

struct EigenSolverTest : public ::testing::Test
{
  typedef TESTMATRIXTYPE MatrixType;
  typedef COMPLEXTESTMATRIXTYPE ComplexMatrixType;

  typedef EigenSolver<MatrixType> EigenSolverType;

  static void produces_correct_results()
  {
    const size_t size = 2;
    // eigenvalues 0 and 2
    const MatrixType symmetric_matrix = XT::Common::from_string<MatrixType>("[1 1; 1 1]");
    // eigenvalues 1 + i and 1 - i
    const MatrixType nonsymmetric_matrix = XT::Common::from_string<MatrixType>("[1 1; -1 1]");
    const ComplexMatrixType nonsymmetric_matrix_complex = XT::Common::from_string<ComplexMatrixType>("[1 1; -1 1]");

    std::vector<std::string> types = EigenSolverType::types();
    if (types.size() == 0)
      DUNE_THROW(Common::Exceptions::results_are_not_as_expected, "Solver has no types!");
    for (auto type : types) {
      out << "solving with type '" << type << "' and options" << std::endl;
      Common::Configuration options = EigenSolverType::options(type);
      options.report(out, "  ");

      const EigenSolverType eigensolver1(symmetric_matrix);
      std::vector<std::complex<double>> eigenvalues1 = eigensolver1.eigenvalues(options);
      for (const auto& eigval : eigenvalues1)
        EXPECT_DOUBLE_EQ(eigval.imag(), 0.);
      std::vector<double> eigenvalues1_real = {eigenvalues1[0].real(), eigenvalues1[1].real()};
      std::vector<typename EigenSolverType::ComplexVectorType> DUNE_UNUSED(eigenvectors1) =
          eigensolver1.eigenvectors(options);
      std::vector<typename EigenSolverType::RealVectorType> eigenvectors1_real =
          eigensolver1.real_eigenvectors(options);
      for (size_t ii = 0; ii < size; ++ii) {
        const auto& eigval = eigenvalues1_real[ii];
        const auto& eigvec = eigenvectors1_real[ii];
        auto zero_vec = eigvec;
        zero_vec *= 0.;
        EXPECT_FALSE(XT::Common::FloatCmp::eq(eigvec, zero_vec));
        auto result = zero_vec;
        symmetric_matrix.mv(eigvec, result);
        for (size_t jj = 0; jj < size; ++jj)
          EXPECT_DOUBLE_EQ(eigval * eigvec[jj], result[jj]);
      } // ii
      std::sort(eigenvalues1_real.begin(), eigenvalues1_real.end());
      EXPECT_DOUBLE_EQ(eigenvalues1_real[0], eigensolver1.min_eigenvalues(1)[0]);
      EXPECT_DOUBLE_EQ(eigenvalues1_real[1], eigensolver1.max_eigenvalues(1)[0]);
      EXPECT_DOUBLE_EQ(0, eigenvalues1_real[0]);
      EXPECT_DOUBLE_EQ(2, eigenvalues1_real[1]);
      std::shared_ptr<typename EigenSolverType::ComplexMatrixType> eigenvectors1_matrix =
          eigensolver1.eigenvectors_as_matrix(options);
      std::shared_ptr<typename EigenSolverType::RealMatrixType> eigenvectors1_realmatrix =
          eigensolver1.real_eigenvectors_as_matrix(options);
      for (size_t ii = 0; ii < size; ++ii) {
        for (size_t jj = 0; jj < size; ++jj) {
          EXPECT_DOUBLE_EQ(eigenvectors1[jj][ii].real(),
                           XT::Common::MatrixAbstraction<typename EigenSolverType::ComplexMatrixType>::get_entry(
                               *eigenvectors1_matrix, ii, jj)
                               .real());
          EXPECT_DOUBLE_EQ(eigenvectors1[jj][ii].imag(),
                           XT::Common::MatrixAbstraction<typename EigenSolverType::ComplexMatrixType>::get_entry(
                               *eigenvectors1_matrix, ii, jj)
                               .imag());
          EXPECT_DOUBLE_EQ(eigenvectors1_real[jj][ii],
                           XT::Common::MatrixAbstraction<typename EigenSolverType::RealMatrixType>::get_entry(
                               *eigenvectors1_realmatrix, ii, jj));
        }
      }

      const EigenSolverType eigensolver2(nonsymmetric_matrix);
      std::vector<std::complex<double>> eigenvalues2 = eigensolver2.eigenvalues();
      std::vector<typename EigenSolverType::ComplexVectorType> eigenvectors2 = eigensolver2.eigenvectors(options);
      for (size_t ii = 0; ii < size; ++ii) {
        const auto& eigval = eigenvalues2[ii];
        const auto& eigvec = eigenvectors2[ii];
        auto zero_vec = eigvec;
        zero_vec *= 0.;
        EXPECT_FALSE(XT::Common::FloatCmp::eq(eigvec, zero_vec));
        auto result = zero_vec;
        nonsymmetric_matrix_complex.mv(eigvec, result);
        for (size_t jj = 0; jj < size; ++jj) {
          EXPECT_DOUBLE_EQ((eigval * eigvec[jj]).real(), result[jj].real());
          EXPECT_DOUBLE_EQ((eigval * eigvec[jj]).imag(), result[jj].imag());
        }
      } // ii
      std::shared_ptr<typename EigenSolverType::ComplexMatrixType> eigenvectors2_matrix =
          eigensolver2.eigenvectors_as_matrix(options);
      for (size_t ii = 0; ii < size; ++ii) {
        for (size_t jj = 0; jj < size; ++jj) {
          EXPECT_DOUBLE_EQ(eigenvectors2[jj][ii].real(),
                           XT::Common::MatrixAbstraction<typename EigenSolverType::ComplexMatrixType>::get_entry(
                               *eigenvectors2_matrix, ii, jj)
                               .real());
          EXPECT_DOUBLE_EQ(eigenvectors2[jj][ii].imag(),
                           XT::Common::MatrixAbstraction<typename EigenSolverType::ComplexMatrixType>::get_entry(
                               *eigenvectors2_matrix, ii, jj)
                               .imag());
        }
      }
      std::sort(eigenvalues2.begin(), eigenvalues2.end(), [](std::complex<double> a, std::complex<double> b) {
        return a.imag() < b.imag();
      });
      for (const auto& eigval : eigenvalues2)
        EXPECT_DOUBLE_EQ(eigval.real(), 1.);
      EXPECT_DOUBLE_EQ(eigenvalues2[0].imag(), -1.);
      EXPECT_DOUBLE_EQ(eigenvalues2[1].imag(), 1.);
    } // types
  } // ... produces_correct_results(...)
}; // struct SolverTest

TEST_F(EigenSolverTest, behaves_correctly)
{
  this->produces_correct_results();
}
