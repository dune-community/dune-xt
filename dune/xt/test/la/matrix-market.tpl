// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019 - 2020)

#include <dune/xt/test/main.hxx>

#include <cstdio>

#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/la/container/matrix-market.hh>
#include <dune/xt/test/la/container.hh>
#include <dune/xt/test/common/float_cmp.hh>

using namespace Dune;
using namespace Dune::XT;

{% for T_NAME, M_TYPE in config.testtypes %}
struct MatrixMarketTest_{{T_NAME}} : public ::testing::Test
{
  using MatrixImp = {{M_TYPE}};
  using M = XT::Common::MatrixAbstraction<MatrixImp>;
  using ScalarType = typename M::ScalarType;
  static constexpr size_t rows = 5;
  static constexpr size_t cols = 4;

  // Writes a matrix to file and checks that reading that file gives the original matrix
  void exporting_and_reimporting_works() const
  {
    MatrixImp mat = M::create(rows, cols);
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        M::set_entry(mat, ii, jj, static_cast<ScalarType>(ii+jj));
    const std::string filename = "testmat_{{T_NAME}}.mtx";
    LA::write_matrix_market(mat, filename);
    const MatrixImp mat2 = LA::read_matrix_market<MatrixImp>(filename);
    DXTC_EXPECT_FLOAT_EQ(mat, mat2);
    std::remove(filename.c_str());
  } // void exporting_and_reimporting_works() const

  // Tests a few matrices from the matrix-market paper (Boisvert, Ronald & Pozo, Roldan & Remington, Karin. (1997). The Matrix Market Exchange Formats: Initial Design.)
  void importing_examples_works() const
  {
    // Figure 1 of the matrix-market paper
    std::string filename1 = "mm_fig1.mtx";
    std::ofstream outputfile1(filename1);
    outputfile1 << "%%MatrixMarket   MATRIX    Coordinate   Real General" << std::endl;
    outputfile1 << "%--------------------------------------------------------------" << std::endl;
    outputfile1 << "%                   Same matrix as in example 1                " << std::endl;
    outputfile1 << "%--------------------------------------------------------------" << std::endl;
    outputfile1 << "%                                                              " << std::endl;
    outputfile1 << "%  See http://math.nist.gov/MartrixMarket for more information." << std::endl;
    outputfile1 << "                                                               " << std::endl;
    outputfile1 << "   5  5             8                                          " << std::endl;
    outputfile1 << "1 1 1.0                                                        " << std::endl;
    outputfile1 << "2 2      10.5                                                  " << std::endl;
    outputfile1 << "3 3            1.5e-2                                          " << std::endl;
    outputfile1 << "4 4                    -2.8E2                                  " << std::endl;
    outputfile1 << "5 5                              12.                           " << std::endl;
    outputfile1 << "     1     4      6                                            " << std::endl;
    outputfile1 << "     4     2      250.5                                        " << std::endl;
    outputfile1 << "     4     5      33.32                                        " << std::endl;
    outputfile1.close();
    const MatrixImp mat_fig1 = LA::read_matrix_market<MatrixImp>(filename1);
    const MatrixImp expected_mat_fig1 = XT::Common::from_string<MatrixImp>("[1 0 0 6 0; 0 10.5 0 0 0; 0 0 .015 0 0; 0 250.5 0 -280 33.32; 0 0 0 0 12]");
    DXTC_EXPECT_FLOAT_EQ(expected_mat_fig1, mat_fig1);
    std::remove(filename1.c_str());

    if (XT::Common::is_complex<ScalarType>::value) {
      // Example 2 of the matrix-market paper
      std::string filename2 = "mm_ex2.mtx";
      std::ofstream outputfile2(filename2);
      outputfile2 << "%%MatrixMarket matrix coordinate complex hermitian" << std::endl;
      outputfile2 << "  " << std::endl;
      outputfile2 << "5 5  7" << std::endl;
      outputfile2 << "1 1 1.0 0" << std::endl;
      outputfile2 << "2 2 10.5 0" << std::endl;
      outputfile2 << "4 2 250.5 22.22" << std::endl;
      outputfile2 << "3 3 1.5e-2 0" << std::endl;
      outputfile2 << "4 4 -2.8e2 0" << std::endl;
      outputfile2 << "5 5 12. 0" << std::endl;
      outputfile2 << "5 4 0 33.32" << std::endl;
      outputfile2.close();
      const MatrixImp mat_ex2 = LA::read_matrix_market<MatrixImp>(filename2);
      const MatrixImp expected_mat_ex2 = XT::Common::from_string<MatrixImp>("[1 0 0 0 0; 0 10.5 0 250.5-22.22i 0; 0 0 .015 0 0; 0 250.5+22.22i 0 -280 -33.32i; 0 0 0 33.32i 12]");
      DXTC_EXPECT_FLOAT_EQ(expected_mat_ex2, mat_ex2);
      std::remove(filename2.c_str());
    }

    // a dense symmetric matrix
    std::string filename3 = "mm_dense_symm.mtx";
    std::ofstream outputfile3(filename3);
    outputfile3 << "%%MatrixMarket matrix array real symmetric" << std::endl;
    outputfile3 << "  " << std::endl;
    outputfile3 << "3 3" << std::endl;
    outputfile3 << "1.0" << std::endl;
    outputfile3 << "2.0" << std::endl;
    outputfile3 << "3" << std::endl;
    outputfile3 << "4." << std::endl;
    outputfile3 << "5" << std::endl;
    outputfile3 << "6" << std::endl;
    outputfile3.close();
    const MatrixImp mat3 = LA::read_matrix_market<MatrixImp>(filename3);
    const MatrixImp expected_mat3 = XT::Common::from_string<MatrixImp>("[1 2 3; 2 4 5; 3 5 6]");
    DXTC_EXPECT_FLOAT_EQ(expected_mat3, mat3);
    std::remove(filename3.c_str());

    // a dense skew-symmetric matrix
    std::string filename4 = "mm_dense_skew-symm.mtx";
    std::ofstream outputfile4(filename4);
    outputfile4 << "%%MatrixMarket matrix array real skew-symmetric" << std::endl;
    outputfile4 << "  " << std::endl;
    outputfile4 << "3 3" << std::endl;
    outputfile4 << "1.0" << std::endl;
    outputfile4 << "2.0" << std::endl;
    outputfile4 << "3" << std::endl;
    outputfile4.close();
    const MatrixImp mat4 = LA::read_matrix_market<MatrixImp>(filename4);
    const MatrixImp expected_mat4 = XT::Common::from_string<MatrixImp>("[0 -1 -2; 1 0 -3; 2 3 0]");
    DXTC_EXPECT_FLOAT_EQ(expected_mat4, mat4);
    std::remove(filename4.c_str());
  } // void importing_examples_works() const
}; // struct MatrixMarketTest


TEST_F(MatrixMarketTest_{{T_NAME}}, exporting_and_reimporting_works)
{
  this->exporting_and_reimporting_works();
}

TEST_F(MatrixMarketTest_{{T_NAME}}, importing_examples_works)
{
  this->importing_examples_works();
}


{% endfor %}
