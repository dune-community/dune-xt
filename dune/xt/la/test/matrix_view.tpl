// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2019)

#include <dune/xt/common/test/main.hxx>
#include <dune/xt/common/test/float_cmp.hh>

#include <dune/xt/la/test/container.hh>
#include <dune/xt/la/container/matrix-view.hh>
#include <dune/xt/la/container/vector-view.hh>

using namespace Dune;
using namespace Dune::XT;

static const size_t dim = 4;
static const size_t rows = 4;
static const size_t cols = 4;

{% for T_NAME, M_TYPE, V_TYPE in config.testtypes %}
struct MatrixViewTest_{{T_NAME}} : public ::testing::Test
{
  typedef {{M_TYPE}} MatrixImp;
  typedef {{V_TYPE}} VectorImp;
  typedef typename Dune::XT::LA::SparsityPatternDefault PatternType;

  void produces_correct_results() const
  {
    typedef typename MatrixImp::ScalarType ScalarType;
    typedef typename MatrixImp::RealType RealType;

    // create test patterns
    PatternType dense_pattern(dim);
    for (size_t ii = 0; ii < dim; ++ii) {
      for (size_t jj = 0; jj < dim; ++jj)
        dense_pattern.inner(ii).push_back(jj);
    }
    PatternType sparse_pattern(dim);
    sparse_pattern.inner(0).push_back(2); //|-, -, x, -|
    sparse_pattern.inner(1).push_back(0); //|x, x, -, -|
    sparse_pattern.inner(1).push_back(1); //|-, -, -, x|
    sparse_pattern.inner(2).push_back(3); //|-, -, -, -|

    // create test matrizes
    MatrixImp testmatrix(dim, dim, dense_pattern); // |0, 1, 2, 3|
    for (size_t ii = 0; ii < dim; ++ii) { // |1, 2, 3, 4|
      for (size_t jj = 0; jj < dim; ++jj) // |2, 3, 4, 5|
        testmatrix.set_entry(ii, jj, ScalarType(ii) + ScalarType(jj)); // |3, 4, 5, 6|
    }
    MatrixImp testmatrix_sparse(dim, dim, sparse_pattern); //|-,   -, 0.5,    -|
    testmatrix_sparse.set_entry(0, 2, ScalarType(0.5)); //|1, 1.5,   -,    -|
    testmatrix_sparse.set_entry(1, 0, ScalarType(1)); //|-,   -,   -, -0.5|
    testmatrix_sparse.set_entry(1, 1, ScalarType(1.5)); //|-,   -,   -,    -|
    testmatrix_sparse.set_entry(2, 3, ScalarType(-0.5));

    // create test vector
    VectorImp vector_countingup(dim); // [0, 1, 2, 3]
    for (size_t ii = 0; ii < dim; ++ii)
      vector_countingup.set_entry(ii, ScalarType(ii));

    // create VectorViews
    LA::ConstVectorView<VectorImp> const_vec_view(vector_countingup, 0, 2);
    LA::VectorView<VectorImp> vec_view(vector_countingup, 1, 3);

    // create MatrixViews
    LA::ConstMatrixView<MatrixImp> const_view_upperleft(testmatrix, 0, 2, 0, 2);
    LA::MatrixView<MatrixImp> view_upperleft(testmatrix, 0, 2, 0, 2);
    LA::ConstMatrixView<MatrixImp> sparse_const_view_upperleft(testmatrix_sparse, 0, 2, 0, 2);
    LA::MatrixView<MatrixImp> sparse_view_upperleft(testmatrix_sparse, 0, 2, 0, 2);
    LA::ConstMatrixView<MatrixImp> const_view_lowerright(testmatrix, 2, 4, 3, 4);
    LA::MatrixView<MatrixImp> view_lowerright(testmatrix, 2, 4, 3, 4);
    LA::ConstMatrixView<MatrixImp> sparse_const_view_lowerright(testmatrix_sparse, 2, 4, 3, 4);
    LA::MatrixView<MatrixImp> sparse_view_lowerright(testmatrix_sparse, 2, 4, 3, 4);

    // test get_entry()
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(5), view_lowerright.get_entry(0, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(5), const_view_lowerright.get_entry(0, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(6), view_lowerright.get_entry(1, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(6), const_view_lowerright.get_entry(1, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(-0.5), sparse_view_lowerright.get_entry(0, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(-0.5), sparse_const_view_lowerright.get_entry(0, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(0), sparse_view_lowerright.get_entry(1, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(0), sparse_const_view_lowerright.get_entry(1, 0));

    // test add_to_entry()
    view_upperleft.add_to_entry(1, 1, 0.5);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(2.5), view_upperleft.get_entry(1, 1));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(2.5), testmatrix.get_entry(1, 1));
    view_upperleft.add_to_entry(1, 1, -0.5);
    view_lowerright.add_to_entry(0, 0, 0.5);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(5.5), view_lowerright.get_entry(0, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(5.5), testmatrix.get_entry(2, 3));
    view_lowerright.add_to_entry(0, 0, -0.5);
    sparse_view_upperleft.add_to_entry(1, 0, 0.5);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(1.5), sparse_view_upperleft.get_entry(1, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(1.5), testmatrix_sparse.get_entry(1, 0));
    sparse_view_upperleft.add_to_entry(1, 0, -0.5);

    // test set_entry()
    view_upperleft.set_entry(1, 1, 0.5);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(0.5), view_upperleft.get_entry(1, 1));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(0.5), testmatrix.get_entry(1, 1));
    view_upperleft.set_entry(1, 1, 2.);
    view_lowerright.set_entry(0, 0, 0.5);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(0.5), view_lowerright.get_entry(0, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(0.5), testmatrix.get_entry(2, 3));
    view_lowerright.set_entry(0, 0, 5.);
    sparse_view_upperleft.set_entry(1, 0, 0.5);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(0.5), sparse_view_upperleft.get_entry(1, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(0.5), testmatrix_sparse.get_entry(1, 0));
    sparse_view_upperleft.set_entry(1, 0, 1.);

    // test scal, operator*
    const auto testmatrix_copy = testmatrix;
    const auto testmatrix_sparse_copy = testmatrix_sparse;
    view_upperleft.scal(ScalarType(5));
    sparse_view_upperleft.scal(ScalarType(5));
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        if (ii < 2 && jj < 2) {
        EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(5)*ScalarType(ii + jj), view_upperleft.get_entry(ii, jj));
        EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(5)*ScalarType(ii + jj), testmatrix.get_entry(ii, jj));
        EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(5)*testmatrix_sparse_copy.get_entry(ii, jj), sparse_view_upperleft.get_entry(ii, jj));
        EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(5)*testmatrix_sparse_copy.get_entry(ii, jj), testmatrix_sparse.get_entry(ii, jj));
        } else {
          EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(ii+jj), testmatrix.get_entry(ii, jj));
          EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(ii, jj), testmatrix_sparse.get_entry(ii, jj));
        }
      }
    }
    view_upperleft.scal(ScalarType(1./5.));
    sparse_view_upperleft.scal(ScalarType(1./5.));

    view_lowerright *= ScalarType(-5);
    sparse_view_lowerright *= ScalarType(-5);
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        if (!((ii == 2 || ii == 3) && jj == 3)) {
          EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_copy.get_entry(ii, jj), testmatrix.get_entry(ii, jj));
          EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(ii, jj), testmatrix_sparse.get_entry(ii, jj));
        }
    EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_copy.get_entry(2, 3) * -5., testmatrix.get_entry(2, 3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_copy.get_entry(2, 3) * -5., view_lowerright.get_entry(0, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(2, 3) * -5., testmatrix_sparse.get_entry(2, 3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(2, 3) * -5., sparse_view_lowerright.get_entry(0, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_copy.get_entry(3, 3) * -5., testmatrix.get_entry(3, 3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_copy.get_entry(3, 3) * -5., view_lowerright.get_entry(1, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(3, 3) * -5., testmatrix_sparse.get_entry(3, 3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(3, 3) * -5., sparse_view_lowerright.get_entry(1, 0));
    view_lowerright *= ScalarType(-1./5.);
    sparse_view_lowerright *= ScalarType(-1./5.);

    // test axpy
    view_upperleft.axpy(ScalarType(4), view_upperleft);
    sparse_view_upperleft.axpy(ScalarType(4), sparse_view_upperleft);
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        if (ii < 2 && jj < 2) {
        EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(5)*ScalarType(ii + jj), view_upperleft.get_entry(ii, jj));
        EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(5)*ScalarType(ii + jj), testmatrix.get_entry(ii, jj));
        EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(5)*testmatrix_sparse_copy.get_entry(ii, jj), sparse_view_upperleft.get_entry(ii, jj));
        EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(5)*testmatrix_sparse_copy.get_entry(ii, jj), testmatrix_sparse.get_entry(ii, jj));
        } else {
          EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(ii+jj), testmatrix.get_entry(ii, jj));
          EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(ii, jj), testmatrix_sparse.get_entry(ii, jj));
        }
      }
    }
    view_upperleft.scal(ScalarType(1./5.));
    sparse_view_upperleft.scal(ScalarType(1./5.));

    view_lowerright.axpy(ScalarType(-6), view_lowerright);
    sparse_view_lowerright.axpy(ScalarType(-6), sparse_view_lowerright);
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        if (!((ii == 2 || ii == 3) && jj == 3)) {
          EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_copy.get_entry(ii, jj), testmatrix.get_entry(ii, jj));
          EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(ii, jj), testmatrix_sparse.get_entry(ii, jj));
        }
    EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_copy.get_entry(2, 3) * -5., testmatrix.get_entry(2, 3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_copy.get_entry(2, 3) * -5., view_lowerright.get_entry(0, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(2, 3) * -5., testmatrix_sparse.get_entry(2, 3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(2, 3) * -5., sparse_view_lowerright.get_entry(0, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_copy.get_entry(3, 3) * -5., testmatrix.get_entry(3, 3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_copy.get_entry(3, 3) * -5., view_lowerright.get_entry(1, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(3, 3) * -5., testmatrix_sparse.get_entry(3, 3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(3, 3) * -5., sparse_view_lowerright.get_entry(1, 0));
    view_lowerright *= ScalarType(-1./5.);
    sparse_view_lowerright *= ScalarType(-1./5.);

    // test mv, mtv
    VectorImp testvec(2), ret(2), ret2(4), ret3(2);
    LA::VectorView<VectorImp> ret2_view(ret2, 2, 4);
    testvec[0] = 0.;
    testvec[1] = 1.;
    view_upperleft.mv(testvec, ret);
    view_upperleft.mv(const_vec_view, ret2_view);
    sparse_view_upperleft.mv(vec_view, ret3);
    for (size_t ii = 0; ii < 2; ++ii) {
      EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(ii+1), ret[ii]);
      EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(ii+1), ret2[ii+2]);
      EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(ii+1), ret2_view[ii]);
      EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(4*ii), ret3[ii]);
    }
    view_upperleft.mtv(testvec, ret);
    view_upperleft.mtv(const_vec_view, ret2_view);
    sparse_view_upperleft.mtv(vec_view, ret3);
    for (size_t ii = 0; ii < 2; ++ii) {
      EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(ii+1), ret[ii]);
      EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(ii+1), ret2[ii+2]);
      EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(ii+1), ret2_view[ii]);
      EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(2+ii), ret3[ii]);
    }

    // test clear_row/col, unit_row/col
    view_upperleft.clear_row(0);
    sparse_view_upperleft.clear_row(0);
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        if (ii < 2 && jj < 2) {
          if (ii == 0) {
            EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(0), view_upperleft.get_entry(ii, jj));
            EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(0), testmatrix.get_entry(ii, jj));
            EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(0), sparse_view_upperleft.get_entry(ii, jj));
            EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(0), testmatrix_sparse.get_entry(ii, jj));
          } else {
            EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(ii + jj), view_upperleft.get_entry(ii, jj));
            EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(ii + jj), testmatrix.get_entry(ii, jj));
            EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(ii, jj), sparse_view_upperleft.get_entry(ii, jj));
            EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(ii, jj), testmatrix_sparse.get_entry(ii, jj));
          }
        } else {
          EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(ii+jj), testmatrix.get_entry(ii, jj));
          EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(ii, jj), testmatrix_sparse.get_entry(ii, jj));
        }
      }
    }
    view_upperleft.unit_row(1);
    sparse_view_upperleft.unit_row(1);
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        if (ii < 2 && jj < 2) {
          if (ii == 0) {
            EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(0), view_upperleft.get_entry(ii, jj));
            EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(0), testmatrix.get_entry(ii, jj));
            EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(0), sparse_view_upperleft.get_entry(ii, jj));
            EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(0), testmatrix_sparse.get_entry(ii, jj));
          } else if (ii == 1) {
            EXPECT_DOUBLE_OR_COMPLEX_EQ(jj == 1 ? ScalarType(1) : ScalarType(0), view_upperleft.get_entry(ii, jj));
            EXPECT_DOUBLE_OR_COMPLEX_EQ(jj == 1 ? ScalarType(1) : ScalarType(0), testmatrix.get_entry(ii, jj));
            EXPECT_DOUBLE_OR_COMPLEX_EQ(jj == 1 ? ScalarType(1) : ScalarType(0), sparse_view_upperleft.get_entry(ii, jj));
            EXPECT_DOUBLE_OR_COMPLEX_EQ(jj == 1 ? ScalarType(1) : ScalarType(0), testmatrix_sparse.get_entry(ii, jj));
          } else {
            EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(ii + jj), view_upperleft.get_entry(ii, jj));
            EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(ii + jj), testmatrix.get_entry(ii, jj));
            EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(ii, jj), sparse_view_upperleft.get_entry(ii, jj));
            EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(ii, jj), testmatrix_sparse.get_entry(ii, jj));
          }
        } else {
          EXPECT_DOUBLE_OR_COMPLEX_EQ(ScalarType(ii+jj), testmatrix.get_entry(ii, jj));
          EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy.get_entry(ii, jj), testmatrix_sparse.get_entry(ii, jj));
        }
      }
    }

    auto testmatrix_copy2 = testmatrix;
    auto testmatrix_sparse_copy2 = testmatrix_sparse;
    view_lowerright.clear_col(0);
    sparse_view_lowerright.clear_col(0);
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        if (!((ii == 2 || ii == 3) && jj == 3)) {
          EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_copy2.get_entry(ii, jj), testmatrix.get_entry(ii, jj));
          EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy2.get_entry(ii, jj), testmatrix_sparse.get_entry(ii, jj));
        }
    EXPECT_DOUBLE_OR_COMPLEX_EQ(0., testmatrix.get_entry(2, 3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(0., view_lowerright.get_entry(0, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(0., testmatrix_sparse.get_entry(2, 3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(0., sparse_view_lowerright.get_entry(0, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(0., testmatrix.get_entry(3, 3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(0., view_lowerright.get_entry(1, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(0., testmatrix_sparse.get_entry(3, 3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(0., sparse_view_lowerright.get_entry(1, 0));

    view_lowerright.unit_col(0);
    sparse_view_lowerright.unit_col(0);
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        if (!((ii == 2 || ii == 3) && jj == 3)) {
          EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_copy2.get_entry(ii, jj), testmatrix.get_entry(ii, jj));
          EXPECT_DOUBLE_OR_COMPLEX_EQ(testmatrix_sparse_copy2.get_entry(ii, jj), testmatrix_sparse.get_entry(ii, jj));
        }
    EXPECT_DOUBLE_OR_COMPLEX_EQ(1., testmatrix.get_entry(2, 3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(1., view_lowerright.get_entry(0, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(1., testmatrix_sparse.get_entry(2, 3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(1., sparse_view_lowerright.get_entry(0, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(0., testmatrix.get_entry(3, 3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(0., view_lowerright.get_entry(1, 0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(0., testmatrix_sparse.get_entry(3, 3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(0., sparse_view_lowerright.get_entry(1, 0));
  } // void produces_correct_results() const
}; // struct MatrixTest

TEST_F(MatrixViewTest_{{T_NAME}}, produces_correct_results)
{
  this->produces_correct_results();
}

{% endfor %}
