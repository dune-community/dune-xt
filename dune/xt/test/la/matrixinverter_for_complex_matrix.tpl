// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014 - 2016)
//   Tobias Leibner  (2014 - 2015, 2017)

#include <dune/xt/common/test/main.hxx> // <- has to come first (includes the config.h)!

#include <dune/xt/la/test/matrixinverter.hh>

{% for T_NAME, TESTMATRIXTYPE in config.testtypes %}
struct MatrixInverterForComplexMatrix_{{T_NAME}}
: public MatrixInverterTest<{{TESTMATRIXTYPE}}>
{
  using BaseType = MatrixInverterTest;
  using typename BaseType::MatrixType;

  /**
       Tests inversion of a simple complex matrix.
        "[1+i  i;
        " i-1 2i]"
       The inverse is
        "[ 0.6-0.2i  -0.3+0.1i;"
        " -0.4-0.2i   0.2-0.4i]"
   */
  MatrixInverterForComplexMatrix_{{T_NAME}}()
  {
    matrix_ = XT::Common::from_string<MatrixType>("[1+1i 0+1i; -1+1i 0+2i]");
    expected_inverse_ = XT::Common::from_string<MatrixType>("[0.6-0.2i -0.3+0.1i; -0.4-0.2i 0.2-0.4i]");
    all_matrices_and_inverse_matrices_are_computed_ = true;
  }

  using BaseType::all_matrices_and_inverse_matrices_are_computed_;
  using BaseType::matrix_;
  using BaseType::expected_inverse_;
}; // struct MatrixInverterForComplexMatrix_{{T_NAME}}

TEST_F(MatrixInverterForComplexMatrix_{{T_NAME}}, exports_correct_types)
{
  exports_correct_types();
}

TEST_F(MatrixInverterForComplexMatrix_{{T_NAME}}, has_types_and_options)
{
  has_types_and_options();
}

TEST_F(MatrixInverterForComplexMatrix_{{T_NAME}}, throws_on_broken_matrix_construction)
{
  throws_on_broken_matrix_construction();
}

TEST_F(MatrixInverterForComplexMatrix_{{T_NAME}}, is_constructible)
{
  is_constructible();
}

TEST_F(MatrixInverterForComplexMatrix_{{T_NAME}}, gives_correct_inverse)
{
  gives_correct_inverse({ {"direct", "1e-13"} });
}

{% endfor %}
