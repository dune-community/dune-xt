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

#include "eigensolver.hh"


/**
 * \sa http://eigen.tuxfamily.org/dox/classEigen_1_1EigenSolver.html
 *
 * \note Ours might be quite a bit off, since we only have the truncated representations from the example.
 *
\code
Example:
MatrixXd A = MatrixXd::Random(6,6);
cout << "Here is a random 6x6 matrix, A:" << endl << A << endl << endl;
EigenSolver<MatrixXd> es(A);
cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;
complex<double> lambda = es.eigenvalues()[0];
cout << "Consider the first eigenvalue, lambda = " << lambda << endl;
VectorXcd v = es.eigenvectors().col(0);
cout << "If v is the corresponding eigenvector, then lambda * v = " << endl << lambda * v << endl;
cout << "... and A * v = " << endl << A.cast<complex<double> >() * v << endl << endl;
MatrixXcd D = es.eigenvalues().asDiagonal();
MatrixXcd V = es.eigenvectors();
cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;

Output:

Here is a random 6x6 matrix, A:
   0.68   -0.33   -0.27  -0.717  -0.687  0.0259
 -0.211   0.536  0.0268   0.214  -0.198   0.678
  0.566  -0.444   0.904  -0.967   -0.74   0.225
  0.597   0.108   0.832  -0.514  -0.782  -0.408
  0.823 -0.0452   0.271  -0.726   0.998   0.275
 -0.605   0.258   0.435   0.608  -0.563  0.0486

The eigenvalues of A are:
  (0.049,1.06)
 (0.049,-1.06)
     (0.967,0)
     (0.353,0)
 (0.618,0.129)
(0.618,-0.129)
The matrix of eigenvectors, V, is:
 (-0.292,-0.454)   (-0.292,0.454)      (-0.0607,0)       (-0.733,0)    (0.59,-0.121)     (0.59,0.121)
  (0.134,-0.104)    (0.134,0.104)       (-0.799,0)        (0.136,0)    (0.334,0.368)   (0.334,-0.368)
  (-0.422,-0.18)    (-0.422,0.18)        (0.192,0)       (0.0563,0)  (-0.335,-0.143)   (-0.335,0.143)
 (-0.589,0.0274) (-0.589,-0.0274)      (-0.0788,0)       (-0.627,0)   (0.322,-0.155)    (0.322,0.155)
  (-0.248,0.132)  (-0.248,-0.132)        (0.401,0)        (0.218,0) (-0.335,-0.0761)  (-0.335,0.0761)
    (0.105,0.18)    (0.105,-0.18)       (-0.392,0)     (-0.00564,0)  (-0.0324,0.103) (-0.0324,-0.103)

Consider the first eigenvalue, lambda = (0.049,1.06)
If v is the corresponding eigenvector, then lambda * v =
  (0.466,-0.331)
   (0.117,0.137)
   (0.17,-0.456)
(-0.0578,-0.622)
 (-0.152,-0.256)
   (-0.186,0.12)
... and A * v =
  (0.466,-0.331)
   (0.117,0.137)
   (0.17,-0.456)
(-0.0578,-0.622)
 (-0.152,-0.256)
   (-0.186,0.12)

Finally, V * D * V^(-1) =
   (0.68,4.44e-16)          (-0.33,0)  (-0.27,-3.33e-16)  (-0.717,2.22e-16)  (-0.687,4.44e-16)         (0.0259,0)
(-0.211,-4.44e-16)  (0.536,-1.05e-16)         (0.0268,0)   (0.214,6.66e-16)  (-0.198,4.44e-16)   (0.678,8.88e-16)
  (0.566,1.55e-15)  (-0.444,2.08e-16)  (0.904,-3.89e-16) (-0.967,-5.55e-16)   (-0.74,1.33e-15)   (0.225,6.66e-16)
  (0.597,1.33e-15)   (0.108,3.61e-16)  (0.832,-5.55e-16) (-0.514,-9.99e-16)  (-0.782,8.88e-16)  (-0.408,4.44e-16)
  (0.823,2.22e-16) (-0.0452,1.67e-16)  (0.271,-1.67e-16) (-0.726,-7.77e-16)          (0.998,0)   (0.275,4.44e-16)
(-0.605,-1.39e-16)   (0.258,4.16e-17)   (0.435,8.33e-17)   (0.608,7.29e-17)  (-0.563,5.83e-16)  (0.0486,2.29e-16)
\endcode
 */
struct EigenSolverForMatrixFromEigensExample : public EigenSolverTest
{
  using BaseType = EigenSolverTest;
  using typename BaseType::MatrixType;
  using typename BaseType::ComplexMatrixType;
  using typename BaseType::EigenValuesType;

  EigenSolverForMatrixFromEigensExample()
  {
    matrix_ = XT::Common::from_string<MatrixType>("[ 0.68  -0.33   -0.27   -0.717  -0.687  0.0259;"
                                                  " -0.211  0.536   0.0268  0.214  -0.198  0.678;"
                                                  "  0.566 -0.444   0.904  -0.967  -0.74   0.225;"
                                                  "  0.597  0.108   0.832  -0.514  -0.782 -0.408;"
                                                  "  0.823 -0.0452  0.271  -0.726   0.998  0.275;"
                                                  " -0.605  0.258   0.435   0.608  -0.563  0.0486]");
    // these values were recorded using lapacke
    expected_eigenvalues_ = XT::Common::from_string<EigenValuesType>("[0.048829697485156004+1.0596562712819182i "
                                                                     " 0.048829697485156004-1.0596562712819182i "
                                                                     " 0.96873863535292992 "
                                                                     " 0.34280071327659634 "
                                                                     " 0.62170062820008098+0.13733227413092269i "
                                                                     " 0.62170062820008098-0.13733227413092269i]");
    expected_eigenvectors_ = XT::Common::from_string<ComplexMatrixType>(
        "["
        "-0.27012632227983902-0.46678683793394776i -0.27012632227983902+0.46678683793394776i -0.058555145136596501 "
        "-0.72934250003563827 -0.59311240459518533 -0.59311240459518533; "
        " 0.13887686287144707-0.097771882383961992i 0.13887686287144707+0.097771882383961992i -0.79902060638186234   "
        "0.14851614296331883 -0.23987153829467844-0.45143010420135704i -0.23987153829467844+0.45143010420135704i; "
        " -0.41364142584806307-0.19917683735545064i -0.41364142584806307+0.19917683735545064i 0.18949926740896209 "
        "0.04700276597877924 0.29529738066681194+0.21971404267746528i 0.29529738066681194-0.21971404267746528i; "
        " -0.58911879728583405 -0.58911879728583405 -0.078813329494931098 -0.63174635015851877 "
        "-0.33835579668284677+0.092082191690773141i -0.33835579668284677-0.092082191690773141i; "
        " -0.2538014452885316+0.12001744521279049i -0.2538014452885316-0.12001744521279049i 0.40135933686903708 "
        "0.21115859793042213 0.30424687115889376+0.14919085125463472i 0.30424687115889376-0.14919085125463472i; "
        " 0.097132120486031412+0.18495968552407643i  0.097132120486031412-0.18495968552407643i -0.39360708201025785 "
        "-0.010090019467593912 0.060511969188928441-0.099829032727048669i 0.060511969188928441+0.099829032727048669i "
        "]");
    all_matrices_and_expected_eigenvalues_and_vectors_are_computed_ = true;
  }

  using BaseType::all_matrices_and_expected_eigenvalues_and_vectors_are_computed_;
  using BaseType::matrix_;
  using BaseType::expected_eigenvalues_;
  using BaseType::expected_eigenvectors_;
}; // struct EigenSolverForMatrixFromEigensExample


TEST_F(EigenSolverForMatrixFromEigensExample, exports_correct_types)
{
  exports_correct_types();
}

TEST_F(EigenSolverForMatrixFromEigensExample, has_types_and_options)
{
  has_types_and_options();
}

TEST_F(EigenSolverForMatrixFromEigensExample, throws_on_broken_matrix_construction)
{
  throws_on_broken_matrix_construction();
}

TEST_F(EigenSolverForMatrixFromEigensExample, allows_broken_matrix_construction_when_checks_disabled)
{
  allows_broken_matrix_construction_when_checks_disabled();
}

TEST_F(EigenSolverForMatrixFromEigensExample, throws_on_inconsistent_given_options)
{
  throws_on_inconsistent_given_options();
}

TEST_F(EigenSolverForMatrixFromEigensExample, is_constructible)
{
  is_constructible();
}

TEST_F(EigenSolverForMatrixFromEigensExample, gives_correct_eigenvalues)
{
  gives_correct_eigenvalues({{"lapacke", "1e-15"}, {"eigen", "1e-14"}});
}

TEST_F(EigenSolverForMatrixFromEigensExample, gives_correct_eigenvalues_in_correct_order)
{
  gives_correct_eigenvalues_in_correct_order({{"lapacke", "1e-15"}, {"eigen", "1e-14"}});
}

TEST_F(EigenSolverForMatrixFromEigensExample, gives_correct_eigendecomposition)
{
  gives_correct_eigendecomposition(1e-14);
}
