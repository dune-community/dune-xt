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


/**
 * \sa http://eigen.tuxfamily.org/dox/classEigen_1_1EigenSolver.html
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

{% for T_NAME, TESTMATRIXTYPE, TESTFIELDTYPE, TESTCOMPLEXMATRIXTYPE, TESTREALMATRIXTYPE in config.testtypes %}
struct EigenSolverForMatrixFromEigensExample_{{T_NAME}}
: public EigenSolverTestForMatricesWithRealEigenvaluesAndVectors<{{TESTMATRIXTYPE}},
        {{TESTFIELDTYPE}},
        {{TESTCOMPLEXMATRIXTYPE}},
        {{TESTREALMATRIXTYPE}}>
{
  using BaseType = EigenSolverTest;
  using typename BaseType::MatrixType;
  using typename BaseType::ComplexMatrixType;
  using typename BaseType::EigenValuesType;

  EigenSolverForMatrixFromEigensExample_{{T_NAME}}()
  {
    // these are the numbers in std::setprecision(17) from the above example
    matrix_ = XT::Common::from_string<MatrixType>( // clang-format off
                "[ 0.68037543430941905 -0.32955448857022196  -0.27043105441631332  -0.7167948892883933  -0.68664182149183084  0.025864788808796835; "
                " -0.21123414636181392  0.53645918962380801   0.026801820391231024  0.21393775251411729 -0.198111211507633    0.67822446938521441; "
                "  0.56619844751721171 -0.44445057839362445   0.90445945034942565  -0.96739885675134085 -0.74041910643708853  0.22527966519132248; "
                "  0.59688006695214657  0.10793991159086103   0.8323901360074013   -0.51422645874052608 -0.78238239594846148 -0.40793676460531392; "
                "  0.82329471587356862 -0.045205896275679502  0.27142345591980188  -0.72553684642796257  0.99784903600711794  0.27510453540603841; "
                " -0.60489726141323208  0.25774184952384882   0.43459385886536617   0.60835350845398084 -0.56348618937818618  0.048574380133568518]");
    expected_eigenvalues_ = XT::Common::from_string<EigenValuesType>("[0.048963596941571319+1.0595752114710908i "
                                                                     " 0.048963596941571319-1.0595752114710908i "
                                                                     " 0.96689928885360044 "
                                                                     " 0.35304325510300372 "
                                                                     " 0.61781064692153476+0.12941286410485134i "
                                                                     " 0.61781064692153476-0.12941286410485134i]");
    expected_eigenvectors_ = XT::Common::from_string<ComplexMatrixType>(
        "[-0.2916493853334598-0.45353428327568535i   -0.2916493853334598+0.45353428327568535i   -0.060739856420097799 -0.732892086787302      0.58989693747516225-0.12134050155506512i   0.58989693747516225+0.12134050155506512i; "
        "  0.13427786119296695-0.10415042620914604i   0.13427786119296695+0.10415042620914604i  -0.79936618427058692   0.13588801812778584    0.33446578429890217+0.36779173930086734i   0.33446578429890217-0.36779173930086734i; "
        " -0.4224933428958097-0.17959111262899077i   -0.4224933428958097+0.17959111262899077i    0.19181897363210443   0.056319404917871921  -0.33512506847695012-0.1433533445630015i   -0.33512506847695012+0.1433533445630015i; "
        " -0.58861043190958418+0.027350002433760786i -0.58861043190958418-0.027350002433760786i -0.078824737066622547 -0.62743714713230025    0.32187174947115971-0.1554720757387448i    0.32187174947115971+0.1554720757387448i; "
        " -0.24767280752759105+0.13196040251638746i  -0.24767280752759105-0.13196040251638746i   0.40098186135709829   0.21799673455433938   -0.33488773557308854-0.076133571019896026i -0.33488773557308854+0.076133571019896026i; "
        "  0.10538083226488951+0.18040258229827225i   0.10538083226488951-0.18040258229827225i  -0.39182912389260799  -0.0056400714249183401 -0.032444861799488058+0.10261186083116479i -0.032444861799488058-0.10261186083116479i]"); // clang-format on
    all_matrices_and_expected_eigenvalues_and_vectors_are_computed_ = true;
  }

  using BaseType::all_matrices_and_expected_eigenvalues_and_vectors_are_computed_;
  using BaseType::matrix_;
  using BaseType::expected_eigenvalues_;
  using BaseType::expected_eigenvectors_;
}; // struct EigenSolverForMatrixFromEigensExample_{{T_NAME}}


TEST_F(EigenSolverForMatrixFromEigensExample_{{T_NAME}}, exports_correct_types)
{
  exports_correct_types();
}

TEST_F(EigenSolverForMatrixFromEigensExample_{{T_NAME}}, has_types_and_options)
{
  has_types_and_options();
}

TEST_F(EigenSolverForMatrixFromEigensExample_{{T_NAME}}, throws_on_broken_matrix_construction)
{
  throws_on_broken_matrix_construction();
}

TEST_F(EigenSolverForMatrixFromEigensExample_{{T_NAME}}, allows_broken_matrix_construction_when_checks_disabled)
{
  allows_broken_matrix_construction_when_checks_disabled();
}

TEST_F(EigenSolverForMatrixFromEigensExample_{{T_NAME}}, throws_on_inconsistent_given_options)
{
  throws_on_inconsistent_given_options();
}

TEST_F(EigenSolverForMatrixFromEigensExample_{{T_NAME}}, is_constructible)
{
  is_constructible();
}

TEST_F(EigenSolverForMatrixFromEigensExample_{{T_NAME}}, gives_correct_eigenvalues)
{
  gives_correct_eigenvalues({ {"lapack", "1e-14"}, {"eigen", "1e-15"}, {"numpy", "1e-14"}, {"shifted_qr", /*we_expect_a_failure: */ "-1"} });
}

TEST_F(EigenSolverForMatrixFromEigensExample_{{T_NAME}}, gives_correct_eigenvalues_in_correct_order)
{
  gives_correct_eigenvalues_in_correct_order({ {"lapack", "1e-14"}, {"eigen", "1e-15"}, {"numpy", "1e-14"}, {"shifted_qr", /*we_expect_a_failure: */ "-1"} });
}

TEST_F(EigenSolverForMatrixFromEigensExample_{{T_NAME}}, gives_correct_eigenvectors_in_correct_order)
{
  gives_correct_eigenvectors_in_correct_order(
      { {"lapack", /*we expect a failure: */ "-1"}, {"eigen", "1e-15"}, {"numpy", /*we expect a failure: */ "-1"}, {"shifted_qr", /*we_expect_a_failure: */ "-1"} });
}

TEST_F(EigenSolverForMatrixFromEigensExample_{{T_NAME}}, gives_correct_eigendecomposition)
{
  gives_correct_eigendecomposition({ {"lapack", "1e-14"}, {"eigen", "1e-14"}, {"numpy", "1e-14"}, {"shifted_qr", /*we_expect_a_failure: */ "-1"} });
}

{% endfor %}
