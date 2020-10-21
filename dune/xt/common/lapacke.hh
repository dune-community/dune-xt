// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2014, 2016 - 2019)
//   René Fritze     (2011 - 2012, 2015 - 2016, 2018 - 2019)
//   Tobias Leibner  (2017 - 2018, 2020)

#ifndef DUNE_XT_COMMON_LAPACKE_HH
#define DUNE_XT_COMMON_LAPACKE_HH

#include <complex>

namespace Dune::XT::Common::Lapacke {


/**
 * \brief If true, calling any of the other methods makes sense.
 */
bool available();


/**
 * \brief Wrapper around LAPACK_ROW_MAJOR
 * \sa    LAPACK_ROW_MAJOR
 */
int row_major();


/**
 * \brief Wrapper around LAPACK_COL_MAJOR
 * \sa    LAPACK_COL_MAJOR
 */
int col_major();


/**
 * \brief Wrapper around LAPACKE_dgeev
 * \sa    LAPACKE_dgeev
 */
int dgeev(int matrix_layout,
          char jobvl,
          char jobvr,
          int n,
          double* a,
          int lda,
          double* wr,
          double* wi,
          double* vl,
          int ldvl,
          double* vr,
          int ldvr);

/**
 * \brief Wrapper around LAPACKE_dgeev_work
 * \sa    LAPACKE_dgeev_work
 */
int dgeev_work(int matrix_layout,
               char jobvl,
               char jobvr,
               int n,
               double* a,
               int lda,
               double* wr,
               double* wi,
               double* vl,
               int ldvl,
               double* vr,
               int ldvr,
               double* work,
               int lwork);


/**
 * \brief Wrapper around LAPACKE_dgeevx
 * \sa    LAPACKE_dgeevx
 */
int dgeevx(int matrix_layout,
           char balanc,
           char jobvl,
           char jobvr,
           char sense,
           int n,
           double* a,
           int lda,
           double* wr,
           double* wi,
           double* vl,
           int ldvl,
           double* vr,
           int ldvr,
           int* ilo,
           int* ihi,
           double* scale,
           double* abnrm,
           double* rconde,
           double* rcondv);

/**
 * \brief Wrapper around LAPACKE_dgeevx_work
 * When calling this routine with lwork = -1, the optimal size for the working array work is returned in work[0]. The
 * size of the integer working array iwork is 2*n-2.
 * \sa    LAPACKE_dgeevx_work
 */
int dgeevx_work(int matrix_layout,
                char balanc,
                char jobvl,
                char jobvr,
                char sense,
                int n,
                double* a,
                int lda,
                double* wr,
                double* wi,
                double* vl,
                int ldvl,
                double* vr,
                int ldvr,
                int* ilo,
                int* ihi,
                double* scale,
                double* abnrm,
                double* rconde,
                double* rcondv,
                double* work,
                int lwork,
                int* iwork);

/**
 * \brief Wrapper around LAPACKE_dgeqp3
 * \sa    LAPACKE_dgeqp3
 */
int dgeqp3(int matrix_layout, int m, int n, double* a, int lda, int* jpvt, double* tau);

/**
 * \brief Wrapper around LAPACKE_dgeqp3_work
 * \sa    LAPACKE_dgeqp3
 */
int dgeqp3_work(int matrix_layout, int m, int n, double* a, int lda, int* jpvt, double* tau, double* work, int lwork);

/**
 * \brief Wrapper around LAPACKE_dgesvd
 * \sa    LAPACKE_dgesvd
 */
int dgesvd(int matrix_layout,
           char jobu,
           char jobvt,
           int m,
           int n,
           double* a,
           int lda,
           double* s,
           double* u,
           int ldu,
           double* vt,
           int ldvt,
           double* superb);

/**
 * \brief Wrapper around LAPACKE_dlamch
 * \sa    LAPACKE_dlamch
 */
double dlamch(char cmach);

/**
 * \brief Wrapper around LAPACKE_dorgqr
 * \sa    LAPACKE_dorgqr
 */
int dorgqr(int matrix_layout, int m, int n, int k, double* a, int lda, const double* tau);

/**
 * \brief Wrapper around LAPACKE_dorgqr_work
 * \sa    LAPACKE_dorgqr_work
 */
int dorgqr_work(int matrix_layout, int m, int n, int k, double* a, int lda, const double* tau, double* work, int lwork);


/**
 * \brief Wrapper around LAPACKE_dormqr
 * \sa    LAPACKE_dormqr
 */
int dormqr(int matrix_layout,
           char side,
           char trans,
           int m,
           int n,
           int k,
           const double* a,
           int lda,
           const double* tau,
           double* c,
           int ldc);

/**
 * \brief Wrapper around LAPACKE_dormqr_work
 * \sa    LAPACKE_dormqr_work
 */
int dormqr_work(int matrix_layout,
                char side,
                char trans,
                int m,
                int n,
                int k,
                const double* a,
                int lda,
                const double* tau,
                double* c,
                int ldc,
                double* work,
                int lwork);


/**
 * \brief Wrapper around LAPACKE_dpotrf
 * \sa    LAPACKE_dpotrf
 */
int dpotrf(int matrix_layout, char uplo, int n, double* a, int lda);

/**
 * \brief Wrapper around LAPACKE_dpotrf_work
 * \sa    LAPACKE_dpotrf_work
 */
int dpotrf_work(int matrix_layout, char uplo, int n, double* a, int lda);


/**
 * \brief Wrapper around LAPACKE_dptcon
 * \sa    LAPACKE_dptcon
 */
int dptcon(int n, const double* d, const double* e, double anorm, double* rcond);


/**
 * \brief Wrapper around LAPACKE_dptcon
 * \sa    LAPACKE_dpocon
 */
int dpocon(int matrix_layout, char uplo, int n, const double* a, int lda, double anorm, double* rcond);


/**
 * \brief Wrapper around LAPACKE_dsygv
 *
 * \sa    LAPACKE_dsygv
 * \sa    https://software.intel.com/en-us/mkl-developer-reference-c-sygv
 */
int dsygv(int matrix_layout, int itype, char jobz, char uplo, int n, double* a, int lda, double* b, int ldb, double* w);


/**
 * \brief Wrapper around LAPACKE_dtrcon
 * \sa    LAPACKE_dtrcon
 */
int dtrcon(int matrix_layout, char norm, char uplo, char diag, int n, const double* a, int lda, double* rcond);


/**
 * \brief Wrapper around LAPACKE_dpttrf
 * \sa    LAPACKE_dpttrf
 */
int dpttrf(int n, double* d, double* e);


/**
 * \brief Wrapper around LAPACKE_dpttrs
 * \sa    LAPACKE_dpttrs
 */
int dpttrs(int matrix_layout, int n, int nrhs, const double* d, const double* e, double* b, int ldb);


/**
 * \brief Wrapper around LAPACKE_zgeqp3
 * \sa    LAPACKE_zgeqp3
 */
int zgeqp3(int matrix_layout, int m, int n, std::complex<double>* a, int lda, int* jpvt, std::complex<double>* tau);


/**
 * \brief Wrapper around LAPACKE_zunmqr
 * \sa    LAPACKE_zunmqr
 */
int zunmqr(int matrix_layout,
           char side,
           char trans,
           int m,
           int n,
           int k,
           const std::complex<double>* a,
           int lda,
           const std::complex<double>* tau,
           std::complex<double>* c,
           int ldc);


/**
 * \brief Wrapper around LAPACKE_zungqr
 * \sa    LAPACKE_zungqr
 */
int zungqr(int matrix_layout, int m, int n, int k, std::complex<double>* a, int lda, const std::complex<double>* tau);


} // namespace Dune::XT::Common::Lapacke

#endif // DUNE_XT_COMMON_LAPACKE_HH
