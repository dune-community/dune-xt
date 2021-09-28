// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2020)
//   René Fritze     (2018 - 2020)
//   Tim Keil        (2020)
//   Tobias Leibner  (2018, 2020)

#include "config.h"

#include <cmath>
// without the following lapacke will include <complex.h>, which will break dune/commontypetraits.hh^^
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

#if HAVE_MKL
#  include <mkl.h>
#elif HAVE_LAPACKE
#  include <lapacke.h>
#endif

#include <dune/xt/common/exceptions.hh>

#include "lapacke.hh"

#if HAVE_MKL || HAVE_LAPACKE
#  define DXTC_LAPACKE_ONLY
#else
#  define DXTC_LAPACKE_ONLY [[maybe_unused]]
#endif

namespace Dune::XT::Common::Lapacke {


bool available()
{
#if HAVE_MKL || HAVE_LAPACKE
  return true;
#else
  return false;
#endif
}


int row_major()
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACK_ROW_MAJOR;
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int col_major()
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACK_COL_MAJOR;
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dgeev(DXTC_LAPACKE_ONLY int matrix_layout,
          DXTC_LAPACKE_ONLY char jobvl,
          DXTC_LAPACKE_ONLY char jobvr,
          DXTC_LAPACKE_ONLY int n,
          DXTC_LAPACKE_ONLY double* a,
          DXTC_LAPACKE_ONLY int lda,
          DXTC_LAPACKE_ONLY double* wr,
          DXTC_LAPACKE_ONLY double* wi,
          DXTC_LAPACKE_ONLY double* vl,
          DXTC_LAPACKE_ONLY int ldvl,
          DXTC_LAPACKE_ONLY double* vr,
          DXTC_LAPACKE_ONLY int ldvr)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dgeev(matrix_layout, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}
int dgeev_work(DXTC_LAPACKE_ONLY int matrix_layout,
               DXTC_LAPACKE_ONLY char jobvl,
               DXTC_LAPACKE_ONLY char jobvr,
               DXTC_LAPACKE_ONLY int n,
               DXTC_LAPACKE_ONLY double* a,
               DXTC_LAPACKE_ONLY int lda,
               DXTC_LAPACKE_ONLY double* wr,
               DXTC_LAPACKE_ONLY double* wi,
               DXTC_LAPACKE_ONLY double* vl,
               DXTC_LAPACKE_ONLY int ldvl,
               DXTC_LAPACKE_ONLY double* vr,
               DXTC_LAPACKE_ONLY int ldvr,
               DXTC_LAPACKE_ONLY double* work,
               DXTC_LAPACKE_ONLY int lwork)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dgeev_work(matrix_layout, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dgeevx(DXTC_LAPACKE_ONLY int matrix_layout,
           DXTC_LAPACKE_ONLY char balanc,
           DXTC_LAPACKE_ONLY char jobvl,
           DXTC_LAPACKE_ONLY char jobvr,
           DXTC_LAPACKE_ONLY char sense,
           DXTC_LAPACKE_ONLY int n,
           DXTC_LAPACKE_ONLY double* a,
           DXTC_LAPACKE_ONLY int lda,
           DXTC_LAPACKE_ONLY double* wr,
           DXTC_LAPACKE_ONLY double* wi,
           DXTC_LAPACKE_ONLY double* vl,
           DXTC_LAPACKE_ONLY int ldvl,
           DXTC_LAPACKE_ONLY double* vr,
           DXTC_LAPACKE_ONLY int ldvr,
           DXTC_LAPACKE_ONLY int* ilo,
           DXTC_LAPACKE_ONLY int* ihi,
           DXTC_LAPACKE_ONLY double* scale,
           DXTC_LAPACKE_ONLY double* abnrm,
           DXTC_LAPACKE_ONLY double* rconde,
           DXTC_LAPACKE_ONLY double* rcondv)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dgeevx(matrix_layout,
                        balanc,
                        jobvl,
                        jobvr,
                        sense,
                        n,
                        a,
                        lda,
                        wr,
                        wi,
                        vl,
                        ldvl,
                        vr,
                        ldvr,
                        ilo,
                        ihi,
                        scale,
                        abnrm,
                        rconde,
                        rcondv);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}

int dgeevx_work(DXTC_LAPACKE_ONLY int matrix_layout,
                DXTC_LAPACKE_ONLY char balanc,
                DXTC_LAPACKE_ONLY char jobvl,
                DXTC_LAPACKE_ONLY char jobvr,
                DXTC_LAPACKE_ONLY char sense,
                DXTC_LAPACKE_ONLY int n,
                DXTC_LAPACKE_ONLY double* a,
                DXTC_LAPACKE_ONLY int lda,
                DXTC_LAPACKE_ONLY double* wr,
                DXTC_LAPACKE_ONLY double* wi,
                DXTC_LAPACKE_ONLY double* vl,
                DXTC_LAPACKE_ONLY int ldvl,
                DXTC_LAPACKE_ONLY double* vr,
                DXTC_LAPACKE_ONLY int ldvr,
                DXTC_LAPACKE_ONLY int* ilo,
                DXTC_LAPACKE_ONLY int* ihi,
                DXTC_LAPACKE_ONLY double* scale,
                DXTC_LAPACKE_ONLY double* abnrm,
                DXTC_LAPACKE_ONLY double* rconde,
                DXTC_LAPACKE_ONLY double* rcondv,
                DXTC_LAPACKE_ONLY double* work,
                DXTC_LAPACKE_ONLY int lwork,
                DXTC_LAPACKE_ONLY int* iwork)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dgeevx_work(matrix_layout,
                             balanc,
                             jobvl,
                             jobvr,
                             sense,
                             n,
                             a,
                             lda,
                             wr,
                             wi,
                             vl,
                             ldvl,
                             vr,
                             ldvr,
                             ilo,
                             ihi,
                             scale,
                             abnrm,
                             rconde,
                             rcondv,
                             work,
                             lwork,
                             iwork);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dgeqp3(DXTC_LAPACKE_ONLY int matrix_layout,
           DXTC_LAPACKE_ONLY int m,
           DXTC_LAPACKE_ONLY int n,
           DXTC_LAPACKE_ONLY double* a,
           DXTC_LAPACKE_ONLY int lda,
           DXTC_LAPACKE_ONLY int* jpvt,
           DXTC_LAPACKE_ONLY double* tau)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dgeqp3(matrix_layout, m, n, a, lda, jpvt, tau);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}

int dgeqp3_work(DXTC_LAPACKE_ONLY int matrix_layout,
                DXTC_LAPACKE_ONLY int m,
                DXTC_LAPACKE_ONLY int n,
                DXTC_LAPACKE_ONLY double* a,
                DXTC_LAPACKE_ONLY int lda,
                DXTC_LAPACKE_ONLY int* jpvt,
                DXTC_LAPACKE_ONLY double* tau,
                DXTC_LAPACKE_ONLY double* work,
                DXTC_LAPACKE_ONLY int lwork)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dgeqp3_work(matrix_layout, m, n, a, lda, jpvt, tau, work, lwork);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}

int dgesvd(DXTC_LAPACKE_ONLY int matrix_layout,
           DXTC_LAPACKE_ONLY char jobu,
           DXTC_LAPACKE_ONLY char jobvt,
           DXTC_LAPACKE_ONLY int m,
           DXTC_LAPACKE_ONLY int n,
           DXTC_LAPACKE_ONLY double* a,
           DXTC_LAPACKE_ONLY int lda,
           DXTC_LAPACKE_ONLY double* s,
           DXTC_LAPACKE_ONLY double* u,
           DXTC_LAPACKE_ONLY int ldu,
           DXTC_LAPACKE_ONLY double* vt,
           DXTC_LAPACKE_ONLY int ldvt,
           DXTC_LAPACKE_ONLY double* superb)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}

double dlamch(DXTC_LAPACKE_ONLY char cmach)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dlamch(cmach);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1.;
#endif
}

int dorgqr(DXTC_LAPACKE_ONLY int matrix_layout,
           DXTC_LAPACKE_ONLY int m,
           DXTC_LAPACKE_ONLY int n,
           DXTC_LAPACKE_ONLY int k,
           DXTC_LAPACKE_ONLY double* a,
           DXTC_LAPACKE_ONLY int lda,
           DXTC_LAPACKE_ONLY const double* tau)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dorgqr(matrix_layout, m, n, k, a, lda, tau);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}

int dorgqr_work(DXTC_LAPACKE_ONLY int matrix_layout,
                DXTC_LAPACKE_ONLY int m,
                DXTC_LAPACKE_ONLY int n,
                DXTC_LAPACKE_ONLY int k,
                DXTC_LAPACKE_ONLY double* a,
                DXTC_LAPACKE_ONLY int lda,
                DXTC_LAPACKE_ONLY const double* tau,
                DXTC_LAPACKE_ONLY double* work,
                DXTC_LAPACKE_ONLY int lwork)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dorgqr_work(matrix_layout, m, n, k, a, lda, tau, work, lwork);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dormqr(DXTC_LAPACKE_ONLY int matrix_layout,
           DXTC_LAPACKE_ONLY char side,
           DXTC_LAPACKE_ONLY char trans,
           DXTC_LAPACKE_ONLY int m,
           DXTC_LAPACKE_ONLY int n,
           DXTC_LAPACKE_ONLY int k,
           DXTC_LAPACKE_ONLY const double* a,
           DXTC_LAPACKE_ONLY int lda,
           DXTC_LAPACKE_ONLY const double* tau,
           DXTC_LAPACKE_ONLY double* c,
           DXTC_LAPACKE_ONLY int ldc)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dormqr(matrix_layout, side, trans, m, n, k, a, lda, tau, c, ldc);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}

int dormqr_work(DXTC_LAPACKE_ONLY int matrix_layout,
                DXTC_LAPACKE_ONLY char side,
                DXTC_LAPACKE_ONLY char trans,
                DXTC_LAPACKE_ONLY int m,
                DXTC_LAPACKE_ONLY int n,
                DXTC_LAPACKE_ONLY int k,
                DXTC_LAPACKE_ONLY const double* a,
                DXTC_LAPACKE_ONLY int lda,
                DXTC_LAPACKE_ONLY const double* tau,
                DXTC_LAPACKE_ONLY double* c,
                DXTC_LAPACKE_ONLY int ldc,
                DXTC_LAPACKE_ONLY double* work,
                DXTC_LAPACKE_ONLY int lwork)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dormqr_work(matrix_layout, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dpotrf(DXTC_LAPACKE_ONLY int matrix_layout,
           DXTC_LAPACKE_ONLY char uplo,
           DXTC_LAPACKE_ONLY int n,
           DXTC_LAPACKE_ONLY double* a,
           DXTC_LAPACKE_ONLY int lda)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dpotrf(matrix_layout, uplo, n, a, lda);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}

int dpotrf_work(DXTC_LAPACKE_ONLY int matrix_layout,
                DXTC_LAPACKE_ONLY char uplo,
                DXTC_LAPACKE_ONLY int n,
                DXTC_LAPACKE_ONLY double* a,
                DXTC_LAPACKE_ONLY int lda)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dpotrf_work(matrix_layout, uplo, n, a, lda);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dptcon(DXTC_LAPACKE_ONLY int n,
           DXTC_LAPACKE_ONLY const double* d,
           DXTC_LAPACKE_ONLY const double* e,
           DXTC_LAPACKE_ONLY double anorm,
           DXTC_LAPACKE_ONLY double* rcond)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dptcon(n, d, e, anorm, rcond);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dpocon(DXTC_LAPACKE_ONLY int matrix_layout,
           DXTC_LAPACKE_ONLY char uplo,
           DXTC_LAPACKE_ONLY int n,
           DXTC_LAPACKE_ONLY const double* a,
           DXTC_LAPACKE_ONLY int lda,
           DXTC_LAPACKE_ONLY double anorm,
           DXTC_LAPACKE_ONLY double* rcond)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dpocon(matrix_layout, uplo, n, a, lda, anorm, rcond);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dsygv(DXTC_LAPACKE_ONLY int matrix_layout,
          DXTC_LAPACKE_ONLY int itype,
          DXTC_LAPACKE_ONLY char jobz,
          DXTC_LAPACKE_ONLY char uplo,
          DXTC_LAPACKE_ONLY int n,
          DXTC_LAPACKE_ONLY double* a,
          DXTC_LAPACKE_ONLY int lda,
          DXTC_LAPACKE_ONLY double* b,
          DXTC_LAPACKE_ONLY int ldb,
          DXTC_LAPACKE_ONLY double* w)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dsygv(matrix_layout, itype, jobz, uplo, n, a, lda, b, ldb, w);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dtrcon(DXTC_LAPACKE_ONLY int matrix_layout,
           DXTC_LAPACKE_ONLY char norm,
           DXTC_LAPACKE_ONLY char uplo,
           DXTC_LAPACKE_ONLY char diag,
           DXTC_LAPACKE_ONLY int n,
           DXTC_LAPACKE_ONLY const double* a,
           DXTC_LAPACKE_ONLY int lda,
           DXTC_LAPACKE_ONLY double* rcond)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dtrcon(matrix_layout, norm, uplo, diag, n, a, lda, rcond);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dpttrf(DXTC_LAPACKE_ONLY int n, DXTC_LAPACKE_ONLY double* d, DXTC_LAPACKE_ONLY double* e)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dpttrf(n, d, e);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dpttrs(DXTC_LAPACKE_ONLY int matrix_layout,
           DXTC_LAPACKE_ONLY int n,
           DXTC_LAPACKE_ONLY int nrhs,
           DXTC_LAPACKE_ONLY const double* d,
           DXTC_LAPACKE_ONLY const double* e,
           DXTC_LAPACKE_ONLY double* b,
           DXTC_LAPACKE_ONLY int ldb)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dpttrs(matrix_layout, n, nrhs, d, e, b, ldb);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int zgeqp3(DXTC_LAPACKE_ONLY int matrix_layout,
           DXTC_LAPACKE_ONLY int m,
           DXTC_LAPACKE_ONLY int n,
           DXTC_LAPACKE_ONLY std::complex<double>* a,
           DXTC_LAPACKE_ONLY int lda,
           DXTC_LAPACKE_ONLY int* jpvt,
           DXTC_LAPACKE_ONLY std::complex<double>* tau)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_zgeqp3(matrix_layout, m, n, a, lda, jpvt, tau);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int zungqr(DXTC_LAPACKE_ONLY int matrix_layout,
           DXTC_LAPACKE_ONLY int m,
           DXTC_LAPACKE_ONLY int n,
           DXTC_LAPACKE_ONLY int k,
           DXTC_LAPACKE_ONLY std::complex<double>* a,
           DXTC_LAPACKE_ONLY int lda,
           DXTC_LAPACKE_ONLY const std::complex<double>* tau)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_zungqr(matrix_layout, m, n, k, a, lda, tau);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int zunmqr(DXTC_LAPACKE_ONLY int matrix_layout,
           DXTC_LAPACKE_ONLY char side,
           DXTC_LAPACKE_ONLY char trans,
           DXTC_LAPACKE_ONLY int m,
           DXTC_LAPACKE_ONLY int n,
           DXTC_LAPACKE_ONLY int k,
           DXTC_LAPACKE_ONLY const std::complex<double>* a,
           DXTC_LAPACKE_ONLY int lda,
           DXTC_LAPACKE_ONLY const std::complex<double>* tau,
           DXTC_LAPACKE_ONLY std::complex<double>* c,
           DXTC_LAPACKE_ONLY int ldc)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_zunmqr(matrix_layout, side, trans, m, n, k, a, lda, tau, c, ldc);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


} // namespace Dune::XT::Common::Lapacke
