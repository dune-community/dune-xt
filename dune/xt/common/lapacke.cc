// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2019)
//   René Fritze     (2018 - 2019)
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
#include <dune/xt/common/unused.hh>

#include "lapacke.hh"

#if HAVE_MKL || HAVE_LAPACKE
#  define DXTC_LAPACKE_ONLY(param) param
#else
#  define DXTC_LAPACKE_ONLY(param) DXTC_UNUSED(param)
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


int dgeev(int DXTC_LAPACKE_ONLY(matrix_layout),
          char DXTC_LAPACKE_ONLY(jobvl),
          char DXTC_LAPACKE_ONLY(jobvr),
          int DXTC_LAPACKE_ONLY(n),
          double* DXTC_LAPACKE_ONLY(a),
          int DXTC_LAPACKE_ONLY(lda),
          double* DXTC_LAPACKE_ONLY(wr),
          double* DXTC_LAPACKE_ONLY(wi),
          double* DXTC_LAPACKE_ONLY(vl),
          int DXTC_LAPACKE_ONLY(ldvl),
          double* DXTC_LAPACKE_ONLY(vr),
          int DXTC_LAPACKE_ONLY(ldvr))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dgeev(matrix_layout, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}
int dgeev_work(int DXTC_LAPACKE_ONLY(matrix_layout),
               char DXTC_LAPACKE_ONLY(jobvl),
               char DXTC_LAPACKE_ONLY(jobvr),
               int DXTC_LAPACKE_ONLY(n),
               double* DXTC_LAPACKE_ONLY(a),
               int DXTC_LAPACKE_ONLY(lda),
               double* DXTC_LAPACKE_ONLY(wr),
               double* DXTC_LAPACKE_ONLY(wi),
               double* DXTC_LAPACKE_ONLY(vl),
               int DXTC_LAPACKE_ONLY(ldvl),
               double* DXTC_LAPACKE_ONLY(vr),
               int DXTC_LAPACKE_ONLY(ldvr),
               double* DXTC_LAPACKE_ONLY(work),
               int DXTC_LAPACKE_ONLY(lwork))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dgeev_work(matrix_layout, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dgeevx(int DXTC_LAPACKE_ONLY(matrix_layout),
           char DXTC_LAPACKE_ONLY(balanc),
           char DXTC_LAPACKE_ONLY(jobvl),
           char DXTC_LAPACKE_ONLY(jobvr),
           char DXTC_LAPACKE_ONLY(sense),
           int DXTC_LAPACKE_ONLY(n),
           double* DXTC_LAPACKE_ONLY(a),
           int DXTC_LAPACKE_ONLY(lda),
           double* DXTC_LAPACKE_ONLY(wr),
           double* DXTC_LAPACKE_ONLY(wi),
           double* DXTC_LAPACKE_ONLY(vl),
           int DXTC_LAPACKE_ONLY(ldvl),
           double* DXTC_LAPACKE_ONLY(vr),
           int DXTC_LAPACKE_ONLY(ldvr),
           int* DXTC_LAPACKE_ONLY(ilo),
           int* DXTC_LAPACKE_ONLY(ihi),
           double* DXTC_LAPACKE_ONLY(scale),
           double* DXTC_LAPACKE_ONLY(abnrm),
           double* DXTC_LAPACKE_ONLY(rconde),
           double* DXTC_LAPACKE_ONLY(rcondv))
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

int dgeevx_work(int DXTC_LAPACKE_ONLY(matrix_layout),
                char DXTC_LAPACKE_ONLY(balanc),
                char DXTC_LAPACKE_ONLY(jobvl),
                char DXTC_LAPACKE_ONLY(jobvr),
                char DXTC_LAPACKE_ONLY(sense),
                int DXTC_LAPACKE_ONLY(n),
                double* DXTC_LAPACKE_ONLY(a),
                int DXTC_LAPACKE_ONLY(lda),
                double* DXTC_LAPACKE_ONLY(wr),
                double* DXTC_LAPACKE_ONLY(wi),
                double* DXTC_LAPACKE_ONLY(vl),
                int DXTC_LAPACKE_ONLY(ldvl),
                double* DXTC_LAPACKE_ONLY(vr),
                int DXTC_LAPACKE_ONLY(ldvr),
                int* DXTC_LAPACKE_ONLY(ilo),
                int* DXTC_LAPACKE_ONLY(ihi),
                double* DXTC_LAPACKE_ONLY(scale),
                double* DXTC_LAPACKE_ONLY(abnrm),
                double* DXTC_LAPACKE_ONLY(rconde),
                double* DXTC_LAPACKE_ONLY(rcondv),
                double* DXTC_LAPACKE_ONLY(work),
                int DXTC_LAPACKE_ONLY(lwork),
                int* DXTC_LAPACKE_ONLY(iwork))
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


int dgeqp3(int DXTC_LAPACKE_ONLY(matrix_layout),
           int DXTC_LAPACKE_ONLY(m),
           int DXTC_LAPACKE_ONLY(n),
           double* DXTC_LAPACKE_ONLY(a),
           int DXTC_LAPACKE_ONLY(lda),
           int* DXTC_LAPACKE_ONLY(jpvt),
           double* DXTC_LAPACKE_ONLY(tau))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dgeqp3(matrix_layout, m, n, a, lda, jpvt, tau);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}

int dgeqp3_work(int DXTC_LAPACKE_ONLY(matrix_layout),
                int DXTC_LAPACKE_ONLY(m),
                int DXTC_LAPACKE_ONLY(n),
                double* DXTC_LAPACKE_ONLY(a),
                int DXTC_LAPACKE_ONLY(lda),
                int* DXTC_LAPACKE_ONLY(jpvt),
                double* DXTC_LAPACKE_ONLY(tau),
                double* DXTC_LAPACKE_ONLY(work),
                int DXTC_LAPACKE_ONLY(lwork))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dgeqp3_work(matrix_layout, m, n, a, lda, jpvt, tau, work, lwork);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}

int dgesvd(int DXTC_LAPACKE_ONLY(matrix_layout),
           char DXTC_LAPACKE_ONLY(jobu),
           char DXTC_LAPACKE_ONLY(jobvt),
           int DXTC_LAPACKE_ONLY(m),
           int DXTC_LAPACKE_ONLY(n),
           double* DXTC_LAPACKE_ONLY(a),
           int DXTC_LAPACKE_ONLY(lda),
           double* DXTC_LAPACKE_ONLY(s),
           double* DXTC_LAPACKE_ONLY(u),
           int DXTC_LAPACKE_ONLY(ldu),
           double* DXTC_LAPACKE_ONLY(vt),
           int DXTC_LAPACKE_ONLY(ldvt),
           double* DXTC_LAPACKE_ONLY(superb))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}

double dlamch(char cmach)
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dlamch(cmach);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1.;
#endif
}

int dorgqr(int DXTC_LAPACKE_ONLY(matrix_layout),
           int DXTC_LAPACKE_ONLY(m),
           int DXTC_LAPACKE_ONLY(n),
           int DXTC_LAPACKE_ONLY(k),
           double* DXTC_LAPACKE_ONLY(a),
           int DXTC_LAPACKE_ONLY(lda),
           const double* DXTC_LAPACKE_ONLY(tau))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dorgqr(matrix_layout, m, n, k, a, lda, tau);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}

int dorgqr_work(int DXTC_LAPACKE_ONLY(matrix_layout),
                int DXTC_LAPACKE_ONLY(m),
                int DXTC_LAPACKE_ONLY(n),
                int DXTC_LAPACKE_ONLY(k),
                double* DXTC_LAPACKE_ONLY(a),
                int DXTC_LAPACKE_ONLY(lda),
                const double* DXTC_LAPACKE_ONLY(tau),
                double* DXTC_LAPACKE_ONLY(work),
                int DXTC_LAPACKE_ONLY(lwork))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dorgqr_work(matrix_layout, m, n, k, a, lda, tau, work, lwork);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dormqr(int DXTC_LAPACKE_ONLY(matrix_layout),
           char DXTC_LAPACKE_ONLY(side),
           char DXTC_LAPACKE_ONLY(trans),
           int DXTC_LAPACKE_ONLY(m),
           int DXTC_LAPACKE_ONLY(n),
           int DXTC_LAPACKE_ONLY(k),
           const double* DXTC_LAPACKE_ONLY(a),
           int DXTC_LAPACKE_ONLY(lda),
           const double* DXTC_LAPACKE_ONLY(tau),
           double* DXTC_LAPACKE_ONLY(c),
           int DXTC_LAPACKE_ONLY(ldc))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dormqr(matrix_layout, side, trans, m, n, k, a, lda, tau, c, ldc);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}

int dormqr_work(int DXTC_LAPACKE_ONLY(matrix_layout),
                char DXTC_LAPACKE_ONLY(side),
                char DXTC_LAPACKE_ONLY(trans),
                int DXTC_LAPACKE_ONLY(m),
                int DXTC_LAPACKE_ONLY(n),
                int DXTC_LAPACKE_ONLY(k),
                const double* DXTC_LAPACKE_ONLY(a),
                int DXTC_LAPACKE_ONLY(lda),
                const double* DXTC_LAPACKE_ONLY(tau),
                double* DXTC_LAPACKE_ONLY(c),
                int DXTC_LAPACKE_ONLY(ldc),
                double* DXTC_LAPACKE_ONLY(work),
                int DXTC_LAPACKE_ONLY(lwork))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dormqr_work(matrix_layout, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dpotrf(int DXTC_LAPACKE_ONLY(matrix_layout),
           char DXTC_LAPACKE_ONLY(uplo),
           int DXTC_LAPACKE_ONLY(n),
           double* DXTC_LAPACKE_ONLY(a),
           int DXTC_LAPACKE_ONLY(lda))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dpotrf(matrix_layout, uplo, n, a, lda);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}

int dpotrf_work(int DXTC_LAPACKE_ONLY(matrix_layout),
                char DXTC_LAPACKE_ONLY(uplo),
                int DXTC_LAPACKE_ONLY(n),
                double* DXTC_LAPACKE_ONLY(a),
                int DXTC_LAPACKE_ONLY(lda))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dpotrf_work(matrix_layout, uplo, n, a, lda);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dptcon(int DXTC_LAPACKE_ONLY(n),
           const double* DXTC_LAPACKE_ONLY(d),
           const double* DXTC_LAPACKE_ONLY(e),
           double DXTC_LAPACKE_ONLY(anorm),
           double* DXTC_LAPACKE_ONLY(rcond))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dptcon(n, d, e, anorm, rcond);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dpocon(int DXTC_LAPACKE_ONLY(matrix_layout),
           char DXTC_LAPACKE_ONLY(uplo),
           int DXTC_LAPACKE_ONLY(n),
           const double* DXTC_LAPACKE_ONLY(a),
           int DXTC_LAPACKE_ONLY(lda),
           double DXTC_LAPACKE_ONLY(anorm),
           double* DXTC_LAPACKE_ONLY(rcond))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dpocon(matrix_layout, uplo, n, a, lda, anorm, rcond);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dsygv(int DXTC_LAPACKE_ONLY(matrix_layout),
          int DXTC_LAPACKE_ONLY(itype),
          char DXTC_LAPACKE_ONLY(jobz),
          char DXTC_LAPACKE_ONLY(uplo),
          int DXTC_LAPACKE_ONLY(n),
          double* DXTC_LAPACKE_ONLY(a),
          int DXTC_LAPACKE_ONLY(lda),
          double* DXTC_LAPACKE_ONLY(b),
          int DXTC_LAPACKE_ONLY(ldb),
          double* DXTC_LAPACKE_ONLY(w))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dsygv(matrix_layout, itype, jobz, uplo, n, a, lda, b, ldb, w);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dtrcon(int DXTC_LAPACKE_ONLY(matrix_layout),
           char DXTC_LAPACKE_ONLY(norm),
           char DXTC_LAPACKE_ONLY(uplo),
           char DXTC_LAPACKE_ONLY(diag),
           int DXTC_LAPACKE_ONLY(n),
           const double* DXTC_LAPACKE_ONLY(a),
           int DXTC_LAPACKE_ONLY(lda),
           double* DXTC_LAPACKE_ONLY(rcond))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dtrcon(matrix_layout, norm, uplo, diag, n, a, lda, rcond);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dpttrf(int DXTC_LAPACKE_ONLY(n), double* DXTC_LAPACKE_ONLY(d), double* DXTC_LAPACKE_ONLY(e))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dpttrf(n, d, e);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int dpttrs(int DXTC_LAPACKE_ONLY(matrix_layout),
           int DXTC_LAPACKE_ONLY(n),
           int DXTC_LAPACKE_ONLY(nrhs),
           const double* DXTC_LAPACKE_ONLY(d),
           const double* DXTC_LAPACKE_ONLY(e),
           double* DXTC_LAPACKE_ONLY(b),
           int DXTC_LAPACKE_ONLY(ldb))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_dpttrs(matrix_layout, n, nrhs, d, e, b, ldb);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int zgeqp3(int DXTC_LAPACKE_ONLY(matrix_layout),
           int DXTC_LAPACKE_ONLY(m),
           int DXTC_LAPACKE_ONLY(n),
           std::complex<double>* DXTC_LAPACKE_ONLY(a),
           int DXTC_LAPACKE_ONLY(lda),
           int* DXTC_LAPACKE_ONLY(jpvt),
           std::complex<double>* DXTC_LAPACKE_ONLY(tau))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_zgeqp3(matrix_layout, m, n, a, lda, jpvt, tau);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int zungqr(int DXTC_LAPACKE_ONLY(matrix_layout),
           int DXTC_LAPACKE_ONLY(m),
           int DXTC_LAPACKE_ONLY(n),
           int DXTC_LAPACKE_ONLY(k),
           std::complex<double>* DXTC_LAPACKE_ONLY(a),
           int DXTC_LAPACKE_ONLY(lda),
           const std::complex<double>* DXTC_LAPACKE_ONLY(tau))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_zungqr(matrix_layout, m, n, k, a, lda, tau);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


int zunmqr(int DXTC_LAPACKE_ONLY(matrix_layout),
           char DXTC_LAPACKE_ONLY(side),
           char DXTC_LAPACKE_ONLY(trans),
           int DXTC_LAPACKE_ONLY(m),
           int DXTC_LAPACKE_ONLY(n),
           int DXTC_LAPACKE_ONLY(k),
           const std::complex<double>* DXTC_LAPACKE_ONLY(a),
           int DXTC_LAPACKE_ONLY(lda),
           const std::complex<double>* DXTC_LAPACKE_ONLY(tau),
           std::complex<double>* DXTC_LAPACKE_ONLY(c),
           int DXTC_LAPACKE_ONLY(ldc))
{
#if HAVE_MKL || HAVE_LAPACKE
  return LAPACKE_zunmqr(matrix_layout, side, trans, m, n, k, a, lda, tau, c, ldc);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing lapacke or the intel mkl, check available() first!");
  return 1;
#endif
}


} // namespace Dune::XT::Common::Lapacke
