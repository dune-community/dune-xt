// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2019)
//   Tobias Leibner  (2020)

#include "config.h"

#include <cmath>

// Due to problems with old cblas versions, we only use the MKL cblas for the moment. See also
// https://github.com/dune-community/dune-xt/pull/198.
#if HAVE_MKL
#  include <mkl.h>
#endif

#include <dune/xt/common/exceptions.hh>

#include "cblas.hh"

#if HAVE_MKL
#  define DXTC_CBLAS_ONLY
#else
#  define DXTC_CBLAS_ONLY [[maybe_unused]]
#endif

namespace Dune::XT::Common::Cblas {


/**
 * \brief If true, calling any of the other methods makes sense.
 */
bool available()
{
#if HAVE_MKL
  return true;
#else
  return false;
#endif
}


int row_major()
{
#if HAVE_MKL
  return CblasRowMajor;
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing CBLAS or the intel mkl, check available() first!");
  return 1;
#endif
}


int col_major()
{
#if HAVE_MKL
  return CblasColMajor;
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing CBLAS or the intel mkl, check available() first!");
  return 1;
#endif
}


int left()
{
#if HAVE_MKL
  return CblasLeft;
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing CBLAS or the intel mkl, check available() first!");
  return 1;
#endif
}


int right()
{
#if HAVE_MKL
  return CblasRight;
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing CBLAS or the intel mkl, check available() first!");
  return 1;
#endif
}


int upper()
{
#if HAVE_MKL
  return CblasUpper;
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing CBLAS or the intel mkl, check available() first!");
  return 1;
#endif
}


int lower()
{
#if HAVE_MKL
  return CblasLower;
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing CBLAS or the intel mkl, check available() first!");
  return 1;
#endif
}


int trans()
{
#if HAVE_MKL
  return CblasTrans;
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing CBLAS or the intel mkl, check available() first!");
  return 1;
#endif
}


int no_trans()
{
#if HAVE_MKL
  return CblasNoTrans;
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing CBLAS or the intel mkl, check available() first!");
  return 1;
#endif
}


int unit()
{
#if HAVE_MKL
  return CblasUnit;
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing CBLAS or the intel mkl, check available() first!");
  return 1;
#endif
}


int non_unit()
{
#if HAVE_MKL
  return CblasNonUnit;
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing CBLAS or the intel mkl, check available() first!");
  return 1;
#endif
}


void dgemv(DXTC_CBLAS_ONLY const int layout,
           DXTC_CBLAS_ONLY const int trans,
           DXTC_CBLAS_ONLY const int m,
           DXTC_CBLAS_ONLY const int n,
           DXTC_CBLAS_ONLY const double alpha,
           DXTC_CBLAS_ONLY const double* a,
           DXTC_CBLAS_ONLY const int lda,
           DXTC_CBLAS_ONLY const double* x,
           DXTC_CBLAS_ONLY const int incx,
           DXTC_CBLAS_ONLY const double beta,
           DXTC_CBLAS_ONLY double* y,
           DXTC_CBLAS_ONLY const int incy)
{
#if HAVE_MKL
  cblas_dgemv(static_cast<CBLAS_LAYOUT>(layout),
              static_cast<CBLAS_TRANSPOSE>(trans),
              m,
              n,
              alpha,
              a,
              lda,
              x,
              incx,
              beta,
              y,
              incy);
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing CBLAS or the intel mkl, check available() first!");
#endif
}


void dtrsm(DXTC_CBLAS_ONLY const int layout,
           DXTC_CBLAS_ONLY const int side,
           DXTC_CBLAS_ONLY const int uplo,
           DXTC_CBLAS_ONLY const int transa,
           DXTC_CBLAS_ONLY const int diag,
           DXTC_CBLAS_ONLY const int m,
           DXTC_CBLAS_ONLY const int n,
           DXTC_CBLAS_ONLY const double alpha,
           DXTC_CBLAS_ONLY const double* a,
           DXTC_CBLAS_ONLY const int lda,
           DXTC_CBLAS_ONLY double* b,
           DXTC_CBLAS_ONLY const int ldb)
{
#if HAVE_MKL
  cblas_dtrsm(static_cast<CBLAS_LAYOUT>(layout),
              static_cast<CBLAS_SIDE>(side),
              static_cast<CBLAS_UPLO>(uplo),
              static_cast<CBLAS_TRANSPOSE>(transa),
              static_cast<CBLAS_DIAG>(diag),
              m,
              n,
              alpha,
              a,
              lda,
              b,
              ldb);
#  ifndef NDEBUG
  for (int ii = 0; ii < m; ++ii)
    if (std::isnan(b[ii]) || std::isinf(b[ii]))
      DUNE_THROW(Dune::MathError, "Triangular solve using cblas_dtrsm failed!");
#  endif
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing CBLAS or the intel mkl, check available() first!");
#endif
}


void dtrsv(DXTC_CBLAS_ONLY const int layout,
           DXTC_CBLAS_ONLY const int uplo,
           DXTC_CBLAS_ONLY const int transa,
           DXTC_CBLAS_ONLY const int diag,
           DXTC_CBLAS_ONLY const int n,
           DXTC_CBLAS_ONLY const double* a,
           DXTC_CBLAS_ONLY const int lda,
           DXTC_CBLAS_ONLY double* x,
           DXTC_CBLAS_ONLY const int incx)
{
#if HAVE_MKL
  cblas_dtrsv(static_cast<CBLAS_LAYOUT>(layout),
              static_cast<CBLAS_UPLO>(uplo),
              static_cast<CBLAS_TRANSPOSE>(transa),
              static_cast<CBLAS_DIAG>(diag),
              n,
              a,
              lda,
              x,
              incx);
#  ifndef NDEBUG
  for (int ii = 0; ii < n; ++ii)
    if (std::isnan(x[ii]) || std::isinf(x[ii]))
      DUNE_THROW(Dune::MathError, "Triangular solve using cblas_dtrsv failed!");
#  endif
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing CBLAS or the intel mkl, check available() first!");
#endif
}


void ztrsm(DXTC_CBLAS_ONLY const int layout,
           DXTC_CBLAS_ONLY const int side,
           DXTC_CBLAS_ONLY const int uplo,
           DXTC_CBLAS_ONLY const int transa,
           DXTC_CBLAS_ONLY const int diag,
           DXTC_CBLAS_ONLY const int m,
           DXTC_CBLAS_ONLY const int n,
           DXTC_CBLAS_ONLY const void* alpha,
           DXTC_CBLAS_ONLY const void* a,
           DXTC_CBLAS_ONLY const int lda,
           DXTC_CBLAS_ONLY void* b,
           DXTC_CBLAS_ONLY const int ldb)
{
#if HAVE_MKL
  cblas_ztrsm(static_cast<CBLAS_LAYOUT>(layout),
              static_cast<CBLAS_SIDE>(side),
              static_cast<CBLAS_UPLO>(uplo),
              static_cast<CBLAS_TRANSPOSE>(transa),
              static_cast<CBLAS_DIAG>(diag),
              m,
              n,
              alpha,
              a,
              lda,
              b,
              ldb);
#  ifndef NDEBUG
  for (int ii = 0; ii < m; ++ii)
    if (std::isnan(std::abs(static_cast<std::complex<double>*>(b)[ii]))
        || std::isinf(std::abs(static_cast<std::complex<double>*>(b)[ii])))
      DUNE_THROW(Dune::MathError, "Triangular solve using cblas_ztrsm failed!");
#  endif
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing CBLAS or the intel mkl, check available() first!");
#endif
}


void ztrsv(DXTC_CBLAS_ONLY const int layout,
           DXTC_CBLAS_ONLY const int uplo,
           DXTC_CBLAS_ONLY const int transa,
           DXTC_CBLAS_ONLY const int diag,
           DXTC_CBLAS_ONLY const int n,
           DXTC_CBLAS_ONLY const void* a,
           DXTC_CBLAS_ONLY const int lda,
           DXTC_CBLAS_ONLY void* x,
           DXTC_CBLAS_ONLY const int incx)
{
#if HAVE_MKL
  cblas_ztrsv(static_cast<CBLAS_LAYOUT>(layout),
              static_cast<CBLAS_UPLO>(uplo),
              static_cast<CBLAS_TRANSPOSE>(transa),
              static_cast<CBLAS_DIAG>(diag),
              n,
              a,
              lda,
              x,
              incx);
#  ifndef NDEBUG
  for (int ii = 0; ii < n; ++ii)
    if (std::isnan(std::abs(static_cast<std::complex<double>*>(x)[ii]))
        || std::isinf(std::abs(static_cast<std::complex<double>*>(x)[ii])))
      DUNE_THROW(Dune::MathError, "Triangular solve using cblas_ztrsv failed!");
#  endif
#else
  DUNE_THROW(Exceptions::dependency_missing, "You are missing CBLAS or the intel mkl, check available() first!");
#endif
}


} // namespace Dune::XT::Common::Cblas
