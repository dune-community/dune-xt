// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2019)

#include "config.h"

#include <cmath>

// Due to probles with old cblas versions, we only use the MKL cblas for the moment. See also
// https://github.com/dune-community/dune-xt/pull/198.
#if HAVE_MKL
#  include <mkl.h>
#endif

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/unused.hh>

#include "cblas.hh"

#if HAVE_MKL
#  define DXTC_CBLAS_ONLY(param) param
#else
#  define DXTC_CBLAS_ONLY(param) DXTC_UNUSED(param)
#endif

namespace Dune {
namespace XT {
namespace Common {
namespace Cblas {


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


void dgemv(const int DXTC_CBLAS_ONLY(layout),
           const int DXTC_CBLAS_ONLY(trans),
           const int DXTC_CBLAS_ONLY(m),
           const int DXTC_CBLAS_ONLY(n),
           const double DXTC_CBLAS_ONLY(alpha),
           const double* DXTC_CBLAS_ONLY(a),
           const int DXTC_CBLAS_ONLY(lda),
           const double* DXTC_CBLAS_ONLY(x),
           const int DXTC_CBLAS_ONLY(incx),
           const double DXTC_CBLAS_ONLY(beta),
           double* DXTC_CBLAS_ONLY(y),
           const int DXTC_CBLAS_ONLY(incy))
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


void dtrsm(const int DXTC_CBLAS_ONLY(layout),
           const int DXTC_CBLAS_ONLY(side),
           const int DXTC_CBLAS_ONLY(uplo),
           const int DXTC_CBLAS_ONLY(transa),
           const int DXTC_CBLAS_ONLY(diag),
           const int DXTC_CBLAS_ONLY(m),
           const int DXTC_CBLAS_ONLY(n),
           const double DXTC_CBLAS_ONLY(alpha),
           const double* DXTC_CBLAS_ONLY(a),
           const int DXTC_CBLAS_ONLY(lda),
           double* DXTC_CBLAS_ONLY(b),
           const int DXTC_CBLAS_ONLY(ldb))
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


void dtrsv(const int DXTC_CBLAS_ONLY(layout),
           const int DXTC_CBLAS_ONLY(uplo),
           const int DXTC_CBLAS_ONLY(transa),
           const int DXTC_CBLAS_ONLY(diag),
           const int DXTC_CBLAS_ONLY(n),
           const double* DXTC_CBLAS_ONLY(a),
           const int DXTC_CBLAS_ONLY(lda),
           double* DXTC_CBLAS_ONLY(x),
           const int DXTC_CBLAS_ONLY(incx))
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


void ztrsm(const int DXTC_CBLAS_ONLY(layout),
           const int DXTC_CBLAS_ONLY(side),
           const int DXTC_CBLAS_ONLY(uplo),
           const int DXTC_CBLAS_ONLY(transa),
           const int DXTC_CBLAS_ONLY(diag),
           const int DXTC_CBLAS_ONLY(m),
           const int DXTC_CBLAS_ONLY(n),
           const void* DXTC_CBLAS_ONLY(alpha),
           const void* DXTC_CBLAS_ONLY(a),
           const int DXTC_CBLAS_ONLY(lda),
           void* DXTC_CBLAS_ONLY(b),
           const int DXTC_CBLAS_ONLY(ldb))
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


void ztrsv(const int DXTC_CBLAS_ONLY(layout),
           const int DXTC_CBLAS_ONLY(uplo),
           const int DXTC_CBLAS_ONLY(transa),
           const int DXTC_CBLAS_ONLY(diag),
           const int DXTC_CBLAS_ONLY(n),
           const void* DXTC_CBLAS_ONLY(a),
           const int DXTC_CBLAS_ONLY(lda),
           void* DXTC_CBLAS_ONLY(x),
           const int DXTC_CBLAS_ONLY(incx))
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


} // namespace Cblas
} // namespace Common
} // namespace XT
} // namespace Dune
