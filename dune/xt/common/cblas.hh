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

#ifndef DUNE_XT_COMMON_BLAS_HH
#define DUNE_XT_COMMON_BLAS_HH

namespace Dune {
namespace XT {
namespace Common {
namespace Cblas {


/**
 * \brief If true, calling any of the other methods makes sense.
 */
bool available();


/**
 * \brief Wrapper around CBlasRowMajor
 * \sa    CblasRowMajor
 */
int row_major();


/**
 * \brief Wrapper around CBlasColMajor
 * \sa    CblasColMajor
 */
int col_major();


/**
 * \brief Wrapper around CBlasLeft
 * \sa    CblasLeft
 */
int left();


/**
 * \brief Wrapper around CBlasRight
 * \sa    CblasRight
 */
int right();


/**
 * \brief Wrapper around CblasUpper
 * \sa    CblasUpper
 */
int upper();


/**
 * \brief Wrapper around CblasLower
 * \sa    CblasLower
 */
int lower();


/**
 * \brief Wrapper around CblasTrans
 * \sa    CblasTrans
 */
int trans();


/**
 * \brief Wrapper around CblasNoTrans
 * \sa    CblasNoTrans
 */
int no_trans();


/**
 * \brief Wrapper around CblasUnit
 * \sa    CblasUnit
 */
int unit();


/**
 * \brief Wrapper around CblasNonUnit
 * \sa    CblasNonUnit
 */
int non_unit();


/**
 * \brief Wrapper around cblas_dgemv
 * \sa    cblas_dgemv
 */
void dgemv(const int layout,
           const int trans,
           const int m,
           const int n,
           const double alpha,
           const double* a,
           const int lda,
           const double* x,
           const int incx,
           const double beta,
           double* y,
           const int incy);


/**
 * \brief Wrapper around cblas_dtrsm
 * \sa    cblas_dtrsm
 */
void dtrsm(const int layout,
           const int side,
           const int uplo,
           const int transa,
           const int diag,
           const int m,
           const int n,
           const double alpha,
           const double* a,
           const int lda,
           double* b,
           const int ldb);


/**
 * \brief Wrapper around cblas_dtrsv
 * \sa    cblas_dtrsv
 */
void dtrsv(const int layout,
           const int uplo,
           const int transa,
           const int diag,
           const int n,
           const double* a,
           const int lda,
           double* x,
           const int incx);


/**
 * \brief Wrapper around cblas_ztrsm
 * \sa    cblas_ztrsm
 */
void ztrsm(const int layout,
           const int side,
           const int uplo,
           const int transa,
           const int diag,
           const int m,
           const int n,
           const void* alpha,
           const void* a,
           const int lda,
           void* b,
           const int ldb);


/**
 * \brief Wrapper around cblas_ztrsv
 * \sa    cblas_ztrsv
 */
void ztrsv(const int layout,
           const int uplo,
           const int transa,
           const int diag,
           const int n,
           const void* a,
           const int lda,
           void* x,
           const int incx);


} // namespace Cblas
} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_COMMON_BLAS_HH
