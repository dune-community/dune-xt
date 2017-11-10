// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2017)

#include "config.h"

#include "lapacke.hh"

#if HAVE_LAPACKE

#include <dune/xt/common/lapacke.hh>

namespace Dune {
namespace XT {
namespace LA {
namespace internal {


UnitMatrix::UnitMatrix(int N)
  : unit_matrix_(new std::vector<double>(N, 0.))
{
}

// need to reset unit_matrix every time because Lapack changes it
double* UnitMatrix::get(int N)
{
  unit_matrix_->resize(N * N);
  std::fill(unit_matrix_->begin(), unit_matrix_->end(), 0.);
  for (int rr = 0; rr < N; ++rr)
    (*unit_matrix_)[rr * (N + 1)] = 1.;
  return unit_matrix_->data();
}

int LapackeWrapper::dggev(char jobvl,
                          char jobvr,
                          int n,
                          double* a,
                          int lda,
                          double* b,
                          int ldb,
                          double* alphar,
                          double* alphai,
                          double* beta,
                          double* vl,
                          int ldvl,
                          double* vr,
                          int ldvr)
{
  return LAPACKE_dggev(LAPACK_ROW_MAJOR, jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr);
}


} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // HAVE_LAPACKE
