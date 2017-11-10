// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2017)

#include "lapacke.hh"

#if HAVE_LAPACKE

#include <dune/xt/common/lapacke.hh>

namespace Dune {
namespace XT {
namespace LA {
namespace internal {


int LapackeWrapper::dgeev(char jobvl,
                          char jobvr,
                          int n,
                          double* a,
                          int lda,
                          double* wr,
                          double* wi,
                          double* vl,
                          int ldvl,
                          double* vr,
                          int ldvr)
{
  return LAPACKE_dgeev(LAPACK_ROW_MAJOR, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr);
}


} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // HAVE_LAPACKE
