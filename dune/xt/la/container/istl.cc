// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include "config.h"

#if HAVE_DUNE_ISTL && HAVE_DUNE_PYBINDXI

#include "istl.hh"

namespace Dune {
namespace XT {
namespace LA {


template class IstlDenseVector<double>;

template class IstlRowMajorSparseMatrix<double>;


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // HAVE_DUNE_ISTL && HAVE_DUNE_PYBINDXI
