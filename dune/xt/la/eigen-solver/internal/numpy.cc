// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#include "config.h"

#if HAVE_DUNE_PYBINDXI
#include <dune/pybindxi/interpreter.hh>
#endif

namespace Dune {
namespace XT {
namespace LA {
namespace internal {


bool numpy_eigensolver_available()
{
#if HAVE_DUNE_PYBINDXI
  try {
    PybindXI::GlobalInterpreter().import_module("numpy.linalg");
  } catch (...) {
    return false;
  }
  return true;
#else // HAVE_DUNE_PYBINDXI
  return false;
#endif
} // ... numpy_eigensolver_available(...)


} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune
