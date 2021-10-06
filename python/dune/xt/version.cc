// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <dune/pybindxi/pybind11.h>

// see https://stackoverflow.com/questions/240353/convert-a-preprocessor-token-to-a-string
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)


PYBIND11_MODULE(_version, m)
{
  m.attr("__version__") = pybind11::str(TOSTRING(DUNE_XT_VERSION));
}


#undef TOSTRING
#undef STRINGIFY