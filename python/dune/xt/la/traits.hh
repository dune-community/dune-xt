// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_XT_LA_TRAITS_HH
#define PYTHON_DUNE_XT_LA_TRAITS_HH

#include <dune/pybindxi/pybind11.h>

namespace Dune {
namespace XT {
namespace LA {
namespace bindings {


struct Common
{};

struct Eigen
{};

struct Istl
{};


struct Dense
{};

struct Sparse
{};


} // namespace bindings
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // PYTHON_DUNE_XT_LA_TRAITS_HH
