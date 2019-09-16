// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_COMMON_EXCEPTIONS_BINDINGS_HH
#define DUNE_XT_COMMON_EXCEPTIONS_BINDINGS_HH

#include <dune/xt/common/deprecated.hh>

namespace pybind11 {
class module;
}

namespace Dune {
namespace XT {
namespace Common {
namespace bindings {


DXT_DEPRECATED_MSG("you no longer need to call this for your module")
void addbind_exceptions(pybind11::module&) {} // ... addbind_exceptions(...)


} // namespace bindings
} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_COMMON_EXCEPTIONS_BINDINGS_HH
