// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2020)
//
// reserved. License: Dual licensed as BSD 2-Clause License
// (http://opensource.org/licenses/BSD-2-Clause)

#warning This header is deprecated, include <dune/xt/grid/grids.hh> instead (31.07.2019)!

#ifndef PYTHON_DUNE_XT_GRID_TYPES_HH
#  define PYTHON_DUNE_XT_GRID_TYPES_HH

#  include <dune/xt/common/deprecated.hh>
#  include <dune/xt/grid/grids.hh>

namespace Dune::XT::Grid::bindings {


using AvailableTypes
    DXT_DEPRECATED_MSG("Use XT::Grid::AvailableGridTypes instead (31.07.2019)!") = Dune::XT::Grid::AvailableGridTypes;


} // namespace Dune::XT::Grid::bindings

#endif // PYTHON_DUNE_XT_GRID_TYPES_HH
