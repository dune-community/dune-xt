// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_XT_GRID_EXCEPTIONS_HH
#define DUNE_XT_GRID_EXCEPTIONS_HH

#include <dune/xt/common/exceptions.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace Exceptions {


class boundary_type_error : public Dune::Exception
{
};


class boundary_info_error : public Dune::Exception
{
};


} // namespace Exceptions
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_EXCEPTIONS_HH
