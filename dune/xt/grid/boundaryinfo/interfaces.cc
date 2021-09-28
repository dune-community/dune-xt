// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014, 2016 - 2017)
//   René Fritze     (2012 - 2013, 2015 - 2016, 2018 - 2020)
//   Tobias Leibner  (2014, 2016, 2020)

#include <config.h>

#include <ostream>

#include "interfaces.hh"

namespace Dune::XT::Grid {


std::ostream& operator<<(std::ostream& out, const BoundaryType& type)
{
  out << type.id();
  return out;
}


} // namespace Dune::XT::Grid
