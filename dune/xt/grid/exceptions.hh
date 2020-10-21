// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014, 2016 - 2018)
//   René Fritze     (2012 - 2013, 2015 - 2016, 2018 - 2019)
//   Tobias Leibner  (2014, 2016, 2020)

#ifndef DUNE_XT_GRID_EXCEPTIONS_HH
#define DUNE_XT_GRID_EXCEPTIONS_HH

#include <dune/xt/common/exceptions.hh>

namespace Dune::XT::Grid::Exceptions {


class boundary_type_error : public Dune::Exception
{};

class boundary_info_error : public Dune::Exception
{};

class not_bound_to_an_element_yet : public Dune::InvalidStateException
{};

class functor_error : public Dune::InvalidStateException
{};

class wrong_dimension : public Dune::Exception
{};

class wrong_codimension : public Dune::Exception
{};


} // namespace Dune::XT::Grid::Exceptions

#endif // DUNE_XT_GRID_EXCEPTIONS_HH
