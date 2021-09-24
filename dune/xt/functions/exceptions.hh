// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018, 2020)
//   René Fritze     (2018 - 2020)
//   Tim Keil        (2018)
//   Tobias Leibner  (2018, 2020)

#ifndef DUNE_XT_FUNCTIONS_EXCEPTIONS_HH
#define DUNE_XT_FUNCTIONS_EXCEPTIONS_HH

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/grid/exceptions.hh>

namespace Dune::XT::Functions::Exceptions {


class wrong_input_given : public Common::Exceptions::wrong_input_given
{};

class not_bound_to_an_element_yet : public Grid::Exceptions::not_bound_to_an_element_yet
{};

class reinterpretation_error : public Dune::Exception
{};

class parameter_error : public Common::Exceptions::parameter_error
{};

class spe10_data_file_missing : public Dune::IOError
{};

class element_function_error : public Dune::Exception
{};

class combined_error : public Dune::Exception
{};

class function_error : public Dune::Exception
{};

class grid_function_error : public Dune::Exception
{};

class flux_function_error : public Dune::Exception
{};


} // namespace Dune::XT::Functions::Exceptions

#endif // DUNE_XT_FUNCTIONS_EXCEPTIONS_HH
