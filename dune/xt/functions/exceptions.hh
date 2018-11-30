// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   TiKeil          (2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_XT_FUNCTIONS_EXCEPTIONS_HH
#define DUNE_XT_FUNCTIONS_EXCEPTIONS_HH

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/grid/exceptions.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace Exceptions {


class wrong_input_given : public Common::Exceptions::wrong_input_given
{
};

class not_bound_to_an_element_yet : public Grid::Exceptions::not_bound_to_an_element_yet
{
};

class reinterpretation_error : public Dune::Exception
{
};

class parameter_error : public Common::Exceptions::parameter_error
{
};

class spe10_data_file_missing : public Dune::IOError
{
};


} // namespace Exceptions
} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_EXCEPTIONS_HH
