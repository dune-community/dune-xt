// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_LA_EXCEPTIONS_HH
#define DUNE_XT_LA_EXCEPTIONS_HH

#include <dune/xt/common/exceptions.hh>

namespace Dune {
namespace XT {
namespace LA {
namespace Exceptions {


class linear_solver_failed : public Dune::Exception
{
};

class linear_solver_failed_bc_data_did_not_fulfill_requirements : public linear_solver_failed
{
};

class linear_solver_failed_bc_it_did_not_converge : public linear_solver_failed
{
};

class linear_solver_failed_bc_it_was_not_set_up_correctly : public linear_solver_failed
{
};

class linear_solver_failed_bc_the_solution_does_not_solve_the_system : public linear_solver_failed
{
};


class eigen_solver_failed : public Dune::Exception
{
};

class eigen_solver_failed_bc_data_did_not_fulfill_requirements : public eigen_solver_failed
{
};

class eigen_solver_failed_bc_eigenvalues_are_not_real_as_requested : public eigen_solver_failed
{
};

class eigen_solver_failed_bc_eigenvalues_are_not_positive_as_requested : public eigen_solver_failed
{
};

class eigen_solver_failed_bc_eigenvectors_are_not_real_as_requested : public eigen_solver_failed
{
};


class matrix_invert_failed : public Dune::Exception
{
};

class matrix_invert_failed_bc_data_did_not_fulfill_requirements : public matrix_invert_failed
{
};

class matrix_invert_failed_bc_it_was_not_set_up_correctly : public matrix_invert_failed
{
};

class matrix_invert_failed_bc_result_contained_inf_or_nan : public matrix_invert_failed
{
};

class matrix_invert_failed_bc_result_is_not_a_left_inverse : public matrix_invert_failed
{
};

class matrix_invert_failed_bc_result_is_not_a_right_inverse : public matrix_invert_failed
{
};


} // namespace Exceptions
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EXCEPTIONS_HH
