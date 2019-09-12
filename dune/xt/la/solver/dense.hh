// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_XT_LA_SOLVER_DENSE_HH
#define DUNE_XT_LA_SOLVER_DENSE_HH

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/vector.hh>

#include <dune/xt/la/algorithms/qr.hh>
#include <dune/xt/la/type_traits.hh>

#include "../solver.hh"

namespace Dune {
namespace XT {
namespace LA {


template <class Matrix>
class SolverOptions<
    Matrix,
    std::enable_if_t<XT::Common::is_matrix<Matrix>::value && !is_matrix<Matrix>::value, SequentialCommunication>>
  : protected internal::SolverUtils
{
public:
  using MatrixType = Matrix;
  using CommunicatorType = void;

  static std::vector<std::string> types()
  {
    return {"qr.householder"};
  }

  static Common::Configuration options(const std::string type = "")
  {
    const std::string tp = !type.empty() ? type : types()[0];
    internal::SolverUtils::check_given(tp, types());
    return Common::Configuration({"type", "post_check_solves_system"}, {tp.c_str(), "1e-5"});
  }
}; // class SolverOptions<...>


template <class Matrix>
class Solver<
    Matrix,
    std::enable_if_t<XT::Common::is_matrix<Matrix>::value && !is_matrix<Matrix>::value, SequentialCommunication>>
  : protected internal::SolverUtils
{
  using M = XT::Common::MatrixAbstraction<Matrix>;

public:
  using MatrixType = Matrix;
  using CommunicatorType = SequentialCommunication;

  Solver(const MatrixType& matrix)
    : matrix_(matrix)
  {}

  Solver(const MatrixType& matrix, const CommunicatorType& /*communicator*/)
    : matrix_(matrix)
  {}

  static std::vector<std::string> types()
  {
    return {"qr.householder"};
  }

  static Common::Configuration options(const std::string type = "")
  {
    const std::string tp = !type.empty() ? type : types()[0];
    internal::SolverUtils::check_given(tp, types());
    return Common::Configuration({"type", "post_check_solves_system"}, {tp.c_str(), "1e-5"});
  }

  template <class VectorType>
  void apply(const VectorType& rhs, VectorType& solution) const
  {
    apply(rhs, solution, types()[0]);
  }

  template <class VectorType>
  void apply(const VectorType& rhs, VectorType& solution, const std::string& type) const
  {
    apply(rhs, solution, options(type));
  }

  template <class VectorType>
  std::enable_if_t<XT::Common::is_vector<VectorType>::value, void>
  apply(const VectorType& rhs, VectorType& solution, const Common::Configuration& opts) const
  {
    if (!opts.has_key("type"))
      DUNE_THROW(Common::Exceptions::configuration_error,
                 "Given options (see below) need to have at least the key 'type' set!\n\n"
                     << opts);
    const auto type = opts.get<std::string>("type");
    internal::SolverUtils::check_given(type, types());
    const Common::Configuration default_opts = options(type);
    // solve
    auto writable_copy_of_matrix_ = matrix_;
    solve_by_qr_decomposition(writable_copy_of_matrix_, solution, rhs);
    // check
    const auto post_check_solves_system_threshold =
        opts.get("post_check_solves_system", default_opts.get<double>("post_check_solves_system"));
    if (post_check_solves_system_threshold > 0) {
      auto tmp = XT::Common::zeros_like(rhs);
      XT::Common::mv(matrix_, solution, tmp);
      tmp -= rhs;
      const auto sup_norm = XT::Common::sup_norm(tmp);
      if (sup_norm > post_check_solves_system_threshold || Common::isnan(sup_norm) || Common::isinf(sup_norm))
        DUNE_THROW(Exceptions::linear_solver_failed_bc_the_solution_does_not_solve_the_system,
                   "The computed solution does not solve the system and you requested checking (see options below)! "
                       << "If you want to disable this check, set 'post_check_solves_system = 0' in the options."
                       << "\n\n"
                       << "  (A * x - b).sup_norm() = " << sup_norm << "\n\n"
                       << "Those were the given options:\n\n"
                       << opts);
    }
  } // ... apply(...)

private:
  const MatrixType& matrix_;
}; // class Solver<...>


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_SOLVER_DENSE_HH
