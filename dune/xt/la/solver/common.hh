// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2014 - 2017)
//   Rene Milk        (2015 - 2018)
//   Tobias Leibner   (2014, 2018)

#ifndef DUNE_XT_LA_SOLVER_COMMON_HH
#define DUNE_XT_LA_SOLVER_COMMON_HH

#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cmath>

#include <dune/xt/common/configuration.hh>

#include <dune/xt/la/algorithms/qr.hh>
#include <dune/xt/la/container/common/matrix/dense.hh>
#include <dune/xt/la/container/common/vector/dense.hh>

#include "../solver.hh"

namespace Dune {
namespace XT {
namespace LA {


template <class S, class CommunicatorType>
class SolverOptions<CommonDenseMatrix<S>, CommunicatorType> : protected internal::SolverUtils
{
public:
  using MatrixType = CommonDenseMatrix<S>;

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
}; // class SolverOptions<CommonDenseMatrix<...>>


template <class S, class CommunicatorType>
class Solver<CommonDenseMatrix<S>, CommunicatorType> : protected internal::SolverUtils
{
public:
  typedef CommonDenseMatrix<S> MatrixType;
  typedef typename MatrixType::RealType R;

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
  } // ... options(...)

  void apply(const CommonDenseVector<S>& rhs, CommonDenseVector<S>& solution) const
  {
    apply(rhs, solution, types()[0]);
  }

  void apply(const CommonDenseVector<S>& rhs, CommonDenseVector<S>& solution, const std::string& type) const
  {
    apply(rhs, solution, options(type));
  }

  void apply(const CommonDenseVector<S>& rhs, CommonDenseVector<S>& solution, const Common::Configuration& opts) const
  {
    if (!opts.has_key("type"))
      DUNE_THROW(Common::Exceptions::configuration_error,
                 "Given options (see below) need to have at least the key 'type' set!\n\n"
                     << opts);
    const auto type = opts.get<std::string>("type");
    internal::SolverUtils::check_given(type, types());
    const Common::Configuration default_opts = options(type);
    // solve
    try {
      auto QR = matrix_;
      solve_by_qr_decomposition(QR, solution, rhs);
    } catch (FMatrixError&) {
      DUNE_THROW(Exceptions::linear_solver_failed_bc_data_did_not_fulfill_requirements,
                 "The dune-common backend reported 'FMatrixError'!\n"
                     << "Those were the given options:\n\n"
                     << opts);
    }
    // check
    const R post_check_solves_system_threshold =
        opts.get("post_check_solves_system", default_opts.get<R>("post_check_solves_system"));
    if (post_check_solves_system_threshold > 0) {
      auto tmp = rhs.copy();
      matrix_.mv(solution, tmp);
      tmp -= rhs;
      const R sup_norm = tmp.sup_norm();
      if (sup_norm > post_check_solves_system_threshold || Common::isnan(sup_norm) || Common::isinf(sup_norm))
        DUNE_THROW(Exceptions::linear_solver_failed_bc_the_solution_does_not_solve_the_system,
                   "The computed solution does not solve the system (although the dune-common backend "
                       << "reported no error) and you requested checking (see options below)! "
                       << "If you want to disable this check, set 'post_check_solves_system = 0' in the options."
                       << "\n\n"
                       << "  (A * x - b).sup_norm() = " << tmp.sup_norm() << "\n\n"
                       << "Those were the given options:\n\n"
                       << opts);
    }
  } // ... apply(...)

private:
  const MatrixType& matrix_;
}; // class Solver< CommonDenseMatrix< ... > >


} // namespace LA
} // namespace XT
} // namespace Dune


// begin: this is what we need for the lib
#if DUNE_XT_WITH_PYTHON_BINDINGS


extern template class Dune::XT::LA::Solver<Dune::XT::LA::CommonDenseMatrix<double>>;


#endif // DUNE_XT_WITH_PYTHON_BINDINGS
// end: this is what we need for the lib


#endif // DUNE_XT_LA_SOLVER_COMMON_HH
