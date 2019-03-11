// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner   (2019)

#ifndef DUNE_XT_LA_SOLVER_VIEW_HH
#define DUNE_XT_LA_SOLVER_VIEW_HH

#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cmath>

#include <dune/xt/common/configuration.hh>

#include <dune/xt/la/algorithms/qr.hh>
#include <dune/xt/la/container/matrix-view.hh>
#include <dune/xt/la/container/vector-view.hh>

#include "../solver.hh"

namespace Dune {
namespace XT {
namespace LA {


template <class MatrixImp, class CommunicatorType>
class SolverOptions<MatrixView<MatrixImp>, CommunicatorType> : protected internal::SolverUtils
{
public:
  using MatrixType = MatrixView<MatrixImp>;

  static std::vector<std::string> types()
  {
    return SolverOptions<MatrixImp, CommunicatorType>::types();
  }

  static Common::Configuration options(const std::string type = "")
  {
    return SolverOptions<MatrixImp, CommunicatorType>::options(type);
  }
}; // class SolverOptions<MatrixView<...>>


template <class MatrixImp, class CommunicatorType>
class Solver<MatrixView<MatrixImp>, CommunicatorType> : protected internal::SolverUtils
{
public:
  typedef MatrixView<MatrixImp> MatrixType;
  typedef typename MatrixType::RealType R;
  using ActualSolver = Solver<MatrixImp, CommunicatorType>;

  Solver(const MatrixType& matrix_view)
    : matrix_view_(matrix_view)
    , matrix_(matrix_view_.rows(), matrix_view_.cols(), matrix_view_.pattern())
    , actual_solver_(matrix_)
  {
    const auto pattern = matrix_view_.pattern();
    for (size_t ii = 0; ii < matrix_view_.rows(); ++ii)
      for (const auto& jj : pattern.inner(ii))
        matrix_.set_entry(ii, jj, matrix_view_.get_entry(ii, jj));
  }

  Solver(const MatrixType& matrix, const CommunicatorType& /*communicator*/)
    : Solver(matrix)
  {}

  static std::vector<std::string> types()
  {
    SolverOptions<MatrixType, CommunicatorType>::types();
  }

  static Common::Configuration options(const std::string type = "")
  {
    SolverOptions<MatrixType, CommunicatorType>::options(type);
  } // ... options(...)

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
  void apply(const VectorType& rhs, VectorType& solution, const Common::Configuration& opts) const
  {
    std::cerr << "Warning: Currently, Solver for MatrixView just copies everything to a new matrix, may be slow!"
              << std::endl;
    actual_solver_.apply(rhs, solution, opts);
  } // ... apply(...)

  template <class VectorType>
  void apply(const VectorView<VectorType>& rhs, VectorView<VectorType>& solution) const
  {
    apply(rhs, solution, types()[0]);
  }

  template <class VectorType>
  void apply(const VectorView<VectorType>& rhs, VectorView<VectorType>& solution, const std::string& type) const
  {
    apply(rhs, solution, options(type));
  }

  template <class VectorType>
  void
  apply(const VectorView<VectorType>& rhs, VectorView<VectorType>& solution, const Common::Configuration& opts) const
  {
    std::cerr << "Warning: Currently, Solver for MatrixView just copies everything to a new matrix, may be slow!"
              << std::endl;
    VectorType actual_rhs(rhs.size()), actual_solution(solution.size());
    for (size_t ii = 0; ii < rhs.size(); ++ii)
      actual_rhs[ii] = rhs[ii];
    apply(actual_rhs, actual_solution, opts);
    for (size_t ii = 0; ii < solution.size(); ++ii)
      actual_solution[ii] = solution[ii];
  } // ... apply(...)

private:
  const MatrixType& matrix_view_;
  MatrixImp matrix_;
  ActualSolver actual_solver_;
}; // class Solver< MatrixView< ... > >


} // namespace LA
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_LA_SOLVER_VIEW_HH
