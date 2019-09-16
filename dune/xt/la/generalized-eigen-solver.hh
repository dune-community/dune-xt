// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017, 2019)
//   René Fritze     (2018)
//   Tobias Leibner  (2017 - 2018)

#ifndef DUNE_XT_LA_GENERALIZED_EIGEN_SOLVER_HH
#define DUNE_XT_LA_GENERALIZED_EIGEN_SOLVER_HH

#include <string>
#include <vector>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/type_traits.hh>

namespace Dune {
namespace XT {
namespace LA {


/**
 * \brief A means to obtain available options at compile time.
 * \note  This class needs to be specialized for each MatrixType, the purpose of this variant is merely to document the
 *        expected functionality.
 */
template <class MatrixType, bool is_matrix = XT::Common::is_matrix<MatrixType>::value>
class GeneralizedEigenSolverOptions
{
  static_assert(AlwaysFalse<MatrixType>::value,
                "Please implement for given MatrixType and add the respective include below!");

  static std::vector<std::string> types();

  static Common::Configuration options(const std::string /*type*/ = "");
}; // class GeneralizedEigenSolverOptions


template <class MatrixType>
std::vector<std::string> generalized_eigen_solver_types(const MatrixType& /*matrix*/)
{
  return GeneralizedEigenSolverOptions<MatrixType>::types();
}


template <class MatrixType>
Common::Configuration generalized_eigen_solver_options(const MatrixType& /*matrix*/, const std::string type = "")
{
  return GeneralizedEigenSolverOptions<MatrixType>::options(type);
}


template <class MatrixImp, bool is_matrix = XT::Common::is_matrix<MatrixImp>::value>
class GeneralizedEigenSolver
{
  static_assert(AlwaysFalse<MatrixImp>::value,
                "Please implement for given MatrixType and add the respective include below!");

public:
  typedef MatrixImp MatrixType;
  typedef double FieldType;
  typedef int RealMatrixType;
  typedef int ComplexMatrixType;

  GeneralizedEigenSolver(const MatrixType& /*matrix*/, const std::string& /*type*/ = "")
  {
    static_assert(AlwaysFalse<MatrixType>::value,
                  "Please implement for given MatrixType and add the respective include below!");
  }

  GeneralizedEigenSolver(const MatrixType& /*matrix*/, const Common::Configuration /*opts*/)
  {
    static_assert(AlwaysFalse<MatrixType>::value,
                  "Please implement for given MatrixType and add the respective include below!");
  }

  const Common::Configuration& options() const;

  const MatrixType& lhs_matrix() const;

  const MatrixType& rhs_matrix() const;

  const std::vector<Common::complex_t<FieldType>>& eigenvalues() const;

  const std::vector<Common::real_t<FieldType>>& real_eigenvalues() const;

  const std::vector<Common::real_t<FieldType>>&
  min_eigenvalues(const size_t /*num_evs*/ = std::numeric_limits<size_t>::max()) const;

  const std::vector<Common::real_t<FieldType>>&
  max_eigenvalues(const size_t /*num_evs*/ = std::numeric_limits<size_t>::max()) const;

  const ComplexMatrixType& eigenvectors() const;

  const RealMatrixType& real_eigenvectors() const;
}; // class GeneralizedEigenSolver


template <class M>
GeneralizedEigenSolver<M>
make_generalized_eigen_solver(const M& lhs_matrix, const M& rhs_matrix, const std::string& type = "")
{
  return GeneralizedEigenSolver<M>(lhs_matrix, rhs_matrix, type);
}


template <class M>
GeneralizedEigenSolver<M>
make_generalized_eigen_solver(const M& lhs_matrix, const M& rhs_matrix, const XT::Common::Configuration& options)
{
  return GeneralizedEigenSolver<M>(lhs_matrix, rhs_matrix, options);
}


} // namespace LA
} // namespace XT
} // namespace Dune

#include "generalized-eigen-solver/default.hh"

#endif // DUNE_XT_LA_GENERALIZED_EIGEN_SOLVER_HH
