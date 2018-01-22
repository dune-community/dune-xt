// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_LA_EIGEN_SOLVER_HH
#define DUNE_XT_LA_EIGEN_SOLVER_HH

#include <complex>
#include <limits>

#include <dune/common/typetraits.hh>

#include <dune/common/dynvector.hh>
#include <dune/common/typetraits.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/container.hh>

namespace Dune {
namespace XT {
namespace LA {


/**
 * \brief A means to obtain available options at compile time.
 * \note  This class needs to be specialized for each MatrixType, the purpose of this variant is merely to document the
 *        expected functionality.
 */
template <class MatrixType>
class EigenSolverOptions
{
  static_assert(AlwaysFalse<MatrixType>::value,
                "Please implement for given MatrixType and add the respective include below!");

  static std::vector<std::string> types();

  static Common::Configuration options(const std::string /*type*/ = "");
}; // class EigenSolverOptions


template <class MatrixType>
std::vector<std::string> eigen_solver_types(const MatrixType& /*matrix*/)
{
  return EigenSolverOptions<MatrixType>::types();
}


template <class MatrixType>
Common::Configuration eigen_solver_options(const MatrixType& /*matrix*/, const std::string type = "")
{
  return EigenSolverOptions<MatrixType>::options(type);
}


template <class MatrixImp>
class EigenSolver
{
  static_assert(AlwaysFalse<MatrixImp>::value,
                "Please implement for given MatrixType and add the respective include below!");

public:
  typedef MatrixImp MatrixType;
  typedef double FieldType;
  typedef int RealMatrixType;
  typedef int ComplexMatrixType;

  EigenSolver(const MatrixType& /*matrix*/, const std::string& /*type*/ = "")
  {
    static_assert(AlwaysFalse<MatrixType>::value,
                  "Please implement for given MatrixType and add the respective include below!");
  }

  EigenSolver(const MatrixType& /*matrix*/, const Common::Configuration /*opts*/)
  {
    static_assert(AlwaysFalse<MatrixType>::value,
                  "Please implement for given MatrixType and add the respective include below!");
  }

  const Common::Configuration& options() const;

  const MatrixType& matrix() const;

  const std::vector<Common::complex_t<FieldType>>& eigenvalues() const;

  const std::vector<Common::real_t<FieldType>>& real_eigenvalues() const;

  const std::vector<Common::real_t<FieldType>>&
  min_eigenvalues(const size_t /*num_evs*/ = std::numeric_limits<size_t>::max()) const;

  const std::vector<Common::real_t<FieldType>>&
  max_eigenvalues(const size_t /*num_evs*/ = std::numeric_limits<size_t>::max()) const;

  const ComplexMatrixType& eigenvectors() const;

  const ComplexMatrixType& eigenvectors_inverse() const;

  const RealMatrixType& real_eigenvectors() const;

  const RealMatrixType& real_eigenvectors_inverse() const;
}; // class EigenSolver


template <class M>
EigenSolver<M> make_eigen_solver(const M& matrix, const std::string& type = "")
{
  return EigenSolver<M>(matrix, type);
}


template <class M>
EigenSolver<M> make_eigen_solver(const M& matrix, const XT::Common::Configuration& options)
{
  return EigenSolver<M>(matrix, options);
}


} // namespace LA
} // namespace XT
} // namespace Dune

#include "eigen-solver/eigen.hh"
#include "eigen-solver/fmatrix.hh"

#endif // DUNE_XT_LA_EIGEN_SOLVER_HH
