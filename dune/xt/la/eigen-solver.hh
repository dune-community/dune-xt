// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_LA_EIGEN_SOLVER_HH
#define DUNE_XT_LA_EIGEN_SOLVER_HH

#include <complex>
#include <limits>

#include <dune/common/dynvector.hh>
#include <dune/common/typetraits.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/container.hh>

namespace Dune {
namespace XT {
namespace LA {


template <class MatrixImp>
class EigenSolver
{
  static_assert(
      AlwaysFalse<MatrixImp>::value,
      "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for this matrix!");

public:
  typedef MatrixImp MatrixType;
  typedef int RealType;
  typedef int ComplexVectorType;
  typedef int RealVectorType;

  EigenSolver(const MatrixType& /*matrix*/)
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
  }

  static std::vector<std::string> types()
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >. "
               "Please include the correct header for your matrix implementation '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
  }

  static Common::Configuration options(const std::string /*type*/ = "")
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::Solver< ... >. "
               "Please include the correct header for your matrix implementation '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
  }

  std::vector<std::complex<RealType>> eigenvalues() const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return std::vector<std::complex<RealType>>();
  }

  std::vector<std::complex<RealType>> eigenvalues(const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return std::vector<std::complex<RealType>>();
  }

  std::vector<std::complex<RealType>> eigenvalues(const Common::Configuration& /*opts*/) const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return std::vector<std::complex<RealType>>();
  }

  std::vector<RealType> min_eigenvalues(const size_t /*num_evs*/ = std::numeric_limits<size_t>::max()) const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return std::vector<RealType>();
  }

  std::vector<RealType> min_eigenvalues(const std::string& /*type*/,
                                        const size_t /*num_evs*/ = std::numeric_limits<size_t>::max()) const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return std::vector<RealType>();
  }

  std::vector<RealType> min_eigenvalues(const Common::Configuration& /*opts*/,
                                        const size_t /*num_evs*/ = std::numeric_limits<size_t>::max()) const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return std::vector<RealType>();
  }

  std::vector<RealType> max_eigenvalues(const size_t /*num_evs*/ = std::numeric_limits<size_t>::max()) const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return std::vector<RealType>();
  }

  std::vector<RealType> max_eigenvalues(const std::string& /*type*/,
                                        const size_t /*num_evs*/ = std::numeric_limits<size_t>::max()) const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return std::vector<RealType>();
  }

  std::vector<RealType> max_eigenvalues(const Common::Configuration& /*opts*/,
                                        const size_t /*num_evs*/ = std::numeric_limits<size_t>::max()) const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return std::vector<RealType>();
  }

  std::vector<ComplexVectorType> eigenvectors() const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return std::vector<ComplexVectorType>();
  }

  std::vector<ComplexVectorType> eigenvectors(const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return std::vector<ComplexVectorType>();
  }

  std::vector<ComplexVectorType> eigenvectors(const Common::Configuration& /*opts*/) const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return std::vector<ComplexVectorType>();
  }

  std::vector<RealVectorType> real_eigenvectors() const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return std::vector<RealVectorType>();
  }

  std::vector<RealVectorType> real_eigenvectors(const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return std::vector<RealVectorType>();
  }

  std::vector<RealVectorType> real_eigenvectors(const Common::Configuration& /*opts*/) const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return std::vector<RealVectorType>();
  }

}; // class EigenSolver


template <class M>
typename std::enable_if<XT::Common::MatrixAbstraction<M>::is_matrix, EigenSolver<M>>::type
make_eigen_solver(const M& matrix)
{
  return EigenSolver<M>(matrix);
}


} // namespace LA
} // namespace XT
} // namespace Dune

#include "eigen-solver/eigen.hh"
#include "eigen-solver/fmatrix.hh"

#endif // DUNE_XT_LA_EIGEN_SOLVER_HH
