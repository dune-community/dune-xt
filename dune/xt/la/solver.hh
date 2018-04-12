// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2014, 2016 - 2017)
//   Rene Milk       (2013, 2015 - 2016, 2018)
//   Tobias Leibner  (2014, 2018)

#ifndef DUNE_XT_LA_SOLVER_HH
#define DUNE_XT_LA_SOLVER_HH

#include <string>
#include <vector>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/parallel/helper.hh>

#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/type_traits.hh>

namespace Dune {
namespace XT {
namespace LA {
namespace internal {


static const constexpr size_t max_size_to_print = 5;


class SolverUtils
{
public:
  static void check_given(const std::string& type, const std::vector<std::string>& opts)
  {
    if (std::find(opts.begin(), opts.end(), type) == opts.end()) {
      std::stringstream ss;
      for (auto opt : opts)
        ss << opt << " ";
      DUNE_THROW(Common::Exceptions::configuration_error,
                 "Given type '" << type << "' not supported (see below for a list of supported ones). "
                                << "Call options() first!\n"
                                << ss.str());
    }
  }
};


} // namespace internal


template <class MatrixImp, class CommunicatorType = SequentialCommunication>
class Solver
{
  static_assert(AlwaysFalse<MatrixImp>::value,
                "This is the unspecialized version of LA::Solver< ... >. Please include "
                "the correct header for your matrix implementation!");

public:
  typedef MatrixImp MatrixType;

  Solver(const MatrixType& /*matrix*/, const CommunicatorType&)
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::Solver< ... >. "
               "Please include the correct header for your matrix implementation '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
  }

  Solver(const MatrixType& /*matrix*/)
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::Solver< ... >. "
               "Please include the correct header for your matrix implementation '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
  }

  static std::vector<std::string> types()
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::Solver< ... >. "
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

  /**
   *  Throws any of the above exceptions, if there was a problem. If none was thrown we beleive that a suitable solution
   *  was found (given the current options).
   */
  template <class RhsType, class SolutionType>
  void apply(const RhsType& /*rhs*/, SolutionType& /*solution*/) const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::Solver< ... >. "
               "Please include the correct header for your matrix implementation '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
  }

  template <class RhsType, class SolutionType>
  void apply(const RhsType& /*rhs*/, SolutionType& /*solution*/, const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::Solver< ... >. "
               "Please include the correct header for your matrix implementation '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
  }

  template <class RhsType, class SolutionType>
  void apply(const RhsType& /*rhs*/, SolutionType& /*solution*/, const Common::Configuration& /*options*/) const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::Solver< ... >. "
               "Please include the correct header for your matrix implementation '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
  }
}; // class Solver


template <class M>
typename std::enable_if<XT::LA::is_matrix<M>::value, Solver<M>>::type make_solver(const M& matrix)
{
  return Solver<M>(matrix);
}


template <class M, class V, class... Args>
typename std::enable_if<XT::LA::is_matrix<M>::value && XT::LA::is_vector<V>::value, void>::type
solve(const M& A, const V& b, V& x, Args&&... args)
{
  make_solver(A).apply(b, x, std::forward<Args>(args)...);
}


template <class M, class C>
typename std::enable_if<XT::LA::is_matrix<M>::value, Solver<M, C>>::type make_solver(const M& matrix, const C& dof_comm)
{
  return Solver<M, C>(matrix, dof_comm);
}


template <class M, class V, class C, class... Args>
typename std::enable_if<XT::LA::is_matrix<M>::value && XT::LA::is_vector<V>::value, void>::type
solve(const M& A, const V& b, V& x, const C& dof_comm, Args&&... args)
{
  make_solver(A, dof_comm).apply(b, x, std::forward<Args>(args)...);
}


template <class M, class V, class... Args>
typename std::enable_if<XT::LA::is_matrix<M>::value && XT::LA::is_vector<V>::value, V>::type
solve(const M& A, const V& b, Args&&... args)
{
  V x(A.cols());
  solve(A, b, x, std::forward<Args>(args)...);
  return x;
}


} // namespace LA
} // namespace XT
} // namespace Dune

#include "solver/common.hh"
#include "solver/eigen.hh"
#include "solver/istl.hh"

#endif // DUNE_XT_LA_SOLVER_HH
