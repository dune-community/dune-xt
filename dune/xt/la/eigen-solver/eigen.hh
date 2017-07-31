// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_LA_EIGEN_SOLVER_EIGEN_HH
#define DUNE_XT_LA_EIGEN_SOLVER_EIGEN_HH

#include <algorithm>
#include <functional>

#include <dune/xt/la/container/eigen/dense.hh>
#include <dune/xt/la/solver.hh>

#include "../eigen-solver.hh"
#include "internal/eigen.hh"

namespace Dune {
namespace XT {
namespace LA {

//#if HAVE_EIGEN


template <class S>
class EigenSolver<EigenDenseMatrix<S>> : protected internal::SolverUtils
{
public:
  typedef EigenDenseMatrix<S> MatrixType;
  typedef typename MatrixType::RealType RealType;
  typedef
      typename XT::LA::Container<std::complex<RealType>, MatrixType::vector_type>::VectorType ComplexEigenvectorType;
  typedef typename XT::LA::Container<RealType, MatrixType::vector_type>::VectorType RealEigenvectorType;

  EigenSolver(const MatrixType& matrix)
    : matrix_(matrix)
  {
  }

  static std::vector<std::string> types()
  {
    return {"default"};
  }

  static Common::Configuration options(const std::string type = "")
  {
    const std::string tp = !type.empty() ? type : types()[0];
    internal::SolverUtils::check_given(tp, types());
    Common::Configuration default_options(
        {"type", "check_for_inf_nan", "check_evs_are_real", "check_evs_are_positive", "check_eigenvectors_are_real"},
        {tp, "1", "0", "0", "0"});
    return default_options;
  }

  std::vector<std::complex<RealType>> eigenvalues() const
  {
    return eigenvalues(types().at(0));
  }

  std::vector<std::complex<RealType>> eigenvalues(const std::string type) const
  {
    return eigenvalues(options(type));
  }

  std::vector<std::complex<RealType>> eigenvalues(const Common::Configuration& opts) const
  {
    if (!opts.has_key("type"))
      DUNE_THROW(Common::Exceptions::configuration_error,
                 "Given options (see below) need to have at least the key 'type' set!\n\n"
                     << opts);
    const auto type = opts.get<std::string>("type");
    internal::SolverUtils::check_given(type, types());
    const Common::Configuration default_opts = options(type);
    // checks
    if (matrix_.rows() != matrix_.cols())
      DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
                 "Matrix has to be square!\n   matrix.rows() = " << matrix_.rows() << "\n   matrix.cols()"
                                                                 << matrix_.cols());
    const bool check_for_inf_nan = opts.get("check_for_inf_nan", default_opts.get<bool>("check_for_inf_nan"));
    if (check_for_inf_nan) {
      for (size_t ii = 0; ii < matrix_.rows(); ++ii) {
        for (size_t jj = 0; jj < matrix_.cols(); ++jj) {
          const S& val = matrix_.backend()(ii, jj);
          if (Common::isnan(val) || Common::isinf(val)) {
            std::stringstream msg;
            msg << "Given matrix contains inf or nan and you requested checking (see options below)!\n"
                << "If you want to disable this check, set 'check_for_inf_nan = 0' in the options.\n\n"
                << "Those were the given options:\n\n"
                << opts;
            if (matrix_.rows() <= internal::max_size_to_print && matrix_.cols() <= internal::max_size_to_print)
              msg << "\nThis was the given matrix:\n\n" << matrix_ << "\n";
            DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements, msg.str());
          }
        }
      }
    }
    // solve
    std::vector<std::complex<RealType>> evs;
    if (type == "default")
      evs = internal::compute_all_eigenvalues_using_eigen(matrix_.backend());
    else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type << "' is not supported, although it was reported by types()!");
    // checks
    if (check_for_inf_nan)
      for (const auto& ev : evs) {
        if (Common::isnan(ev) || Common::isinf(ev)) {
          std::stringstream msg;
          msg << "The computed eigenvalues contain inf or nan and you requested checking (see options "
              << "below)!\n"
              << "If you want to disable this check, set 'check_for_inf_nan = 0' in the options.\n\n"
              << "Those were the given options:\n\n"
              << opts;
          if (matrix_.rows() <= internal::max_size_to_print && matrix_.cols() <= internal::max_size_to_print)
            msg << "\nThis was the given matrix:\n\n" << matrix_ << "\n\nThese were the computed eigenvalues:";
          for (const auto& print_ev : evs)
            msg << "\n  " << print_ev;
          msg << "\n";
          DUNE_THROW(Exceptions::eigen_solver_failed, msg.str());
        }
      }
    const double check_evs_are_real = opts.get("check_evs_are_real", default_opts.get<double>("check_evs_are_real"));
    const double check_evs_are_positive =
        opts.get("check_evs_are_positive", default_opts.get<double>("check_evs_are_positive"));
    if (check_evs_are_real > 0 || check_evs_are_positive)
      for (size_t ii = 0; ii < evs.size(); ++ii)
        if (std::abs(evs[ii].imag()) > check_evs_are_real)
          DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvalues_are_not_real_as_requested,
                     "Not all eigenvalues are purely real!\n   evs[" << ii << "].imag(): " << evs[ii].imag());
    if (check_evs_are_positive)
      for (size_t ii = 0; ii < evs.size(); ++ii)
        if (evs[ii].real() < check_evs_are_real)
          DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvalues_are_not_real_as_requested,
                     "Not all eigenvalues are positive!\n   evs[" << ii << "].real(): " << evs[ii].real());
    return evs;
  } // ... eigenvalues(...)

  std::vector<RealType> min_eigenvalues(const size_t num_evs = std::numeric_limits<size_t>::max()) const
  {
    return min_eigenvalues(types().at(0), num_evs);
  }

  std::vector<RealType> min_eigenvalues(const std::string type,
                                        const size_t num_evs = std::numeric_limits<size_t>::max()) const
  {
    return min_eigenvalues(options(type), num_evs);
  }

  std::vector<RealType> min_eigenvalues(const Common::Configuration& opts,
                                        const size_t num_evs = std::numeric_limits<size_t>::max()) const
  {
    Common::Configuration options_to_ensure_real_evs = opts;
    if (!options_to_ensure_real_evs.has_key("check_evs_are_real"))
      options_to_ensure_real_evs["check_evs_are_real"] = "1e-15";
    std::vector<RealType> evs;
    try {
      for (const auto& ev : eigenvalues(opts))
        evs.push_back(ev.real());
    } catch (Exceptions::eigen_solver_failed_bc_eigenvalues_are_not_real_as_requested& ee) {
      DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvalues_are_not_real_as_requested,
                 "Computing the minimum eigenvalue does not make sense for complex eigenvalues!\nYou can configure the "
                 "check by adapting check_evs_are_real in the options (see below).\n   This was the "
                 "original error:\n\n"
                     << ee.what()
                     << "\n\nThese were the options:\n"
                     << options);
    }
    std::sort(evs.begin(), evs.end());
    std::vector<RealType> sorted_evs;
    for (size_t ii = 0; ii < std::min(evs.size(), num_evs); ++ii)
      sorted_evs.push_back(evs[ii]);
    return sorted_evs;
  } // ... min_eigenvalues(...)

  std::vector<RealType> max_eigenvalues(const size_t num_evs = std::numeric_limits<size_t>::max()) const
  {
    return max_eigenvalues(types().at(0), num_evs);
  }

  std::vector<RealType> max_eigenvalues(const std::string type,
                                        const size_t num_evs = std::numeric_limits<size_t>::max()) const
  {
    return max_eigenvalues(options(type), num_evs);
  }

  std::vector<RealType> max_eigenvalues(const Common::Configuration& opts,
                                        const size_t num_evs = std::numeric_limits<size_t>::max()) const
  {
    Common::Configuration options_to_ensure_real_evs = opts;
    if (!options_to_ensure_real_evs.has_key("check_evs_are_real"))
      options_to_ensure_real_evs["check_evs_are_real"] = "1e-15";
    std::vector<RealType> evs;
    try {
      for (const auto& ev : eigenvalues(opts))
        evs.push_back(ev.real());
    } catch (Exceptions::eigen_solver_failed_bc_eigenvalues_are_not_real_as_requested& ee) {
      DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvalues_are_not_real_as_requested,
                 "Computing the maximum eigenvalue does not make sense for complex eigenvalues!\nYou can configure the "
                 "check by adapting check_evs_are_real in the options (see below).\n   This was the "
                 "original error:\n\n"
                     << ee.what()
                     << "\n\nThese were the options:\n"
                     << options);
    }
    std::sort(evs.begin(), evs.end(), std::greater<RealType>());
    std::vector<RealType> sorted_evs;
    for (size_t ii = 0; ii < std::min(evs.size(), num_evs); ++ii)
      sorted_evs.push_back(evs[ii]);
    return sorted_evs;
  } // ... max_eigenvalues(...)

  std::vector<ComplexEigenvectorType> eigenvectors() const
  {
    return eigenvectors(types().at(0));
  }

  std::vector<ComplexEigenvectorType> eigenvectors(const std::string& type) const
  {
    return eigenvectors(options(type));
  }

  std::vector<ComplexEigenvectorType> eigenvectors(const Common::Configuration& opts) const
  {
    if (!opts.has_key("type"))
      DUNE_THROW(Common::Exceptions::configuration_error,
                 "Given options (see below) need to have at least the key 'type' set!\n\n"
                     << opts);
    const auto type = opts.get<std::string>("type");
    internal::SolverUtils::check_given(type, types());
    const Common::Configuration default_opts = options(type);
    // checks
    if (matrix_.rows() != matrix_.cols())
      DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
                 "Matrix has to be square!\n   matrix.rows() = " << matrix_.rows() << "\n   matrix.cols()"
                                                                 << matrix_.cols());
    const bool check_for_inf_nan = opts.get("check_for_inf_nan", default_opts.get<bool>("check_for_inf_nan"));
    if (check_for_inf_nan) {
      for (size_t ii = 0; ii < matrix_.rows(); ++ii) {
        for (size_t jj = 0; jj < matrix_.cols(); ++jj) {
          const S& val = matrix_.backend()(ii, jj);
          if (Common::isnan(val) || Common::isinf(val)) {
            std::stringstream msg;
            msg << "Given matrix contains inf or nan and you requested checking (see options below)!\n"
                << "If you want to disable this check, set 'check_for_inf_nan = 0' in the options.\n\n"
                << "Those were the given options:\n\n"
                << opts;
            if (matrix_.rows() <= internal::max_size_to_print && matrix_.cols() <= internal::max_size_to_print)
              msg << "\nThis was the given matrix:\n\n" << matrix_ << "\n";
            DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements, msg.str());
          }
        }
      }
    }
    // solve
    std::vector<ComplexEigenvectorType> evs;
    if (type == "default")
      evs = internal::compute_all_eigenvectors_using_eigen(matrix_.backend());
    else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type << "' is not supported, although it was reported by types()!");
    // checks
    if (check_for_inf_nan)
      for (const auto& ev : evs) {
        for (const auto& element : ev) {
          if (Common::isnan(element) || Common::isinf(element)) {
            std::stringstream msg;
            msg << "The computed eigenvectors contain inf or nan and you requested checking (see options "
                << "below)!\n"
                << "If you want to disable this check, set 'check_for_inf_nan = 0' in the options.\n\n"
                << "Those were the given options:\n\n"
                << opts;
            if (matrix_.rows() <= internal::max_size_to_print && matrix_.cols() <= internal::max_size_to_print)
              msg << "\nThis was the given matrix:\n\n" << matrix_ << "\n\nThese were the computed eigenvectors:";
            for (const auto& print_ev : evs)
              msg << "\n  " << print_ev;
            msg << "\n";
            DUNE_THROW(Exceptions::eigen_solver_failed, msg.str());
          }
        }
      }
    const double check_eigenvectors_are_real =
        opts.get("check_eigenvectors_are_real", default_opts.get<double>("check_eigenvectors_are_real"));
    if (check_eigenvectors_are_real > 0)
      for (size_t ii = 0; ii < evs.size(); ++ii)
        for (size_t jj = 0; jj < evs[ii].size(); ++jj)
          if (std::abs(evs[ii][jj].imag()) > check_eigenvectors_are_real)
            DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvectors_are_not_real_as_requested,
                       "Not all eigenvectors are purely real!\n   evs[" << ii << "][" << jj << ".imag(): "
                                                                        << evs[ii][jj].imag());
    return evs;
  } // ... eigenvectors(...)

  std::vector<RealEigenvectorType> real_eigenvectors() const
  {
    return real_eigenvectors(types().at(0));
  }

  std::vector<RealEigenvectorType> real_eigenvectors(const std::string& type) const
  {
    return real_eigenvectors(options(type));
  }

  std::vector<RealEigenvectorType> real_eigenvectors(const Common::Configuration& opts) const
  {
    Common::Configuration options_to_ensure_real_evs = opts;
    if (!options_to_ensure_real_evs.has_key("check_eigenvectors_are_real"))
      options_to_ensure_real_evs["check_eigenvectors_are_real"] = "1e-15";
    const auto evs = eigenvectors(options_to_ensure_real_evs);
    std::vector<RealEigenvectorType> ret;
    for (size_t ii = 0; ii < evs.size(); ++ii) {
      ret.emplace_back(RealEigenvectorType(evs[ii].size()));
      for (size_t jj = 0; jj < evs[ii].size(); ++jj)
        ret.back()[jj] = evs[ii][jj].real();
    }
    return ret;
  } // ... real_eigenvectors(...)

private:
  const MatrixType& matrix_;
}; // class EigenSolver<EigenDenseMatrix<S>>


//#else // HAVE_EIGEN


// template <class S>
// class EigenSolver<EigenDenseMatrix<S>>
//{
//  static_assert(AlwaysFalse<S>::value, "You are missing eigen!");
//};


//#endif // HAVE_EIGEN

} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_EIGEN_HH
