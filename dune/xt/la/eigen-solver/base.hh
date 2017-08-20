// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_LA_EIGEN_SOLVER_BASE_HH
#define DUNE_XT_LA_EIGEN_SOLVER_BASE_HH

#include <algorithm>
#include <functional>

#include <dune/xt/la/solver.hh>

namespace Dune {
namespace XT {
namespace LA {


template <class Traits>
class EigenSolverBase : protected internal::SolverUtils
{
public:
  typedef typename Traits::MatrixType MatrixType;
  typedef typename Traits::RealType RealType;
  typedef typename Traits::ComplexType ComplexType;
  typedef typename Traits::RealMatrixType RealMatrixType;
  typedef typename Traits::ComplexMatrixType ComplexMatrixType;
  typedef typename Traits::RealVectorType RealVectorType;
  typedef typename Traits::ComplexVectorType ComplexVectorType;
  typedef typename Traits::derived_type derived_type;
  typedef XT::Common::MatrixAbstraction<MatrixType> MatrixAbstractionType;
  typedef XT::Common::MatrixAbstraction<RealMatrixType> RealMatrixAbstractionType;
  typedef XT::Common::MatrixAbstraction<ComplexMatrixType> ComplexMatrixAbstractionType;

  EigenSolverBase(const MatrixType& matrix)
    : matrix_(matrix)
  {
  }

  virtual ~EigenSolverBase()
  {
  }

  static std::vector<std::string> types()
  {
    return derived_type::types();
  }

  static Common::Configuration options(const std::string type = "")
  {
    return derived_type::options(type);
  }

  virtual void get_eigenvalues(std::vector<ComplexType>& evs, const std::string& type) const = 0;
  virtual void get_eigenvectors(std::vector<ComplexVectorType>& evs, const std::string& type) const = 0;

  virtual std::vector<std::complex<RealType>> eigenvalues() const
  {
    return eigenvalues(types().at(0));
  }

  virtual std::vector<std::complex<RealType>> eigenvalues(const std::string type) const
  {
    return eigenvalues(options(type));
  }

  virtual std::string parse_opts(const XT::Common::Configuration& opts) const
  {
    if (!opts.has_key("type"))
      DUNE_THROW(Common::Exceptions::configuration_error,
                 "Given options (see below) need to have at least the key 'type' set!\n\n"
                     << opts);
    const auto type = opts.get<std::string>("type");
    internal::SolverUtils::check_given(type, types());
    return type;
  }

  virtual std::vector<std::complex<RealType>> eigenvalues(const Common::Configuration& opts) const
  {
    const auto type = parse_opts(opts);
    // check
    const Common::Configuration default_opts = options(type);
    check_input(opts, default_opts);
    // solve
    std::vector<ComplexType> evs;
    get_eigenvalues(evs, type);
    // check
    check_eigenvalues(evs, opts, default_opts);
    // return
    return evs;
  } // ... eigenvalues(...)

  template <class... Args>
  std::vector<RealType> min_eigenvalues(Args&&... args) const
  {
    return sorted_eigenvalues(true, std::forward<Args>(args)...);
  } // ... min_eigenvalues(...)

  template <class... Args>
  std::vector<RealType> max_eigenvalues(Args&&... args) const
  {
    return sorted_eigenvalues(false, std::forward<Args>(args)...);
  } // ... min_eigenvalues(...)

  virtual std::vector<ComplexVectorType> eigenvectors() const
  {
    return eigenvectors(types().at(0));
  }

  virtual std::vector<ComplexVectorType> eigenvectors(const std::string& type) const
  {
    return eigenvectors(options(type));
  }

  virtual std::vector<ComplexVectorType> eigenvectors(const Common::Configuration& opts) const
  {
    const auto type = parse_opts(opts);
    // check
    const Common::Configuration default_opts = options(type);
    check_input(opts, default_opts);
    // solve
    const size_t dim = RealMatrixAbstractionType::rows(matrix_);
    std::vector<ComplexVectorType> evs(dim, XT::Common::VectorAbstraction<ComplexVectorType>::create(dim, 0.));
    get_eigenvectors(evs, type);
    // check
    check_eigenvectors(evs, opts, default_opts);
    // return
    return evs;
  } // ... eigenvectors(...)

  virtual std::vector<RealVectorType> real_eigenvectors() const
  {
    return real_eigenvectors(types().at(0));
  }

  virtual std::vector<RealVectorType> real_eigenvectors(const std::string& type) const
  {
    return real_eigenvectors(options(type));
  }

  virtual std::vector<RealVectorType> real_eigenvectors(const Common::Configuration& opts) const
  {
    Common::Configuration options_to_ensure_real_evs = opts;
    if (!options_to_ensure_real_evs.has_key("check_eigenvectors_are_real"))
      options_to_ensure_real_evs["check_eigenvectors_are_real"] = "1e-15";
    const auto evs = eigenvectors(options_to_ensure_real_evs);
    std::vector<RealVectorType> ret(evs.size(), RealVectorType(evs[0].size()));
    for (size_t ii = 0; ii < evs.size(); ++ii)
      for (size_t jj = 0; jj < evs[ii].size(); ++jj)
        ret[ii][jj] = evs[ii][jj].real();
    return ret;
  } // ... real_eigenvectors(...)

  virtual std::shared_ptr<ComplexMatrixType> eigenvectors_as_matrix() const
  {
    return eigenvectors_as_matrix(types().at(0));
  }

  virtual std::shared_ptr<ComplexMatrixType> eigenvectors_as_matrix(const std::string& type) const
  {
    return eigenvectors_as_matrix(options(type));
  }

  virtual std::shared_ptr<ComplexMatrixType> eigenvectors_as_matrix(const Common::Configuration& opts) const
  {
    const size_t rows = MatrixAbstractionType::rows(matrix_);
    const size_t cols = MatrixAbstractionType::cols(matrix_);
    const auto evs = eigenvectors(opts);
    auto ret = std::make_shared<ComplexMatrixType>(ComplexMatrixAbstractionType::create(rows, cols));
    for (size_t ii = 0; ii < MatrixAbstractionType::rows(matrix_); ++ii)
      for (size_t jj = 0; jj < MatrixAbstractionType::cols(matrix_); ++jj)
        ComplexMatrixAbstractionType::set_entry(*ret, ii, jj, evs[jj][ii]);
    return ret;
  } // ... eigenvectors_as_matrix(...)

  virtual std::shared_ptr<RealMatrixType> real_eigenvectors_as_matrix() const
  {
    return real_eigenvectors_as_matrix(types().at(0));
  }

  virtual std::shared_ptr<RealMatrixType> real_eigenvectors_as_matrix(const std::string& type) const
  {
    return real_eigenvectors_as_matrix(options(type));
  }

  virtual std::shared_ptr<RealMatrixType> real_eigenvectors_as_matrix(const Common::Configuration& opts) const
  {
    const size_t rows = MatrixAbstractionType::rows(matrix_);
    const size_t cols = MatrixAbstractionType::cols(matrix_);
    const auto evs = real_eigenvectors(opts);
    auto ret = std::make_shared<RealMatrixType>(RealMatrixAbstractionType::create(rows, cols));
    for (size_t ii = 0; ii < MatrixAbstractionType::rows(matrix_); ++ii)
      for (size_t jj = 0; jj < MatrixAbstractionType::cols(matrix_); ++jj)
        RealMatrixAbstractionType::set_entry(*ret, ii, jj, evs[jj][ii]);
    return ret;
  } // ... real_eigenvectors_as_matrix(...)

protected:
  virtual void check_input(const Common::Configuration& opts, const Common::Configuration& default_opts) const
  {
    const size_t rows = MatrixAbstractionType::rows(matrix_);
    const size_t cols = MatrixAbstractionType::cols(matrix_);
    if (rows != cols)
      DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
                 "Matrix has to be square!\n   matrix.rows() = " << rows << "\n   matrix.cols()" << cols);
    const bool check_for_inf_nan = opts.get("check_for_inf_nan", default_opts.get<bool>("check_for_inf_nan"));
    if (check_for_inf_nan) {
      for (size_t ii = 0; ii < rows; ++ii) {
        for (size_t jj = 0; jj < cols; ++jj) {
          const auto val = MatrixAbstractionType::get_entry(matrix_, ii, jj);
          if (Common::isnan(val) || Common::isinf(val)) {
            std::stringstream msg;
            msg << "Given matrix contains inf or nan and you requested checking (see options below)!\n"
                << "If you want to disable this check, set 'check_for_inf_nan = 0' in the options.\n\n"
                << "Those were the given options:\n\n"
                << opts;
            if (rows <= internal::max_size_to_print)
              msg << "\nThis was the given matrix:\n\n" << XT::Common::to_string(matrix_) << "\n";
            DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements, msg.str());
          } // if (isnan || isinf)
        } // jj
      } // ii
    } // if (check_for_inf_nan)
  } // ... check_input(...)

  virtual void check_eigenvalues(const std::vector<ComplexType>& evs,
                                 const Common::Configuration& opts,
                                 const Common::Configuration& default_opts) const
  {
    const bool check_for_inf_nan = opts.get("check_for_inf_nan", default_opts.get<bool>("check_for_inf_nan"));
    if (check_for_inf_nan) {
      const size_t rows = MatrixAbstractionType::rows(matrix_);
      for (const auto& ev : evs) {
        if (Common::isnan(ev) || Common::isinf(ev)) {
          std::stringstream msg;
          msg << "The computed eigenvalues contain inf or nan and you requested checking (see options "
              << "below)!\n"
              << "If you want to disable this check, set 'check_for_inf_nan = 0' in the options.\n\n"
              << "Those were the given options:\n\n"
              << opts;
          if (rows <= internal::max_size_to_print)
            msg << "\nThis was the given matrix:\n\n"
                << XT::Common::to_string(matrix_) << "\n\nThese were the computed eigenvalues:";
          for (const auto& print_ev : evs)
            msg << "\n  " << print_ev;
          msg << "\n";
          DUNE_THROW(Exceptions::eigen_solver_failed, msg.str());
        } // if (isnan || isinf)
      } // loop over evs
    } // if (check_for_inf_nan)
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
        if (evs[ii].real() < check_evs_are_positive)
          DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvalues_are_not_positive_as_requested,
                     "Not all eigenvalues are positive!\n   evs[" << ii << "].real(): " << evs[ii].real());

  } // ... check_eigenvalues(...);

  virtual void check_eigenvectors(const std::vector<ComplexVectorType>& evs,
                                  const Common::Configuration& opts,
                                  const Common::Configuration& default_opts) const
  {
    const bool check_for_inf_nan = opts.get("check_for_inf_nan", default_opts.get<bool>("check_for_inf_nan"));
    if (check_for_inf_nan) {
      const size_t rows = MatrixAbstractionType::rows(matrix_);
      for (const auto& ev : evs) {
        for (const auto& element : ev) {
          if (Common::isnan(element) || Common::isinf(element)) {
            std::stringstream msg;
            msg << "The computed eigenvectors contain inf or nan and you requested checking (see options "
                << "below)!\n"
                << "If you want to disable this check, set 'check_for_inf_nan = 0' in the options.\n\n"
                << "Those were the given options:\n\n"
                << opts;
            if (rows <= internal::max_size_to_print)
              msg << "\nThis was the given matrix:\n\n" << matrix_ << "\n\nThese were the computed eigenvectors:";
            for (const auto& print_ev : evs)
              msg << "\n  " << print_ev;
            msg << "\n";
            DUNE_THROW(Exceptions::eigen_solver_failed, msg.str());
          } // if (isnan || isinf)
        } // element
      } // ev
    } // if (check_for_inf_nan)
    const double check_eigenvectors_are_real =
        opts.get("check_eigenvectors_are_real", default_opts.get<double>("check_eigenvectors_are_real"));
    if (check_eigenvectors_are_real > 0)
      for (size_t ii = 0; ii < evs.size(); ++ii)
        for (size_t jj = 0; jj < evs[ii].size(); ++jj)
          if (std::abs(evs[ii][jj].imag()) > check_eigenvectors_are_real)
            DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvectors_are_not_real_as_requested,
                       "Not all eigenvectors are purely real!\n   evs[" << ii << "][" << jj << ".imag(): "
                                                                        << evs[ii][jj].imag());
  } // ... check_eigenvectors(...)

  std::vector<RealType> sorted_eigenvalues(bool ascending,
                                           const size_t num_evs = std::numeric_limits<size_t>::max()) const
  {
    return sorted_eigenvalues(ascending, types().at(0), num_evs);
  }

  std::vector<RealType> sorted_eigenvalues(bool ascending,
                                           const std::string type,
                                           const size_t num_evs = std::numeric_limits<size_t>::max()) const
  {
    return sorted_eigenvalues(ascending, options(type), num_evs);
  }

  std::vector<RealType> sorted_eigenvalues(bool ascending,
                                           const Common::Configuration& opts,
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
                 "Computing the " + std::string(ascending ? "minimum" : "maximum")
                         + " eigenvalue does not make sense "
                           "for complex eigenvalues!\nYou can configure the check by adapting check_evs_are_real in "
                           "the options "
                           "(see below).\n This was the original error:\n\n"
                     << ee.what()
                     << "\n\nThese were the options:\n"
                     << options);
    }
    std::sort(evs.begin(), evs.end(), [ascending](RealType a, RealType b) { return ascending ? a < b : b < a; });
    evs.resize(std::min(evs.size(), num_evs));
    return evs;
  } // ... sorted_eigenvalues(...)

  const MatrixType& matrix_;
}; // class EigenSolverBase<...>


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_BASE_HH
