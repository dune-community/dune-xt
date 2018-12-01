// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017 - 2018)

#ifndef DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_BASE_HH
#define DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_BASE_HH

#include <algorithm>
#include <functional>
#include <memory>
#include <numeric>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/vector.hh>
#include <dune/xt/common/matrix.hh>

#include <dune/xt/la/container/common/vector/dense.hh>
#include <dune/xt/la/container/conversion.hh>
#include <dune/xt/la/container/matrix-interface.hh>
#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/matrix-inverter.hh>

namespace Dune {
namespace XT {
namespace LA {


// forward
template <class MatrixType, bool is_matrix>
class EigenSolverOptions;


namespace internal {


static inline void ensure_eigen_solver_type(const std::string& type, const std::vector<std::string>& available_types)
{
  bool contained = false;
  for (const auto& tp : available_types)
    if (type == tp)
      contained = true;
  if (!contained)
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_it_was_not_set_up_correctly,
               "Given type '" << type << "' is not one of the available types: " << available_types);
} // ... ensure_type(...)


static inline Common::Configuration default_eigen_solver_options()
{
  Common::Configuration opts;
  opts["compute_eigenvalues"] = "true";
  opts["compute_eigenvectors"] = "true";
  opts["check_for_inf_nan"] = "true";
  opts["real_tolerance"] = "1e-15"; // is only used if the respective assert_... is negative
  opts["assert_real_eigenvalues"] = "-1"; // if positive, this is the check tolerance
  opts["assert_positive_eigenvalues"] = "-1"; // if positive, this is the check tolerance
  opts["assert_negative_eigenvalues"] = "-1"; // if positive, this is the check tolerance
  opts["assert_real_eigenvectors"] = "-1"; // if positive, this is the check tolerance
  opts["assert_eigendecomposition"] = "1e-10"; // if positive, this is the check tolerance
  opts["assert_real_eigendecomposition"] = "-1"; // if positive, this is the check tolerance
  return opts;
} // ... default_eigen_solver_options(...)


/**
 * \sa default_eigen_solver_options()
 * \note If the provided options contain a subtree "matrix-inverter" that one is forwarded on eigenvector inversion.
 */
template <class MatrixImp, class FieldImp, class RealMatrixImp, class ComplexMatrixImp>
class EigenSolverBase
{
  static_assert(is_matrix<MatrixImp>::value || XT::Common::is_matrix<MatrixImp>::value, "");
  static_assert(is_matrix<RealMatrixImp>::value || XT::Common::is_matrix<RealMatrixImp>::value, "");
  static_assert(is_matrix<ComplexMatrixImp>::value || XT::Common::is_matrix<ComplexMatrixImp>::value, "");
  static_assert((is_matrix<MatrixImp>::value && is_matrix<RealMatrixImp>::value && is_matrix<ComplexMatrixImp>::value)
                    || (XT::Common::is_matrix<MatrixImp>::value && XT::Common::is_matrix<RealMatrixImp>::value
                        && XT::Common::is_matrix<ComplexMatrixImp>::value),
                "");
  using ThisType = EigenSolverBase<MatrixImp, FieldImp, RealMatrixImp, ComplexMatrixImp>;

public:
  using MatrixType = MatrixImp;
  using RealType = XT::Common::field_t<FieldImp>;
  using ComplexType = XT::Common::complex_t<RealType>;
  using RealMatrixType = RealMatrixImp;
  using ComplexMatrixType = ComplexMatrixImp;

  EigenSolverBase(const MatrixType& matrix, const std::string& type = "")
    : EigenSolverBase(matrix, EigenSolverOptions<MatrixType, true>::options(type))
  {}

  EigenSolverBase(const MatrixType& matrix, const Common::Configuration opts)
    : matrix_(matrix)
    , stored_options_(opts)
    , options_(&stored_options_)
    , computed_(false)
    , disable_checks_(options_->get<bool>("disable_checks", false))
  {
    pre_checks();
  }

  EigenSolverBase(const MatrixType& matrix, Common::Configuration* opts)
    : matrix_(matrix)
    , options_(opts)
    , computed_(false)
    , disable_checks_(options_->get<bool>("disable_checks", false))
  {
    pre_checks();
  }

  EigenSolverBase(const ThisType& other) = default;

  EigenSolverBase(ThisType&& source) = default;

  virtual ~EigenSolverBase() = default;

  const Common::Configuration& options() const
  {
    return *options_;
  }

  const MatrixType& matrix() const
  {
    return matrix_;
  }

protected:
  /**
   * \brief     Does the actual computation.
   * \attention The implementor has to fill the appropriate members!
   * \note      The implementor can assume that the given options_ contain a valid 'type' and all default keys.
   * \nte       The implementor does not need to guard against multiple calls of this method.
   */
  virtual void compute() const = 0;

public:
  const std::vector<ComplexType>& eigenvalues() const
  {
    compute_and_check();
    if (eigenvalues_)
      return *eigenvalues_;
    else if (options_->get<bool>("compute_eigenvalues"))
      DUNE_THROW(Common::Exceptions::internal_error, "The eigenvalues_ member is not filled after calling compute()!");
    else
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong,
                 "Do not call eigenvalues() if 'compute_eigenvalues' is false!\n\nThese were the given options:\n\n"
                     << *options_);
  } // ... eigenvalues(...)

  const std::vector<RealType>& real_eigenvalues() const
  {
    compute_and_check();
    if (eigenvalues_) {
      if (!real_eigenvalues_)
        compute_real_eigenvalues();
    } else if (options_->get<bool>("compute_eigenvalues"))
      DUNE_THROW(Common::Exceptions::internal_error, "The eigenvalues_ member is not filled after calling compute()!");
    else
      DUNE_THROW(
          Common::Exceptions::you_are_using_this_wrong,
          "Do not call real_eigenvalues() if 'compute_eigenvalues' is false!\n\nThese were the given options:\n\n"
              << *options_);
    assert(real_eigenvalues_ && "These have to exist after compute_real_eigenvalues()!");
    return *real_eigenvalues_;
  } // ... real_eigenvalues(...)

  std::vector<RealType> min_eigenvalues(const size_t num_evs = std::numeric_limits<size_t>::max()) const
  {
    compute_and_check();
    if (eigenvalues_) {
      if (!real_eigenvalues_)
        compute_real_eigenvalues();
    } else if (options_->get<bool>("compute_eigenvalues"))
      DUNE_THROW(Common::Exceptions::internal_error, "The eigenvalues_ member is not filled after calling compute()!");
    else
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong,
                 "Do not call min_eigenvalues() if 'compute_eigenvalues' is false!\n\nThese were the given options:\n\n"
                     << *options_);
    assert(real_eigenvalues_ && "These have to exist after compute_real_eigenvalues()!");
    std::vector<RealType> evs = *real_eigenvalues_;
    std::sort(evs.begin(), evs.end(), [](const RealType& a, const RealType& b) { return a < b; });
    evs.resize(std::min(evs.size(), num_evs));
    return evs;
  } // ... min_eigenvalues(...)

  std::vector<RealType> max_eigenvalues(const size_t num_evs = std::numeric_limits<size_t>::max()) const
  {
    compute_and_check();
    if (eigenvalues_) {
      if (!real_eigenvalues_)
        compute_real_eigenvalues();
    } else if (options_->get<bool>("compute_eigenvalues"))
      DUNE_THROW(Common::Exceptions::internal_error, "The eigenvalues_ member is not filled after calling compute()!");
    else
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong,
                 "Do not call max_eigenvalues() if 'compute_eigenvalues' is false!\n\nThese were the given options:\n\n"
                     << *options_);
    assert(real_eigenvalues_ && "These have to exist after compute_real_eigenvalues()!");
    std::vector<RealType> evs = *real_eigenvalues_;
    std::sort(evs.begin(), evs.end(), [](const RealType& a, const RealType& b) { return a > b; });
    evs.resize(std::min(evs.size(), num_evs));
    return evs;
  } // ... max_eigenvalues(...)

  const ComplexMatrixType& eigenvectors() const
  {
    compute_and_check();
    if (eigenvectors_)
      return *eigenvectors_;
    else if (options_->get<bool>("compute_eigenvectors"))
      DUNE_THROW(Common::Exceptions::internal_error, "The eigenvectors_ member is not filled after calling compute()!");
    else
      DUNE_THROW(Common::Exceptions::you_are_using_this_wrong,
                 "Do not call eigenvectors() if 'compute_eigenvectors' is false!\n\nThese were the given options:\n\n"
                     << *options_);
  } // ... eigenvectors(...)

  const ComplexMatrixType& eigenvectors_inverse() const
  {
    compute_and_check();
    if (!eigenvectors_) {
      if (options_->get<bool>("compute_eigenvectors"))
        DUNE_THROW(Common::Exceptions::internal_error,
                   "The eigenvectors_ member is not filled after calling compute()!");
      else
        DUNE_THROW(Common::Exceptions::you_are_using_this_wrong,
                   "Do not call eigenvectors_inverse() if 'compute_eigenvectors' is false!\n\nThese were the given "
                   "options:\n\n"
                       << *options_);
    }
    invert_eigenvectors();
    assert(eigenvectors_inverse_ && "This must not happen after calling invert_eigenvectors()!");
    if (!disable_checks_) {
      const double check_eigendecomposition = options_->get<double>("assert_eigendecomposition");
      if (check_eigendecomposition > 0)
        complex_eigendecomposition_helper<>::check(
            *this, check_eigendecomposition > 0 ? check_eigendecomposition : options_->get<double>("real_tolerance"));
    }
    return *eigenvectors_inverse_;
  } // ... eigenvectors(...)

  const RealMatrixType& real_eigenvectors() const
  {
    compute_and_check();
    if (eigenvectors_) {
      if (!real_eigenvectors_)
        compute_real_eigenvectors();
    } else if (options_->get<bool>("compute_eigenvectors"))
      DUNE_THROW(Common::Exceptions::internal_error, "The eigenvectors_ member is not filled after calling compute()!");
    else
      DUNE_THROW(
          Common::Exceptions::you_are_using_this_wrong,
          "Do not call real_eigenvectors() if 'compute_eigenvectors' is false!\n\nThese were the given options:\n\n"
              << *options_);
    assert(real_eigenvectors_ && "These have to exist after compute_real_eigenvectors()!");
    return *real_eigenvectors_;
  } // ... real_eigenvectors(...)

  const MatrixType& real_eigenvectors_inverse() const
  {
    compute_and_check();
    if (!real_eigenvectors_) {
      if (options_->get<double>("assert_real_eigendecomposition") > 0
          || options_->get<double>("assert_real_eigenvectors") > 0)
        DUNE_THROW(Common::Exceptions::internal_error,
                   "The real_eigenvectors_ member is not filled after calling compute()!");
      else
        DUNE_THROW(Common::Exceptions::you_are_using_this_wrong,
                   "Do not call real_eigenvectors_inverse() after providing these options:\n\n"
                       << *options_);
    }
    invert_real_eigenvectors();
    assert(real_eigenvectors_inverse_ && "This must not happen after calling invert_real_eigenvectors()!");
    if (!disable_checks_) {
      const double assert_real_eigendecomposition = options_->get<double>("assert_real_eigendecomposition");
      if (assert_real_eigendecomposition > 0.)
        assert_eigendecomposition(matrix_,
                                  *real_eigenvalues_,
                                  *real_eigenvectors_,
                                  *real_eigenvectors_inverse_,
                                  assert_real_eigendecomposition);
    }
    return *real_eigenvectors_inverse_;
  } // ... eigenvectors(...)

protected:
  void compute_and_check() const
  {
    if (!computed_) {
      compute();
      post_checks();
    }
    computed_ = true;
  }

  void pre_checks()
  {
    if (!disable_checks_) {
      // check options
      if (!options_->has_key("type"))
        DUNE_THROW(Exceptions::eigen_solver_failed_bc_it_was_not_set_up_correctly,
                   "Missing 'type' in given options:\n\n"
                       << *options_);
      internal::ensure_eigen_solver_type(options_->get<std::string>("type"),
                                         EigenSolverOptions<MatrixType, true>::types());
      const Common::Configuration default_opts =
          EigenSolverOptions<MatrixType, true>::options(options_->get<std::string>("type"));
      for (const std::string& default_key : default_opts.getValueKeys()) {
        if (!options_->has_key(default_key))
          (*options_)[default_key] = default_opts.get<std::string>(default_key);
      }
      if (options_->get<double>("real_tolerance") <= 0)
        DUNE_THROW(Exceptions::eigen_solver_failed_bc_it_was_not_set_up_correctly,
                   "It does not make sense to enforce a non-positive tolerance!");
      if (options_->get<double>("assert_positive_eigenvalues") > 0
          && options_->get<double>("assert_negative_eigenvalues") > 0)
        DUNE_THROW(Exceptions::eigen_solver_failed_bc_it_was_not_set_up_correctly,
                   "It does not make sense to assert positive and negative eigenvalues!");
      if (!options_->get<bool>("compute_eigenvalues")
          && (options_->get<double>("assert_real_eigenvalues") > 0
              || options_->get<double>("assert_positive_eigenvalues") > 0
              || options_->get<double>("assert_negative_eigenvalues") > 0
              || options_->get<double>("assert_eigendecomposition") > 0
              || options_->get<double>("assert_real_eigendecomposition") > 0))
        (*options_)["compute_eigenvalues"] = "true";
      if (options_->get<double>("assert_real_eigenvalues") <= 0
          && (options_->get<double>("assert_positive_eigenvalues") > 0
              || options_->get<double>("assert_negative_eigenvalues") > 0
              || options_->get<double>("assert_real_eigendecomposition") > 0))
        (*options_)["assert_real_eigenvalues"] = (*options_)["real_tolerance"];
      // check matrix
      check_size(matrix_);
      if (options_->get<bool>("check_for_inf_nan") && contains_inf_or_nan(matrix_)) {
        DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
                   "Given matrix contains inf or nan and you requested checking. To disable this check set "
                   "'check_for_inf_nan' to false in the options."
                       << "\n\nThese were the given options:\n\n"
                       << *options_ << "\nThis was the given matrix:\n\n"
                       << matrix_);
      }
    }
  } // ... pre_checks(...)

  void post_checks() const
  {
    if (!disable_checks_) {
      if (options_->get<bool>("compute_eigenvalues") && !eigenvalues_)
        DUNE_THROW(Common::Exceptions::internal_error,
                   "The eigenvalues_ member is not filled after calling compute()!");
      if (options_->get<bool>("compute_eigenvectors") && !eigenvectors_)
        DUNE_THROW(Common::Exceptions::internal_error,
                   "The eigenvectors_ member is not filled after calling compute()!");
      if (options_->get<bool>("check_for_inf_nan")) {
        if (eigenvalues_ && contains_inf_or_nan(*eigenvalues_))
          DUNE_THROW(Exceptions::eigen_solver_failed_bc_result_contained_inf_or_nan,
                     "Computed eigenvalues contain inf or nan and you requested checking. To disable this check set "
                     "'check_for_inf_nan' to false in the options."
                         << "\n\nThese were the given options:\n\n"
                         << *options_ << "\nThese are the computed eigenvalues:\n\n"
                         << *eigenvalues_);
        if (eigenvectors_ && contains_inf_or_nan(*eigenvectors_))
          DUNE_THROW(Exceptions::eigen_solver_failed_bc_result_contained_inf_or_nan,
                     "Computed eigenvectors contain inf or nan and you requested checking. To disable this check set "
                     "'check_for_inf_nan' to false in the options."
                         << "\n\nThese were the given options:\n\n"
                         << *options_ << "\nThese are the computed eigenvectors:\n\n"
                         << *eigenvectors_);
      }
      const double assert_real_eigenvalues = options_->get<double>("assert_real_eigenvalues");
      const double assert_positive_eigenvalues = options_->get<double>("assert_positive_eigenvalues");
      const double assert_negative_eigenvalues = options_->get<double>("assert_negative_eigenvalues");
      const double check_real_eigendecomposition = options_->get<double>("assert_real_eigendecomposition");
      if (assert_real_eigenvalues > 0 || assert_positive_eigenvalues > 0 || assert_negative_eigenvalues > 0
          || check_real_eigendecomposition > 0)
        compute_real_eigenvalues();
      if (assert_positive_eigenvalues > 0) {
        assert(real_eigenvalues_ && "This must not happen after compute_real_eigenvalues()!");
        for (const auto& ev : *real_eigenvalues_) {
          if (ev < assert_positive_eigenvalues)
            DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvalues_are_not_positive_as_requested,
                       "These were the given options:\n\n"
                           << *options_ << "\nThese are the computed eigenvectors:\n\n"
                           << *eigenvectors_);
        }
      }
      if (assert_negative_eigenvalues > 0) {
        assert(real_eigenvalues_ && "This must not happen after compute_real_eigenvalues()!");
        for (const auto& ev : *real_eigenvalues_) {
          if (ev > -1 * assert_negative_eigenvalues)
            DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvalues_are_not_negative_as_requested,
                       "These were the given options:\n\n"
                           << *options_ << "\nThese are the computed eigenvectors:\n\n"
                           << *eigenvectors_);
        }
      }
      if (options_->get<double>("assert_real_eigenvectors") > 0 || check_real_eigendecomposition > 0)
        compute_real_eigenvectors();
      const double check_eigendecomposition = options_->get<double>("assert_eigendecomposition");
      if (check_eigendecomposition > 0)
        complex_eigendecomposition_helper<>::check(*this, check_eigendecomposition);
      if (check_real_eigendecomposition > 0) {
        invert_real_eigenvectors();
        assert_eigendecomposition(matrix_,
                                  *real_eigenvalues_,
                                  *real_eigenvectors_,
                                  *real_eigenvectors_inverse_,
                                  check_real_eigendecomposition);
      }
    }
  } // ... post_checks(...)

  void compute_real_eigenvalues() const
  {
    assert(eigenvalues_ && "This should not happen!");
    if (!real_eigenvalues_) {
      real_eigenvalues_ = std::make_unique<std::vector<RealType>>(eigenvalues_->size());
      for (size_t ii = 0; ii < eigenvalues_->size(); ++ii)
        (*real_eigenvalues_)[ii] = (*eigenvalues_)[ii].real();

      if (!disable_checks_) {
        const double assert_real_eigenvalues = options_->get<double>("assert_real_eigenvalues");
        const double tolerance =
            (assert_real_eigenvalues > 0) ? assert_real_eigenvalues : options_->get<double>("real_tolerance");
        for (size_t ii = 0; ii < eigenvalues_->size(); ++ii) {
          if (std::abs((*eigenvalues_)[ii].imag()) > tolerance)
            DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvalues_are_not_real_as_requested,
                       "These were the given options:\n\n"
                           << *options_ << "\nThese are the computed eigenvalues:\n\n"
                           << *eigenvalues_);
        } // ii
      } // if (!disable_checks_)
    } // if (!real_eigenvalues_)
  } // ... compute_real_eigenvalues(...)

  template <bool is_common_matrix = XT::Common::is_matrix<MatrixType>::value, class T = MatrixType>
  struct real_eigenvectors_helper
  {};

  template <class T>
  struct real_eigenvectors_helper<true, T>
  {
    static void compute(const ThisType& self, const double& tolerance)
    {
      using RM = XT::Common::MatrixAbstraction<RealMatrixType>;
      using CM = XT::Common::MatrixAbstraction<ComplexMatrixType>;
      const size_t rows = CM::rows(*self.eigenvectors_);
      const size_t cols = CM::cols(*self.eigenvectors_);
      self.real_eigenvectors_ = std::make_unique<RealMatrixType>(RM::create(rows, cols));
      bool is_complex = false;
      for (size_t ii = 0; ii < rows; ++ii) {
        for (size_t jj = 0; jj < cols; ++jj) {
          const auto complex_value = CM::get_entry(*self.eigenvectors_, ii, jj);
          if (std::abs(complex_value.imag()) > tolerance) {
            is_complex = true;
            ii = rows;
            break;
          }
          RM::set_entry(*self.real_eigenvectors_, ii, jj, complex_value.real());
        } // jj
      } // ii

      if (is_complex) {
        // try to get real eigenvectors from the complex ones. If both the matrix and the eigenvalues are real, the
        // eigenvectors can also be chosen real. If there is a imaginary eigenvector, both real and imaginary part are
        // eigenvectors (if non-zero) to the same eigenvalue. So to get real eigenvectors, sort the eigenvectors by
        // eigenvalues and, separately for each eigenvalue, perform a Gram-Schmidt process with all real and imaginary
        // parts of the eigenvectors
        self.compute_real_eigenvalues();

        // form groups of equal eigenvalues
        struct Cmp
        {
          bool operator()(const RealType& a, const RealType& b) const
          {
            return XT::Common::FloatCmp::lt(a, b);
          }
        };
        std::vector<std::vector<size_t>> eigenvalue_groups;
        std::vector<size_t> eigenvalue_multiplicity;
        std::set<RealType, Cmp> eigenvalues_done;
        for (size_t jj = 0; jj < rows; ++jj) {
          const auto curr_eigenvalue = (*self.real_eigenvalues_)[jj];
          if (!eigenvalues_done.count(curr_eigenvalue)) {
            std::vector<size_t> curr_group;
            curr_group.push_back(jj);
            eigenvalue_multiplicity.push_back(1);
            for (size_t kk = jj + 1; kk < rows; ++kk) {
              if (XT::Common::FloatCmp::eq(curr_eigenvalue, (*self.real_eigenvalues_)[kk])) {
                curr_group.push_back(kk);
                ++(eigenvalue_multiplicity.back());
              }
            } // kk
            eigenvalue_groups.push_back(curr_group);
            eigenvalues_done.insert(curr_eigenvalue);
          }
        } // jj

        // For each eigenvalue, calculate a orthogonal basis of the n-dim real eigenspace from the 2n real
        // and imaginary parts of the complex eigenvectors
        for (size_t kk = 0; kk < eigenvalue_groups.size(); ++kk) {
          const auto& group = eigenvalue_groups[kk];
          typedef typename XT::LA::CommonDenseVector<RealType> RealVectorType;
          std::vector<RealVectorType> input_vectors(2 * eigenvalue_multiplicity[kk], RealVectorType(rows, 0.));
          size_t index = 0;
          for (const auto& jj : group) {
            for (size_t ll = 0; ll < cols; ++ll) {
              input_vectors[index][ll] = CM::get_entry(*self.eigenvectors_, ll, jj).real();
              input_vectors[index + 1][ll] = CM::get_entry(*self.eigenvectors_, ll, jj).imag();
            }
            index += 2;
          } // jj

          // orthonormalize
          for (size_t ii = 0; ii < input_vectors.size(); ++ii) {
            auto& v_i = input_vectors[ii];
            for (size_t jj = 0; jj < ii; ++jj) {
              const auto& v_j = input_vectors[jj];
              const auto vj_vj = v_j.dot(v_j);
              if (XT::Common::FloatCmp::eq(vj_vj, 0.))
                continue;
              const auto vj_vi = v_j.dot(v_i);
              for (size_t rr = 0; rr < rows; ++rr)
                v_i[rr] -= vj_vi / vj_vj * v_j[rr];
            } // jj
            RealType l2_norm = std::sqrt(std::accumulate(
                v_i.begin(), v_i.end(), 0., [](const RealType& a, const RealType& b) { return a + b * b; }));
            if (XT::Common::FloatCmp::ne(l2_norm, 0.))
              v_i *= 1. / l2_norm;
          } // ii
          // copy eigenvectors back to eigenvectors matrix
          index = 0;
          for (size_t ii = 0; ii < input_vectors.size(); ++ii) {
            if (XT::Common::FloatCmp::ne(input_vectors[ii], RealVectorType(rows, 0.))) {
              if (index >= eigenvalue_multiplicity[kk]) {
                DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvectors_are_not_real_as_requested,
                           "Eigenvectors are complex and calculating real eigenvectors failed!"
                               << "These were the given options:\n\n"
                               << *self.options_ << "\n\nThis was the given matrix: " << std::setprecision(17)
                               << self.matrix_ << "\nThese are the computed eigenvectors:\n\n"
                               << std::setprecision(17) << *self.eigenvectors_);
              }
              for (size_t rr = 0; rr < rows; ++rr)
                RM::set_entry(*self.real_eigenvectors_, rr, group[index], input_vectors[ii].get_entry(rr));
              index++;
            } // if (input_vectors[ii] != 0)
          } // ii
          if (index < eigenvalue_multiplicity[kk]) {
            DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvectors_are_not_real_as_requested,
                       "Eigenvectors are complex and calculating real eigenvectors failed!"
                           << "These were the given options:\n\n"
                           << *self.options_ << "\n\nThis was the given matrix: " << std::setprecision(17)
                           << self.matrix_ << "\nThese are the computed eigenvectors:\n\n"
                           << std::setprecision(17) << *self.eigenvectors_);
          }
        } // kk
      } // if(is_complex)
    } // static void compute(...)
  }; // real_eigenvectors_helper<true, ...>

  template <class T>
  struct real_eigenvectors_helper<false, T>
  {
    static void compute(ThisType& self, const double& tolerance)
    {
      if (RealMatrixType::sparse) {
        const size_t rows = self.eigenvectors_->rows();
        const size_t cols = self.eigenvectors_->cols();
        const auto pattern = self.eigenvectors_->pattern();
        self.real_eigenvectors_ = std::make_unique<RealMatrixType>(rows, cols, pattern);
        for (size_t ii = 0; ii < rows; ++ii)
          for (size_t jj : pattern.inner(ii)) {
            const auto complex_value = self.eigenvectors_->get_entry(ii, jj);
            if (std::abs(complex_value.imag()) > tolerance)
              DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvectors_are_not_real_as_requested,
                         "These were the given options:\n\n"
                             << *self.options_ << "\nThese are the computed eigenvectors:\n\n"
                             << std::setprecision(17) << *self.eigenvectors_);
            self.real_eigenvectors_->set_entry(ii, jj, complex_value.real());
          }
      } else {
        const size_t rows = self.eigenvectors_->rows();
        const size_t cols = self.eigenvectors_->cols();
        self.real_eigenvectors_ = std::make_unique<RealMatrixType>(rows, cols);
        for (size_t ii = 0; ii < rows; ++ii)
          for (size_t jj = 0; jj < cols; ++jj) {
            const auto complex_value = self.eigenvectors_->get_entry(ii, jj);
            if (std::abs(complex_value.imag()) > tolerance)
              DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvectors_are_not_real_as_requested,
                         "These were the given options:\n\n"
                             << *self.options_ << "\nThese are the computed eigenvectors:\n\n"
                             << std::setprecision(17) << *self.eigenvectors_);
            self.real_eigenvectors_->set_entry(ii, jj, complex_value.real());
          }
      }
    }
  }; // real_eigenvectors_helper<false, ...>

  void compute_real_eigenvectors() const
  {
    assert(eigenvectors_ && "This should not happen!");
    if (!real_eigenvectors_) {
      const double assert_real_eigenvectors = options_->get<double>("assert_real_eigenvectors");
      const double tolerance =
          (assert_real_eigenvectors > 0) ? assert_real_eigenvectors : options_->get<double>("real_tolerance");
      real_eigenvectors_helper<>::compute(*this, tolerance);
    }
  }

  void invert_eigenvectors() const
  {
    assert(eigenvectors_ && "This must not happen when you call this function!");
    try {
      if (options_->has_sub("matrix-inverter")) {
        eigenvectors_inverse_ =
            std::make_unique<ComplexMatrixType>(invert_matrix(*eigenvectors_, options_->sub("matrix-inverter")));
      } else
        eigenvectors_inverse_ = std::make_unique<ComplexMatrixType>(invert_matrix(*eigenvectors_));
    } catch (const Exceptions::matrix_invert_failed& ee) {
      DUNE_THROW(Exceptions::eigen_solver_failed,
                 "The computed matrix of eigenvectors is not invertible!"
                     << "\n\nmatrix = " << std::setprecision(17) << matrix_ << "\n\noptions: " << *options_
                     << "\n\neigenvectors = " << std::setprecision(17) << *eigenvectors_
                     << "\n\nThis was the original error: " << ee.what());
    }
  } // ... invert_eigenvectors(...)

  void invert_real_eigenvectors() const
  {
    assert(real_eigenvectors_ && "This must not happen when you call this function!");
    try {
      if (options_->has_sub("matrix-inverter")) {
        real_eigenvectors_inverse_ =
            std::make_unique<MatrixType>(invert_matrix(*real_eigenvectors_, options_->sub("matrix-inverter")));
      } else
        real_eigenvectors_inverse_ = std::make_unique<MatrixType>(invert_matrix(*real_eigenvectors_));
    } catch (const Exceptions::matrix_invert_failed& ee) {
      DUNE_THROW(Exceptions::eigen_solver_failed,
                 "The computed matrix of real eigenvectors is not invertible!"
                     << "\n\nmatrix = " << std::setprecision(17) << matrix_ << "\n\noptions: " << *options_
                     << "\n\nreal_eigenvectors = " << std::setprecision(17) << *real_eigenvectors_
                     << "\n\nThis was the original error: " << ee.what());
    }
  } // ... invert_real_eigenvectors(...)

  template <class A, class B, class C, class D>
  void assert_eigendecomposition(const std::unique_ptr<A>& mat,
                                 const B& eigenvalues,
                                 const C& eigenvectors,
                                 const D& eigenvectors_inverse,
                                 const double& tolerance) const
  {
    assert_eigendecomposition(*mat, eigenvalues, eigenvectors, eigenvectors_inverse, tolerance);
  }

  template <class A, class B, class C, class D>
  void assert_eigendecomposition(const A& mat,
                                 const B& eigenvalues,
                                 const C& eigenvectors,
                                 const D& eigenvectors_inverse,
                                 const double& tolerance) const
  {
    using M = Common::MatrixAbstraction<C>;
    const size_t rows = Common::get_matrix_rows(mat);
    const size_t cols = Common::get_matrix_cols(mat);
    auto eigenvalue_matrix = M::create(rows, cols, typename M::ScalarType(0.));
    for (size_t ii = 0; ii < rows; ++ii)
      Common::set_matrix_entry(eigenvalue_matrix, ii, ii, eigenvalues[ii]);
    const auto decomposition_error = eigenvectors * eigenvalue_matrix * eigenvectors_inverse - mat;
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        if (std::abs(Common::get_matrix_entry(decomposition_error, ii, jj)) > tolerance)
          DUNE_THROW(Exceptions::eigen_solver_failed_bc_result_is_not_an_eigendecomposition,
                     "\n\nmatrix = " << std::setprecision(17) << matrix_ << "\n\noptions: " << *options_
                                     << "\n\neigenvalues (lambda)= " << std::setprecision(17) << eigenvalues
                                     << "\n\neigenvectors (T) = " << std::setprecision(17) << eigenvectors
                                     << "\n\n(T * (lambda * T^-1)) - matrix = " << decomposition_error);
  } // ... assert_eigendecomposition(...)

  template <bool upcast_required = !std::is_same<MatrixType, ComplexMatrixType>::value, bool anything = true>
  struct complex_eigendecomposition_helper;

  template <bool anything>
  struct complex_eigendecomposition_helper<true, anything>
  {
    static void check(const ThisType& self, const double& tolerance)
    {
      self.invert_eigenvectors();
      self.assert_eigendecomposition(Dune::XT::LA::convert_to<ComplexMatrixType>(self.matrix_),
                                     *self.eigenvalues_,
                                     *self.eigenvectors_,
                                     *self.eigenvectors_inverse_,
                                     tolerance);
    }
  };

  template <bool anything>
  struct complex_eigendecomposition_helper<false, anything>
  {
    static void check(const ThisType& self, const double& tolerance)
    {
      self.invert_eigenvectors();
      self.assert_eigendecomposition(
          self.matrix_, *self.eigenvalues_, *self.eigenvectors_, *self.eigenvectors_inverse_, tolerance);
    }
  };

  template <class M>
  void check_size(const MatrixInterface<M>& mat) const
  {
    if (mat.rows() != mat.cols())
      DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
                 "Matrix has to be square, is " << mat.rows() << "x" << mat.cols() << "!");
  }

  template <class M>
  typename std::enable_if<XT::Common::is_matrix<M>::value && !is_matrix<M>::value, void>::type
  check_size(const M& mat) const
  {
    using Mat = XT::Common::MatrixAbstraction<M>;
    if (Mat::rows(mat) != Mat::cols(mat))
      DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
                 "Matrix has to be square, is " << Mat::rows(mat) << "x" << Mat::cols(mat) << "!");
  }

  template <class T>
  bool contains_inf_or_nan(const std::vector<T>& vec) const
  {
    for (const auto& element : vec)
      if (XT::Common::isinf(element) || XT::Common::isnan(element))
        return true;
    return false;
  }

  template <class M>
  bool contains_inf_or_nan(const MatrixInterface<M>& mat) const
  {
    return !mat.valid();
  }

  template <class M>
  typename std::enable_if<XT::Common::is_matrix<M>::value && !is_matrix<M>::value, bool>::type
  contains_inf_or_nan(const M& mat) const
  {
    using Mat = XT::Common::MatrixAbstraction<M>;
    for (size_t ii = 0; ii < Mat::rows(mat); ++ii)
      for (size_t jj = 0; jj < Mat::cols(mat); ++jj) {
        const auto value = Mat::get_entry(mat, ii, jj);
        if (XT::Common::isinf(value) || XT::Common::isnan(value))
          return true;
      }
    return false;
  } // ... contains_inf_or_nan(...)

  const MatrixType& matrix_;
  mutable Common::Configuration stored_options_;
  mutable Common::Configuration* options_;
  mutable bool computed_;
  mutable std::unique_ptr<std::vector<ComplexType>> eigenvalues_;
  mutable std::unique_ptr<std::vector<RealType>> real_eigenvalues_;
  mutable std::unique_ptr<ComplexMatrixType> eigenvectors_;
  mutable std::unique_ptr<RealMatrixType> real_eigenvectors_;
  mutable std::unique_ptr<ComplexMatrixType> eigenvectors_inverse_;
  mutable std::unique_ptr<RealMatrixType> real_eigenvectors_inverse_;
  const bool disable_checks_;
}; // class EigenSolverBase


} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_BASE_HH
