// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2017 - 2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_XT_LA_MATRIX_INVERTER_BASE_HH
#define DUNE_XT_LA_MATRIX_INVERTER_BASE_HH

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/math.hh>
#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/type_traits.hh>

#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/la/container/matrix-interface.hh>
#include <dune/xt/la/exceptions.hh>

namespace Dune {
namespace XT {
namespace LA {


// forward
template <class MatrixType, bool is_matrix>
class MatrixInverterOptions;


namespace internal {


static inline void ensure_matrix_inverter_type(const std::string& type, const std::vector<std::string>& available_types)
{
  bool contained = false;
  for (const auto& tp : available_types)
    if (type == tp)
      contained = true;
  if (!contained)
    DUNE_THROW(Exceptions::matrix_invert_failed_bc_it_was_not_set_up_correctly,
               "Given type '" << type << "' is not one of the available types: " << available_types);
} // ... ensure_type(...)


static inline Common::Configuration default_matrix_inverter_options()
{
  Common::Configuration opts;
  opts["delay_computation"] = "false";
  opts["check_for_inf_nan"] = "true";
  opts["post_check_is_left_inverse"] = "1e-10";
  opts["post_check_is_right_inverse"] = "1e-10";
  return opts;
}


/**
 * \todo Refactor just like EigenSolverBase (no delay_computation)!
 */
template <class MatrixImp>
class MatrixInverterBase
{
  static_assert(is_matrix<MatrixImp>::value || XT::Common::is_matrix<MatrixImp>::value, "");

public:
  using MatrixType = MatrixImp;

  /**
   * \attention The implementor has to call compute() in the ctor if (delay_computation == true).
   */
  MatrixInverterBase(const MatrixType& matrix, const std::string& type = "")
    : matrix_(matrix)
    , options_(MatrixInverterOptions<MatrixType, true>::options(type))
    , inverse_(nullptr)
  {
    pre_checks();
  }

  /**
   * \attention The implementor has to call compute() in the ctor if (delay_computation == true).
   */
  MatrixInverterBase(const MatrixType& matrix, const Common::Configuration opts)
    : matrix_(matrix)
    , options_(opts)
    , inverse_(nullptr)
  {
    pre_checks();
  }

  virtual ~MatrixInverterBase() = default;

  const Common::Configuration& options() const
  {
    return options_;
  }

  const MatrixType& matrix() const
  {
    return matrix_;
  }

  /**
   * \brief     Does the actual computation.
   * \attention The implementor has to fill inverse_!
   * \attention The implementor is supposed to call post_check() at the end of the computation!
   * \note      The implementor can assume that the given options_ contain a valid 'type'.
   */
  virtual void compute() = 0;

  const MatrixType& inverse() const
  {
    if (!inverse_)
      compute();
    if (!inverse_)
      DUNE_THROW(Common::Exceptions::internal_error, "The inverse_ member is not filled after calling compute()!");
    return *inverse_;
  }

  MatrixType inverse()
  {
    if (!inverse_)
      compute();
    if (!inverse_)
      DUNE_THROW(Common::Exceptions::internal_error, "The inverse_ member is not filled after calling compute()!");
    return *inverse_;
  }

protected:
  void pre_checks()
  {
    // check options
    if (!options_.has_key("type"))
      DUNE_THROW(Exceptions::matrix_invert_failed_bc_it_was_not_set_up_correctly,
                 "Missing 'type' in given options!"
                     << "\n\nThese were the given options:\n\n"
                     << options_);
    internal::ensure_matrix_inverter_type(options_.get<std::string>("type"),
                                          MatrixInverterOptions<MatrixType, true>::types());
    const Common::Configuration default_opts =
        MatrixInverterOptions<MatrixType, true>::options(options_.get<std::string>("type"));
    // check matrix
    if (options_.get("check_for_inf_nan", default_opts.get<bool>("check_for_inf_nan"))) {
      if (contains_inf_or_nan(matrix_))
        DUNE_THROW(Exceptions::matrix_invert_failed_bc_data_did_not_fulfill_requirements,
                   "Given matrix contains inf or nan and you requested checking. To disable this check set "
                   "'check_for_inf_nan' to false in the options."
                       << "\n\nThese were the given options:\n\n"
                       << options_ << "\nThis was the given matrix:\n\n"
                       << matrix_);
    }
  } // ... pre_checks(...)

  void post_checks() const
  {
    if (!inverse_)
      DUNE_THROW(Common::Exceptions::internal_error, "The inverse_ member is not filled after calling compute()!");
    const Common::Configuration default_opts =
        MatrixInverterOptions<MatrixType, true>::options(options_.get<std::string>("type"));
    if (options_.get("check_for_inf_nan", default_opts.get<bool>("check_for_inf_nan"))) {
      if (contains_inf_or_nan(*inverse_))
        DUNE_THROW(Exceptions::matrix_invert_failed_bc_result_contained_inf_or_nan,
                   "Computed inverse contains inf or nan and you requested checking. To disable this check set "
                   "'check_for_inf_nan' to false in the options."
                       << "\n\nThese were the given options:\n\n"
                       << options_ << "\nThis was the given matrix:\n\n"
                       << matrix_ << "\n\nThis was the computed inverse:\n\n"
                       << *inverse_);
    }
    const auto left_inverse_check =
        options_.get("post_check_is_left_inverse", default_opts.get<double>("post_check_is_left_inverse"));
    if (left_inverse_check > 0) {
      const auto eye = eye_matrix<MatrixType>(Common::get_matrix_cols(matrix_), Common::get_matrix_cols(matrix_));
      if (sup_norm(*inverse_ * matrix_ - eye) > left_inverse_check)
        DUNE_THROW(Exceptions::matrix_invert_failed_bc_result_is_not_a_left_inverse,
                   "Computed inverse is not a left inverse and you requested checking. To disable this check set "
                   "'post_check_is_left_inverse' to 0 in the options."
                       << "\n\nThe error is ||M_inv * M - Identity||_L_\\infty = "
                       << sup_norm(*inverse_ * matrix_ - eye) << "\n\nThese were the given options:\n\n"
                       << options_ << "\nThis was the given matrix M:\n\n"
                       << matrix_ << "\n\nThis is its computed inverse M_inv:\n\n"
                       << *inverse_ << "\n\nThis is M_inv * M\n\n"
                       << *inverse_ * matrix_);
    }
    const auto right_inverse_check =
        options_.get("post_check_is_right_inverse", default_opts.get<double>("post_check_is_right_inverse"));
    if (right_inverse_check > 0) {
      const auto eye = eye_matrix<MatrixType>(Common::get_matrix_rows(matrix_), Common::get_matrix_rows(matrix_));
      if (sup_norm(matrix_ * *inverse_ - eye) > right_inverse_check)
        DUNE_THROW(Exceptions::matrix_invert_failed_bc_result_is_not_a_right_inverse,
                   "Computed inverse is not a right inverse and you requested checking. To disable this check set "
                   "'post_check_is_right_inverse' to 0 in the options."
                       << "\n\nThe error is ||M_inv * M - Identity||_L_\\infty = "
                       << sup_norm(matrix_ * *inverse_ - eye) << "\n\nThese were the given options:\n\n"
                       << options_ << "\nThis was the given matrix M:\n\n"
                       << matrix_ << "\n\nThis is its computed inverse M_inv:\n\n"
                       << *inverse_ << "\n\nThis is M * M_inv\n\n"
                       << matrix_ * *inverse_);
    }
  } // ... post_checks(...)

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

  template <class M>
  double sup_norm(const MatrixInterface<M>& mat) const
  {
    return mat.sup_norm();
  }

  template <class M>
  typename std::enable_if<XT::Common::is_matrix<M>::value && !is_matrix<M>::value, double>::type
  sup_norm(const M& mat) const
  {
    using std::abs;
    using std::max;
    using Mat = XT::Common::MatrixAbstraction<M>;
    auto norm = abs(Mat::get_entry(mat, 0, 0));
    for (size_t ii = 0; ii < Mat::rows(mat); ++ii)
      for (size_t jj = 0; jj < Mat::cols(mat); ++jj)
        norm = std::max(norm, abs(Mat::get_entry(mat, ii, jj)));
    return norm;
  } // ... sup_norm(...)

  const MatrixType& matrix_;
  const Common::Configuration options_;
  mutable std::unique_ptr<MatrixType> inverse_;
}; // class MatrixInverterBase


} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_LA_MATRIX_INVERTER_BASE_HH
