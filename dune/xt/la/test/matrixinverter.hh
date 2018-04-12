// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2018)

#ifndef DUNE_XT_LA_TEST_MATRIXINVERTER_HH
#define DUNE_XT_LA_TEST_MATRIXINVERTER_HH

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/unused.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/xt/la/container/conversion.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/la/matrix-inverter.hh>

using namespace Dune;
using namespace Dune::XT;
using namespace Dune::XT::LA;


template <class MatrixImp>
struct MatrixInverterTest : public ::testing::Test
{
  using MatrixType = MatrixImp;
  typedef MatrixInverter<MatrixType> MatrixInverterType;
  typedef MatrixInverterOptions<MatrixType> MatrixInverterOpts;

  using M = Common::MatrixAbstraction<MatrixType>;

  MatrixInverterTest()
    : all_matrices_and_inverse_matrices_are_computed_(false)
    , broken_matrix_(M::create(M::has_static_size ? M::static_rows : 1, M::has_static_size ? M::static_cols : 1))
    , unit_matrix_(
          eye_matrix<MatrixType>(M::has_static_size ? M::static_rows : 1, M::has_static_size ? M::static_cols : 1))
  {
    M::set_entry(broken_matrix_, 0, 0, std::numeric_limits<typename M::S>::infinity());
  }

  static void exports_correct_types()
  {
    const bool MatrixType_is_correct = std::is_same<typename MatrixInverterType::MatrixType, MatrixType>::value;
    EXPECT_TRUE(MatrixType_is_correct);
  } // ... exports_correct_types()

  static void has_types_and_options()
  {
    std::vector<std::string> types = MatrixInverterOpts::types();
    EXPECT_GT(types.size(), 0);
    Common::Configuration opts = MatrixInverterOpts::options();
    EXPECT_TRUE(opts.has_key("type"));
    EXPECT_EQ(types[0], opts.get<std::string>("type"));
    for (const auto& tp : types) {
      Common::Configuration tp_opts = MatrixInverterOpts::options(tp);
      EXPECT_TRUE(tp_opts.has_key("type"));
      EXPECT_EQ(tp, tp_opts.get<std::string>("type"));
    }
  } // ... has_types_and_options(...)

  void throws_on_broken_matrix_construction() const
  {
    try {
      MatrixInverterType DXTC_UNUSED(default_inverter)(broken_matrix_);
      FAIL() << "Expected LA::Exceptions::eigen_inverter_failed_bc_data_did_not_fulfill_requirements";
    } catch (const LA::Exceptions::matrix_invert_failed_bc_data_did_not_fulfill_requirements& /*ee*/) {
    } catch (...) {
      FAIL() << "Expected LA::Exceptions::eigen_inverter_failed_bc_data_did_not_fulfill_requirements";
    }
    for (const auto& tp : MatrixInverterOpts::types()) {
      try {
        MatrixInverterType DXTC_UNUSED(default_opts_inverter)(broken_matrix_, tp);
        FAIL() << "Expected LA::Exceptions::eigen_inverter_failed_bc_data_did_not_fulfill_requirements";
      } catch (const LA::Exceptions::matrix_invert_failed_bc_data_did_not_fulfill_requirements& /*ee*/) {
      } catch (...) {
        FAIL() << "Expected LA::Exceptions::eigen_inverter_failed_bc_data_did_not_fulfill_requirements";
      }
      try {
        MatrixInverterType DXTC_UNUSED(inverter)(broken_matrix_, MatrixInverterOpts::options(tp));
        FAIL() << "Expected LA::Exceptions::eigen_inverter_failed_bc_data_did_not_fulfill_requirements";
      } catch (const LA::Exceptions::matrix_invert_failed_bc_data_did_not_fulfill_requirements& /*ee*/) {
      } catch (...) {
        FAIL() << "Expected LA::Exceptions::eigen_inverter_failed_bc_data_did_not_fulfill_requirements";
      }
    }
  } // ... throws_on_broken_matrix_construction(...)

  void is_constructible() const
  {
    ASSERT_TRUE(all_matrices_and_inverse_matrices_are_computed_);
    MatrixInverterType DXTC_UNUSED(default_inverter)(matrix_);
    for (const auto& tp : MatrixInverterOpts::types()) {
      MatrixInverterType DXTC_UNUSED(default_opts_inverter)(matrix_, tp);
      MatrixInverterType DXTC_UNUSED(inverter)(matrix_, MatrixInverterOpts::options(tp));
    }
  } // ... is_constructible(...)

  void gives_correct_inverse(const Common::Configuration& tolerances = {}) const
  {
    ASSERT_TRUE(all_matrices_and_inverse_matrices_are_computed_);
    for (const auto& tp : MatrixInverterOpts::types()) {
      const double tolerance = tolerances.get(tp, 1e-15);
      try {
        MatrixInverterType inverter(matrix_, tp);
        const auto& actual_inverse = inverter.inverse();
        EXPECT_EQ(Common::get_matrix_rows(matrix_), Common::get_matrix_rows(actual_inverse));
        EXPECT_EQ(Common::get_matrix_rows(matrix_), Common::get_matrix_cols(actual_inverse));
        if (tolerance > 0) {
          EXPECT_TRUE(Common::FloatCmp::eq(actual_inverse, expected_inverse_, tolerance))
              << "\n\nactual inverse: " << actual_inverse << "\n\nexpected inverse: " << expected_inverse_
              << "\n\ntolerance: " << tolerance << "\n\ntype: " << tp;
        } else {
          // negative tolerance: we expect a failure
          EXPECT_FALSE(Common::FloatCmp::eq(actual_inverse, expected_inverse_))
              << "\n\nTHIS IS A GOOD THING! UPDATE THE EXPECTATIONS IN tolerances!\n\n"
              << "\n\nactual inverse: " << actual_inverse << "\n\nexpected inverse: " << expected_inverse_
              << "\n\ntype: " << tp;
        }
      } catch (const Dune::MathError&) {
        if (tolerance > 0) {
          FAIL() << "Dune::MathError thrown when trying to get inverse!"
                 << "\n\ntype: " << tp << "\n\ntolerance: " << tolerance;
        }
      }
    }
  } // ... gives_correct_eigenvectors_in_correct_order(...)

  bool all_matrices_and_inverse_matrices_are_computed_;
  MatrixType broken_matrix_;
  MatrixType unit_matrix_;
  MatrixType matrix_;
  MatrixType expected_inverse_;
}; // ... struct MatrixInverterTest


#endif // DUNE_XT_LA_TEST_MATRIXINVERTER_HH
