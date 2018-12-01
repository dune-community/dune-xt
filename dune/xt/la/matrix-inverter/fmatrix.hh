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

#ifndef DUNE_XT_LA_MATRIX_INVERTER_FMATRIX_HH
#define DUNE_XT_LA_MATRIX_INVERTER_FMATRIX_HH

#include <dune/common/typetraits.hh>

#include <dune/xt/common/fmatrix.hh>

#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/container/conversion.hh>
#include <dune/xt/la/container/eigen/dense.hh>
#include <dune/xt/la/matrix-inverter.hh>

#include "internal/base.hh"
#include "internal/eigen.hh"

namespace Dune {
namespace XT {
namespace LA {


template <class K, int ROWS, int COLS>
class MatrixInverterOptions<FieldMatrix<K, ROWS, COLS>, true>
{
public:
  static std::vector<std::string> types()
  {
    return
    {
      "direct"
#if HAVE_EIGEN
          ,
          "moore_penrose"
#endif
    };
  }

  static Common::Configuration options(const std::string type = "")
  {
    const std::string actual_type = type.empty() ? types()[0] : type;
    internal::ensure_matrix_inverter_type(actual_type, types());
    Common::Configuration opts = internal::default_matrix_inverter_options();
    opts["type"] = actual_type;
    return opts;
  }
}; // class MatrixInverterOptions<FieldMatrix<K, ROWS, COLS>>

template <class K, int ROWS, int COLS>
class MatrixInverterOptions<Common::FieldMatrix<K, ROWS, COLS>, true>
  : public MatrixInverterOptions<FieldMatrix<K, ROWS, COLS>>
{};

template <class K, int ROWS, int COLS>
class MatrixInverter<FieldMatrix<K, ROWS, COLS>, true> : public internal::MatrixInverterBase<FieldMatrix<K, ROWS, COLS>>
{
  using BaseType = internal::MatrixInverterBase<FieldMatrix<K, ROWS, COLS>>;

public:
  using MatrixType = typename BaseType::MatrixType;

  template <class... Args>
  explicit MatrixInverter(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
    const auto type = options_.template get<std::string>("type");
    const XT::Common::Configuration default_opts = MatrixInverterOptions<MatrixType>::options(type);
    if (!options_.get("delay_computation", default_opts.get<bool>("delay_computation")))
      compute();
  }

  void compute() override final
  {
    const auto type = options_.template get<std::string>("type");
    if (type == "direct") {
      inverse_ = std::make_unique<MatrixType>(matrix_);
      auto inverse_xt = std::make_unique<XT::Common::FieldMatrix<K, ROWS, COLS>>(*inverse_);
      try {
        inverse_xt->invert();
        *inverse_ = *inverse_xt;
      } catch (const FMatrixError& ee) {
        if (std::strcmp(ee.what(), "matrix is singular") != 0)
          DUNE_THROW(Exceptions::matrix_invert_failed_bc_data_did_not_fulfill_requirements,
                     "This was the original error:\n\n"
                         << ee.what());
        DUNE_THROW(Exceptions::matrix_invert_failed, "This was the original error:\n\n" << ee.what());
      }
    }
#if HAVE_EIGEN
    else if (type == "moore_penrose") {
      inverse_ = std::make_unique<MatrixType>(convert_to<MatrixType>(EigenDenseMatrix<K>(
          internal::compute_moore_penrose_inverse_using_eigen(convert_to<EigenDenseMatrix<K>>(matrix_).backend()))));
    }
#endif
    else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type
                                << "' is none of MatrixInverterOptions<FieldMatrix<K, ROWS, COLS>>::types(), and "
                                   "internal::MatrixInverterBase promised to check this!"
                                << "\n\nThese are the available types:\n\n"
                                << MatrixInverterOptions<MatrixType>::types());

    this->post_checks();
  } // ... compute(...)

protected:
  using BaseType::inverse_;
  using BaseType::matrix_;
  using BaseType::options_;
}; // class MatrixInverter<FieldMatrix<...>>


template <class K, int ROWS, int COLS>
class MatrixInverter<Common::FieldMatrix<K, ROWS, COLS>, true> : public MatrixInverter<FieldMatrix<K, ROWS, COLS>>
{
public:
  template <class... Args>
  explicit MatrixInverter(Args&&... args)
    : MatrixInverter<FieldMatrix<K, ROWS, COLS>>(std::forward<Args>(args)...)
  {}
};


} // namespace Dune
} // namespace XT
} // namespace LA


#endif // DUNE_XT_LA_MATRIX_INVERTER_FMATRIX_HH
